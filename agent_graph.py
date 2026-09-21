import sys
sys.dont_write_bytecode = True
import os
from dotenv import load_dotenv
from langchain_core.tools import tool
from langchain_google_genai import ChatGoogleGenerativeAI
from langgraph.prebuilt import create_react_agent

from kegg_client import (
    search_kegg_pathway,
    get_kegg_pathway_description,
    get_kegg_pathway_details,
    get_kegg_organism_code,
    search_kegg_disease,
    get_disease_linked_genes_and_pathways,
    get_kegg_pathways,
    get_kegg_diseases,
    convert_ncbi_to_kegg_id,
    translate_kegg_genes_bulk,
    search_kegg_compound,
    get_kegg_compound_details,
    get_reactome_pathway_description,
)
from ncbi_client import get_ncbi_gene_id, get_ncbi_gene_summary
from clinvar_client import fetch_top_mutations
from pubtator_client import search_pubtator_papers, fetch_pubtator_annotations
from wiki_client import fetch_wikipedia_summary
from reactome_client import get_reactome_pathway_data

load_dotenv()


# ---- GLOBAL PROCESSING THRESHOLDS ---- #
# These control how much data is fetched and returned to the AI at once.
# You can safely tweak these numbers here without touching any other code.
MAX_PATHWAY_MATCHES = 10     # Max number of pathways to deep-dive per query (to avoid IP bans)
MAX_GENES_RETURNED = 500     # Max total genes to return across all matched pathways
MAX_COMPOUNDS_RETURNED = 500 # Max total compounds to return across all matched pathways
MAX_DISEASES_RETURNED = 100  # Max total diseases to return



# ---- PAGINATION HELPER ---- #

def paginate_results(data_list, page, chunk_size=200):
    # data_list: The full list of results to paginate.
    # page: The current page number (1-indexed).
    # chunk_size: How many items to return per page.
    start = (page - 1) * chunk_size
    end = start + chunk_size
    chunk = data_list[start:end]
    has_more = end < len(data_list)
    return {
        "data": chunk,
        "current_page": page,
        "has_more": has_more,
        "next_page_to_call": page + 1 if has_more else None,
        "total_items": len(data_list)
    }


# ---- TOOL DEFINITIONS ---- #
# Each function below is decorated with @tool so the ReAct agent LLM can
# decide when to call it on its own, based on the user's query.

@tool
def tool_search_kegg_pathway(pathway_name: str, fields_needed: list = None, page: int = 1, database_preference: str = "both"):
    """
    Searches for a biological pathway by name.
    Returns ONLY the fields specified in fields_needed to keep the response focused.
    Available fields: 'genes', 'compounds', 'drugs', 'modules', 'related_pathways', 'description'.
    Genes and compounds are fetched from both KEGG and Reactome simultaneously and merged.
    IMPORTANT: The returned genes and compounds have [Source: KEGG] and/or [Source: Reactome] tags. You MUST preserve and display these exact tags to the user so they know where the data came from!
    Use this when the user asks about a pathway (e.g. glycolysis, TCA cycle, etc.)
    Always specify only the fields_needed that are relevant to the user's question.
    database_preference: Optional. Can be 'both', 'kegg', or 'reactome'. If the user asks for data from a specific database, restrict the search using this.
    If the result contains 'has_more: True', you MUST call this tool again with page = next_page_to_call.
    """
    # pathway_name: The pathway name to search for (e.g. "glycolysis").
    # fields_needed: A list of strings specifying which data fields to return.
    #                Example: ["genes", "drugs"] if user asks about genes and drugs only.
    # page: The page number for paginated gene/compound lists (default is 1).
    print(f"\n[Tool: search_kegg_pathway] Searching for: '{pathway_name}', Fields: {fields_needed}, DB: {database_preference}, Page: {page}")
    
    # Helper for smart deduplication and tagging
    def merge_and_tag(kegg_list, reactome_list):
        merged = {}
        for item in kegg_list:
            base_name = item.split(" (ID:")[0].split(" [")[0].strip().lower()
            merged[base_name] = {"display": item, "tags": {"KEGG"}}
        for item in reactome_list:
            base_name = item.split(" (ID:")[0].split(" [")[0].strip().lower()
            if base_name in merged:
                merged[base_name]["tags"].add("Reactome")
            else:
                merged[base_name] = {"display": item, "tags": {"Reactome"}}
                
        result = []
        for data in merged.values():
            tags = list(data["tags"])
            tags.sort()
            tag_str = ", ".join(tags)
            result.append(f"{data['display']} [Source: {tag_str}]")
        return sorted(result)

    if not fields_needed:
        fields_needed = []

    matches = search_kegg_pathway(pathway_name)
    has_kegg = len(matches) > 0

    if not has_kegg:
        # KEGG found nothing at all, try Reactome as a last resort
        reactome_data = get_reactome_pathway_data(pathway_name)
        if not reactome_data or not reactome_data.get("reactome_id"):
            return {"error": f"No pathway found for '{pathway_name}' in KEGG or Reactome."}
        return {
            "pathway_name_searched": pathway_name,
            "reactome_pathway_id": reactome_data.get("reactome_id"),
            "genes": paginate_results(reactome_data.get("genes", []), page),
            "compounds": paginate_results(reactome_data.get("compounds", []), page),
        }

    # List all matched pathway names for the AI to see
    result = {
        "all_matched_pathways": [m["name"] for m in matches],
        "total_pathways_found": len(matches),
        "pathways_deep_dived": min(len(matches), MAX_PATHWAY_MATCHES),
    }

    # Aggregate fields across all matches up to the MAX_PATHWAY_MATCHES threshold
    org_code = get_kegg_organism_code("human")
    needs_details = any(f in fields_needed for f in ["genes", "compounds", "drugs", "modules", "related_pathways"])
    needs_description = "description" in fields_needed and page == 1

    combined_kegg_genes = []
    combined_kegg_compounds = []
    combined_reactome_genes = []
    combined_reactome_compounds = []
    
    combined_drugs = []
    combined_modules = []
    combined_rel_pathways = []
    descriptions = []

    # Loop through matches up to the safety threshold
    for i, match in enumerate(matches[:MAX_PATHWAY_MATCHES]):
        pid = match["id"]
        pname = match.get("name", "")
        
        # Only print aggregation logs on the first page to avoid terminal spam
        if page == 1:
            print(f"  Deep-diving pathway {i+1}/{min(len(matches), MAX_PATHWAY_MATCHES)}: {pname}")

        if needs_details:
            if database_preference in ["both", "kegg"]:
                kegg_details = get_kegg_pathway_details(pid, org_code)
                combined_kegg_genes.extend(kegg_details.get("genes", []))
                combined_kegg_compounds.extend(kegg_details.get("compounds", []))
                combined_drugs.extend(kegg_details.get("drugs", []))
                combined_modules.extend(kegg_details.get("modules", []))
                combined_rel_pathways.extend(kegg_details.get("rel_pathways", []))
                
            if database_preference in ["both", "reactome"]:
                reactome_data = get_reactome_pathway_data(pname)
                combined_reactome_genes.extend(reactome_data.get("genes", []))
                combined_reactome_compounds.extend(reactome_data.get("compounds", []))

        if needs_description:
            desc_text = ""
            
            # Fetch KEGG if allowed
            if database_preference in ["both", "kegg"]:
                kegg_desc = get_kegg_pathway_description(pid)
                if kegg_desc:
                    desc_text = f"{kegg_desc} [Source: KEGG]"
            
            # Fetch Reactome if allowed and no description yet
            if not desc_text and database_preference in ["both", "reactome"]:
                reactome_desc = get_reactome_pathway_description(pname)
                if reactome_desc:
                    desc_text = f"{reactome_desc} [Source: Reactome]"
            
            # Fallback to Wikipedia if still nothing
            if not desc_text:
                wiki_desc = fetch_wikipedia_summary(pname)
                if wiki_desc:
                    desc_text = f"{wiki_desc} [Source: Wikipedia]"
            
            if desc_text:
                descriptions.append(f"[{pname}]: {desc_text}")

    # Merge and Tag (Deduplicate)
    combined_genes = merge_and_tag(combined_kegg_genes, combined_reactome_genes)
    combined_compounds = merge_and_tag(combined_kegg_compounds, combined_reactome_compounds)
    
    combined_drugs = list(dict.fromkeys(combined_drugs))
    combined_modules = list(dict.fromkeys(combined_modules))
    combined_rel_pathways = list(dict.fromkeys(combined_rel_pathways))

    # Apply global thresholds to enforce memory safety before pagination
    combined_genes = combined_genes[:MAX_GENES_RETURNED]
    combined_compounds = combined_compounds[:MAX_COMPOUNDS_RETURNED]

    if page == 1:
        print(f"  Aggregated: {len(combined_genes)} genes, {len(combined_compounds)} compounds, {len(combined_drugs)} drugs.")

    # Paginate the massive combined lists
    if "genes" in fields_needed:
        result["genes"] = paginate_results(combined_genes, page)

    if "compounds" in fields_needed:
        result["compounds"] = paginate_results(combined_compounds, page)

    # Drugs, modules, related_pathways are only returned on the first page to save memory
    if page == 1:
        if "drugs" in fields_needed:
            result["drugs"] = combined_drugs
        if "modules" in fields_needed:
            result["modules"] = combined_modules
        if "related_pathways" in fields_needed:
            result["related_pathways"] = combined_rel_pathways
        if needs_description and descriptions:
            result["descriptions"] = descriptions

    return result


@tool
def tool_search_compound(compound_name: str):
    """
    Searches KEGG for a specific chemical compound or metabolite (e.g. 'glucose', 'ATP', 'NADH').
    Returns the exact chemical formula, mass, related biological reactions, and pathways it participates in.
    Use this when the user asks for elaboration on a specific chemical compound, molecule, or metabolite.
    """
    # compound_name: The name of the chemical compound (e.g. "glucose").
    print(f"\n[Tool: search_compound] Searching for: '{compound_name}'")
    matches = search_kegg_compound(compound_name)
    if not matches:
        return {"error": f"No compound found for '{compound_name}' in KEGG."}

    # Grab the top match (most relevant)
    top_match = matches[0]
    compound_id = top_match["id"]
    
    details = get_kegg_compound_details(compound_id)
    return details


@tool
def tool_search_kegg_disease(disease_name: str, page: int = 1):
    """
    Searches the KEGG disease database for a disease by name.
    Returns the disease ID, name, a Wikipedia overview, and the list of
    linked genes and pathways (for humans).
    Use this when the user asks about a disease and wants its linked genes or pathways.
    If the result contains 'has_more: True', you MUST call this tool again with page = next_page_to_call.
    """
    # disease_name: The disease name to search for (e.g. "Alzheimer disease").
    # page: The page number for paginated results (default is 1).
    print(f"\n[Tool: search_kegg_disease] Searching for: '{disease_name}', Page: {page}")
    matches = search_kegg_disease(disease_name)
    if not matches:
        return {"error": f"No KEGG disease found for '{disease_name}'."}

    # Deduplicate matches
    unique_matches = {d["id"]: d for d in matches}.values()
    
    # Take up to the global safety limit
    safe_matches = list(unique_matches)[:MAX_DISEASES_RETURNED]
    
    results = []
    # If the AI asks for page 1, fetch 1-20, page 2 fetches 21-40, etc.
    # To save API calls, we ONLY process the chunk that is being paginated.
    chunk_size = 5 # Process 5 diseases per page so API isn't overloaded
    start = (page - 1) * chunk_size
    end = start + chunk_size
    
    current_chunk = safe_matches[start:end]
    has_more = end < len(safe_matches)
    
    for d_match in current_chunk:
        d_id = d_match["id"]
        wiki = fetch_wikipedia_summary(d_match["name"])
        linked = get_disease_linked_genes_and_pathways(d_id, "human")
        
        # Don't truncate genes anymore! Give them all!
        raw_gene_ids = linked.get("gene_ids", [])
        translated_genes = translate_kegg_genes_bulk(raw_gene_ids)
        
        results.append({
            "disease_id": d_id,
            "disease_name": d_match["name"],
            "overview": wiki,
            "linked_genes": translated_genes,
            "linked_pathway_ids": linked.get("pathway_ids", []),
        })
        
    return {
        "total_diseases_found": len(safe_matches),
        "data": results,
        "current_page": page,
        "has_more": has_more,
        "next_page_to_call": page + 1 if has_more else None
    }


@tool
def tool_get_ncbi_gene_info(gene_symbol: str):
    """
    Fetches the NCBI gene ID, full description, official gene summary,
    organism, and all KEGG pathways and diseases linked to the given gene symbol.
    Use this when the user asks for information about a specific gene.
    """
    # gene_symbol: The official gene symbol (e.g. "TP53", "BRCA1").
    print(f"\n[Tool: get_ncbi_gene_info] Fetching info for gene: '{gene_symbol}'")
    ncbi_id = get_ncbi_gene_id(gene_symbol, "human")
    if not ncbi_id:
        return {"error": f"No NCBI Gene ID found for '{gene_symbol}'."}

    summary = get_ncbi_gene_summary(ncbi_id)
    kegg_id = convert_ncbi_to_kegg_id(ncbi_id)

    pathways = []
    diseases = []
    if kegg_id:
        pathways = get_kegg_pathways(kegg_id)
        diseases = get_kegg_diseases(kegg_id)

    return {
        "ncbi_gene_id": ncbi_id,
        "kegg_gene_id": kegg_id,
        "symbol": summary.get("symbol", ""),
        "description": summary.get("description", ""),
        "summary": summary.get("summary", ""),
        "organism": summary.get("organism", ""),
        "kegg_pathways": pathways[:10],
        "kegg_diseases": diseases[:10],
    }


@tool
def tool_fetch_clinvar_mutations(gene_symbols: list):
    """
    Fetches known pathogenic mutations for a list of human genes from ClinVar.
    The returned data includes the mutation Name, RS_ID, Type, Phenotypes (diseases), Origin, and Significance.
    Use this when the user asks about mutations, variants, or genetic changes in genes.
    """
    # gene_symbols: A list of official human gene symbols (e.g. ["BRCA1", "TP53"]).
    print(f"\n[Tool: fetch_clinvar_mutations] Fetching mutations for: {gene_symbols}")
    mutations = fetch_top_mutations(gene_symbols)
    if not mutations:
        return {"result": f"No pathogenic mutations found in ClinVar for {gene_symbols}."}
    return {"mutations_found": mutations}


@tool
def tool_search_pubmed(search_query: str, page: int = 1):
    """
    Searches PubMed (via PubTator3) for relevant research papers on a given topic.
    Returns paper titles, authors, journal names, publication years, PMIDs, and advanced NLP annotations 
    including Genes, Chemicals, Diseases, Species, Cell Lines, Variants/Mutations, and structural Relations.
    Use this when the user asks for research papers, publications, or literature.
    If the result contains 'has_more: True', you MUST call this tool again with page = next_page_to_call.
    """
    # search_query: The biomedical topic or terms to search for (e.g. "TP53 cancer mutations").
    # page: The page number for paginated results (default is 1).
    print(f"\n[Tool: search_pubmed] Searching for: '{search_query}', Page: {page}")
    papers = search_pubtator_papers(search_query)
    if not papers:
        return {"result": "No papers found for this query."}

    results = []
    for paper in papers:
        pmid = paper.get("pmid")
        annotations = fetch_pubtator_annotations(pmid)
        results.append({
            "title": paper.get("title"),
            "authors": paper.get("authors"),
            "journal": paper.get("journal"),
            "year": paper.get("year"),
            "pmid": pmid,
            "genes_mentioned": annotations.get("Genes", []),
            "diseases_mentioned": annotations.get("Diseases", []),
            "chemicals_mentioned": annotations.get("Chemicals", []),
            "species_mentioned": annotations.get("Species", []),
            "cell_lines_mentioned": annotations.get("CellLines", []),
            "variants_mentioned": annotations.get("Variants", []),
            "relations": annotations.get("Relations", []),
        })
    return paginate_results(results, page)


# ---- AGENT BUILD FUNCTION ---- #

def build_ncbi_kegg_graph():
    # Builds and returns the dynamic ReAct tool-calling agent.
    # The LLM will decide which tools to call and in what order based on the user query.

    api = os.getenv("LLM_API")
    llm = ChatGoogleGenerativeAI(
        model="gemini-flash-lite-latest",
        temperature=0.2,
        api_key=api
    )

    tools = [
        tool_search_kegg_pathway,
        tool_search_kegg_disease,
        tool_get_ncbi_gene_info,
        tool_fetch_clinvar_mutations,
        tool_search_pubmed,
        tool_search_compound,
    ]

    system_prompt = (
        "You are an expert biomedical research assistant with access to NCBI, KEGG, ClinVar, and PubMed databases via tools. "
        "CRITICAL RULE 1: If the user asks a question that is CLEARLY NOT related to biology, medicine, genetics, pathways, or diseases "
        "(for example, asking about actors, movies, general knowledge, or programming), you MUST refuse to answer. "
        "Reply exactly with: 'I am a specialized biomedical agent. I can only answer questions related to biology, genetics, diseases, and medical research.' "
        "Do NOT answer the off-topic question under any circumstances. "
        "HOWEVER, if the query contains ANY biological terms, abbreviations, or references (like 'glycolysis', 'CO IDs', 'kegg', 'TP53', 'Human Kegg ID'), you MUST accept it as a valid biomedical query. Do NOT falsely reject it. "
        "CRITICAL RULE 2: Do NOT use LaTeX math formatting under any circumstances (e.g., absolutely no $\\text{CO}_2$ or $\\text{NADH}$). Use plain text like CO2 and NADH. "
        "When the user asks a valid biomedical question, analyze it carefully and call ONLY the tools that are necessary to answer it. "
        "CRITICAL RULE 3: You MUST base your answers SOLELY on the information returned by your tools. Do NOT answer questions using your pre-trained internal knowledge. If your tools do not return relevant information, you must reply that you cannot find the answer in the databases. "
        "CRITICAL RULE 4 (PAGINATION): Tools return data in chunks (pages). When a tool result contains 'has_more: True', "
        "you are STRICTLY FORBIDDEN from making another tool call until you have output the Markdown text summarizing the current chunk! "
        "You MUST first write out the data from the current chunk as a clean Markdown section. "
        "Only AFTER you have printed the text for the current chunk are you allowed to call the SAME tool again with the 'next_page_to_call' number. "
        "Keep looping and writing chunk by chunk until 'has_more' is False. Do NOT queue up multiple page calls at once. "
        "Do NOT call tools that are unrelated to the question. "
        "If the user asks for compounds in a pathway, use the pathway tool only. "
        "If the user asks for mutations of genes found in a pathway, first call the pathway tool to get the genes, "
        "then call the ClinVar tool for each relevant gene. "
        "After gathering the data, synthesize a clean, structured Markdown answer. "
        "Do NOT include introductory or concluding filler sentences. Jump straight into the facts. "
        "Do NOT add a 'Suggestions for Further Reading' section unless the user specifically asks for papers. "
        "Focus strictly on what the user asked. If they asked only for compounds, list only compounds. "
        "If they asked only for mutations, list only mutations. "
        "CRITICAL RULE 5: Do NOT arbitrarily truncate or summarize lists! If a tool returns a specific number of items (whether it is mutations, genes, compounds, pathways, or literature annotations), you MUST list EVERY SINGLE ITEM in your Markdown response. Do not drop items."
    )

    agent = create_react_agent(llm, tools, prompt=system_prompt)
    return agent


# Detailed Script Workings:
# 1. Each tool_ function wraps the underlying API client functions from kegg_client.py,
#    ncbi_client.py, clinvar_client.py, and pubtator_client.py using the @tool decorator.
# 2. The @tool decorator gives each function a clear docstring that the LLM reads to
#    decide which tool is relevant to the user's query.
# 3. build_ncbi_kegg_graph() initializes the Gemini LLM, binds all tools to it,
#    and creates a ReAct agent using langgraph.prebuilt.create_react_agent.
# 4. The agent runs a Reasoning -> Acting -> Observation loop automatically:
#    it reads the query, calls tools, reads their output, and decides what to do next.
# 5. No hardcoded edges, static nodes, or manual state routing are needed anymore.
# 6. gene_query.py calls build_ncbi_kegg_graph() and invokes the agent the same way as before.
