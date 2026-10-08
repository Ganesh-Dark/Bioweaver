import sys
sys.dont_write_bytecode = True
import os
from dotenv import load_dotenv
from langchain_core.tools import tool
from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_groq import ChatGroq
from langchain_huggingface import ChatHuggingFace, HuggingFaceEndpoint
from langgraph.prebuilt import create_react_agent

from clients.kegg_client import (
    search_kegg_pathway,
    get_kegg_pathway_description,
    get_kegg_disease_description,
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
    search_kegg_ko,
    get_kegg_ko_details,
    search_kegg_enzyme,
    get_kegg_enzyme_details,
    search_kegg_module,
    get_kegg_module_details,
    search_kegg_reaction,
    get_kegg_reaction_details,
    get_pathway_linked_data,
)

from clients.ncbi_client import get_ncbi_gene_id, get_ncbi_gene_summary
from clients.clinvar_client import fetch_top_mutations
from clients.hf_router_client import route_clinvar_query
from clients.wiki_client import fetch_wikipedia_summary
from clients.reactome_client import get_reactome_pathway_data
from analytics import kegg_analytics as pa
from analytics import reactome_analytics as ra
from analytics import compound_resolver as cr

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
        "total_items": len(data_list)
    }


# ---- TOOL DEFINITIONS ---- #
# Each function below is decorated with @tool so the ReAct agent LLM can
# decide when to call it on its own, based on the user's query.

@tool
def tool_search_kegg_pathway(pathway_name: str, fields_needed: list = None, page: int = 1, database_preference: str = "kegg", just_list: bool = False):
    """
    Searches biological pathway by name.
    CRITICAL: If the user only asks for a 'list' of pathways, you MUST set just_list=True.
    Returns ONLY the fields specified in fields_needed.
    Available fields: 'genes', 'compounds', 'drugs', 'modules', 'related_pathways', 'description'.
    IMPORTANT: By default, this ONLY searches KEGG to ensure 100% accurate results without hallucinated fuzzy matches. 
    DO NOT set database_preference to 'both' or 'reactome' UNLESS the user EXPLICITLY asks for Reactome related data!
    Always specify only the fields_needed that are relevant to the user's question. If the user asks for general information (e.g., 'tell me about X'), ONLY fetch the 'description' field. Do not fetch 'genes' or 'compounds' unless explicitly requested.
    CRITICAL: If the 'description' field is fetched, you MUST output the exact description provided by the tool word-for-word. Do NOT paraphrase the description on your own.
    If the result contains 'has_more: True', DO NOT automatically fetch the next page unless the user explicitly asks for more results.
    """
    # pathway_name: The pathway name to search for (e.g. "glycolysis").
    # fields_needed: A list of strings specifying which data fields to return.
    #                Example: ["genes", "drugs"] if user asks about genes and drugs only.
    # page: The page number for paginated gene/compound lists (default is 1).
    print(f"\n[Tool: search_kegg_pathway] Searching for: '{pathway_name}', Fields: {fields_needed}, DB: {database_preference}, Page: {page}, Just List: {just_list}")

    if just_list:
        fields_needed = []
        
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

    if database_preference == "reactome":
        matches = []
    else:
        matches = search_kegg_pathway(pathway_name)

    has_kegg = len(matches) > 0

    if just_list:
        result = {
            "all_matched_pathways": [{"id": m["id"], "name": m["name"]} for m in matches],
            "total_pathways_found": len(matches),
        }
        if database_preference in ["both", "reactome"]:
            from clients.reactome_client import get_reactome_pathways_list
            reactome_matches = get_reactome_pathways_list(pathway_name)
            if reactome_matches:
                result["reactome_matched_pathways"] = reactome_matches
                result["total_pathways_found"] += len(reactome_matches)
        return result

    if database_preference == "reactome" or (not has_kegg and database_preference in ["both"]):
        # KEGG found nothing at all, try Reactome as a last resort
        reactome_data = get_reactome_pathway_data(pathway_name)
        if isinstance(reactome_data, dict) and "ambiguous_options" in reactome_data:
            return {"error": f"Reactome returned multiple possible matches for '{pathway_name}'. Please present these options to the user clearly, and ask them to reply with the exact Name or ID of the pathway they want:", "options": reactome_data["ambiguous_options"]}
            
        if not reactome_data or not reactome_data.get("reactome_id"):
            return {"error": f"No pathway found for '{pathway_name}' in KEGG or Reactome."}
            
        result = {
            "pathway_name_searched": pathway_name,
            "reactome_pathway_id": reactome_data.get("reactome_id"),
            "genes": paginate_results(reactome_data.get("genes", []), page),
            "compounds": paginate_results(reactome_data.get("compounds", []), page),
        }
        
        if "description" in fields_needed and page == 1:
            reactome_desc = get_reactome_pathway_description(pathway_name)
            if reactome_desc and not isinstance(reactome_desc, dict):
                result["description"] = f"{reactome_desc} [Source: Reactome]"
                
        return result

    # List all matched pathway IDs and names for the AI to see
    # Both ID and name are returned so the AI can use the ID in tool_get_pathway_linked_data
    result = {
        "all_matched_pathways": [{"id": m["id"], "name": m["name"]} for m in matches],
        "total_pathways_found": len(matches),
    }
        
    result["pathways_deep_dived"] = min(len(matches), MAX_PATHWAY_MATCHES)

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
                if isinstance(reactome_data, dict) and "ambiguous_options" in reactome_data:
                    return {"error": f"Reactome returned multiple possible matches for '{pname}'. Please present these options to the user clearly, and ask them to reply with the exact Name or ID of the pathway they want:", "options": reactome_data["ambiguous_options"]}
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
                
                if isinstance(reactome_desc, dict) and "ambiguous_options" in reactome_desc:
                    return {"error": f"Reactome returned multiple possible matches for '{pname}'. Please present these options to the user clearly, and ask them to reply with the exact Name or ID of the pathway they want:", "options": reactome_desc["ambiguous_options"]}
                    
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
def tool_search_compound(compound_name: str, include_pathways: bool = False, include_reactions: bool = False, include_modules: bool = False, include_enzymes: bool = False):
    """
    Searches KEGG for a specific chemical compound or metabolite (e.g. 'glucose', 'ATP').
    Returns the exact chemical formula, mass, and alternative names.
    By default, it DOES NOT return associated pathways, modules, enzymes, or reactions to save space.
    ONLY set the include_* flags to True if the user EXPLICITLY asks for them (e.g. "what pathways use glucose?").
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
    
    if not include_pathways:
        details.pop("pathways", None)
    if not include_reactions:
        details.pop("reactions", None)
    if not include_modules:
        details.pop("modules", None)
    if not include_enzymes:
        details.pop("enzymes", None)
        
    return details


@tool
def tool_search_kegg_disease(disease_name: str, page: int = 1):
    """
    Searches the KEGG disease database for a disease by name.
    Returns the disease ID, name, a Wikipedia overview, and the list of
    linked genes and pathways (for humans).
    Use this when the user asks about a disease and wants its linked genes or pathways.
    If the result contains 'has_more: True', DO NOT automatically call this tool again unless the user specifically asked you to fetch ALL results.
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
    # Hard cap: only process the first 5 results per call so the AI cannot loop endlessly
    chunk_size = 5
    start = (page - 1) * chunk_size
    end = start + chunk_size
    current_chunk = safe_matches[start:end]
    
    results = []
    for d_match in current_chunk:
        d_id = d_match["id"]
        
        kegg_desc = get_kegg_disease_description(d_id)
        if kegg_desc:
            overview_text = f"{kegg_desc} [Source: KEGG]"
        else:
            wiki = fetch_wikipedia_summary(d_match["name"])
            overview_text = f"{wiki} [Source: Wikipedia (Fallback used because KEGG lacked a description)]" if wiki else ""
            
        linked = get_disease_linked_genes_and_pathways(d_id, "human")
        
        raw_gene_ids = linked.get("gene_ids", [])
        translated_genes = translate_kegg_genes_bulk(raw_gene_ids)
        
        results.append({
            "disease_id": d_id,
            "disease_name": d_match["name"],
            "overview": overview_text,
            "linked_genes": translated_genes,
            "linked_pathway_ids": linked.get("pathway_ids", []),
        })
        
    return {
        "total_diseases_found": len(safe_matches),
        "showing_top": len(results),
        "data": results,
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
def tool_fetch_clinvar_mutations(raw_user_query="", default_targets=None):
    """
    Fetches known pathogenic mutations for specific genes or RS IDs from ClinVar.
    The returned data includes the mutation Name, RS_ID, Type, Phenotypes (diseases), Origin, and Significance.
    Use this when the user asks about mutations, variants, or genetic changes.
    
    Parameters:
    raw_user_query: The exact phrasing the user used (e.g., "What are the copy number losses for TP53?").
    default_targets: A list of core gene symbols or RS IDs to use as a fallback if the router fails to extract them (e.g., ["LRRK2", "SNCA"] or ["121918399"]).
    """
    if default_targets is None:
        default_targets = []
        
    print(f"\n[Tool: fetch_clinvar_mutations] Routing query: '{raw_user_query}'")
    hf_config = route_clinvar_query(raw_user_query, default_targets=default_targets)
    
    targets = hf_config.get("targets")
    if not targets and default_targets:
        targets = default_targets
        
    if not targets:
        return {"error": "Could not identify target genes or RS IDs from the query."}
        
    target_type = hf_config.get("target_type", "gene")
    mutation_type = hf_config.get("mutation_type", "single nucleotide variant")
    fallback = hf_config.get("fallback_allowed", True)
    
    mutations = fetch_top_mutations(
        targets,
        target_type=target_type,
        mutation_type=mutation_type,
        fallback_allowed=fallback,
        limit=5
    )
    
    if not mutations:
        return {"result": f"No pathogenic mutations found in ClinVar for {targets}."}
    return {"mutations_found": mutations}





@tool
def tool_search_ko(query: str, fetch_details: bool = False, filter_field: str = "", filter_value: str = ""):
    """
    Searches for a KEGG Orthology (KO) term.
    Set fetch_details=True ONLY if query is a specific KO ID (e.g. 'K00001') to get full details.
    IMPORTANT: Even if full details are fetched, ONLY output basic info (ID, Name, Class) to the user unless they explicitly asked for genes, pathways, or other massive fields.
    Use filter_field (e.g. 'pathway') and filter_value to strictly filter the initial matches.
    """
    # query: The search term
    # fetch_details: If True, fetches full details. Use ONLY for a specific ID.
    # filter_field: Optional field name to filter by (e.g. 'pathway')
    # filter_value: Optional value that must be present in the filter_field
    print(f"\n[Tool: search_ko] Searching for: '{query}'")
    matches = search_kegg_ko(query, filter_field, filter_value)
    if not matches:
        return {"error": f"No KO found for '{query}'"}
    
    is_id = query.upper().startswith("K") and len(query) > 1 and query[1].isdigit()
    if fetch_details or is_id:
        return get_kegg_ko_details(matches[0]["id"])
    return {"matches": matches[:10]}

@tool
def tool_search_enzyme(query: str, fetch_details: bool = False, filter_field: str = "", filter_value: str = ""):
    """
    Searches for a KEGG Enzyme.
    Set fetch_details=True ONLY if query is a specific EC number (e.g. '1.1.1.1') to get full details.
    IMPORTANT: Even if full details are fetched, ONLY output basic info (ID, Name, Class, Reaction) to the user unless they explicitly asked for genes, pathways, or other massive fields.
    Use filter_field (e.g. 'products') and filter_value to strictly filter the initial matches.
    """
    # query: The search term
    # fetch_details: If True, fetches full details. Use ONLY for a specific EC number.
    # filter_field: Optional field name to filter by (e.g. 'products')
    # filter_value: Optional value that must be present in the filter_field
    print(f"\n[Tool: search_enzyme] Searching for: '{query}'")
    matches = search_kegg_enzyme(query, filter_field, filter_value)
    if not matches:
        return {"error": f"No enzyme found for '{query}'"}
    
    is_id = query.replace(".", "").isdigit()
    if fetch_details or is_id:
        return get_kegg_enzyme_details(matches[0]["id"])
    return {"matches": matches[:10]}

@tool
def tool_search_module(query: str, fetch_details: bool = False, filter_field: str = "", filter_value: str = ""):
    """
    Searches for a KEGG Module.
    Set fetch_details=True ONLY if query is a specific Module ID (e.g. 'M00001') to get full details.
    IMPORTANT: Even if full details are fetched, ONLY output basic info (ID, Name, Class) to the user unless they explicitly asked for genes, pathways, or other massive fields.
    Use filter_field and filter_value to strictly filter the initial matches.
    """
    # query: The search term
    # fetch_details: If True, fetches full details. Use ONLY for a specific Module ID.
    # filter_field: Optional field name to filter by
    # filter_value: Optional value that must be present in the filter_field
    print(f"\n[Tool: search_module] Searching for: '{query}'")
    matches = search_kegg_module(query, filter_field, filter_value)
    if not matches:
        return {"error": f"No module found for '{query}'"}
    
    is_id = query.upper().startswith("M") and len(query) > 1 and query[1].isdigit()
    if fetch_details or is_id:
        return get_kegg_module_details(matches[0]["id"])
    return {"matches": matches[:10]}

@tool
def tool_search_reaction(query: str, fetch_details: bool = False, filter_field: str = "", filter_value: str = ""):
    """
    Searches for a KEGG Reaction.
    Set fetch_details=True ONLY if query is a specific Reaction ID (e.g. 'R00001') to get full details.
    IMPORTANT: Even if full details are fetched, ONLY output basic info (ID, Name, Equation) to the user unless they explicitly asked for genes, pathways, or other massive fields.
    Use filter_field (e.g. 'products') and filter_value to strictly filter the initial matches.
    
    CRITICAL: Do NOT use this tool if you are already searching within a specific pathway 
    using tool_pathway_trace_compound.
    """
    # query: The search term
    # fetch_details: If True, fetches full details. Use ONLY for a specific Reaction ID.
    # filter_field: Optional field name to filter by (e.g. 'products')
    # filter_value: Optional value that must be present in the filter_field
    print(f"\n[Tool: search_reaction] Searching for: '{query}'")
    matches = search_kegg_reaction(query, filter_field, filter_value)
    if not matches:
        return {"error": f"No reaction found for '{query}'"}
    
    is_id = query.upper().startswith("R") and len(query) > 1 and query[1].isdigit()
    if fetch_details or is_id:
        return get_kegg_reaction_details(matches[0]["id"])
    return {"matches": matches[:10]}


@tool
def tool_get_pathway_linked_data(pathway_id: str, link_type: str):
    """
    Use this tool when the user asks for ALL reactions, ALL enzymes (EC numbers), or ALL compounds
    belonging to a specific pathway (e.g. 'list all reactions in glycolysis').
    It calls the KEGG Link API to get the IDs and then cross-references them with local files for full details.
    link_type must be one of:
      'rn'  -> returns all Reactions (with equations) in the pathway
      'ec'  -> returns all Enzyme EC numbers in the pathway
      'cpd' -> returns all Compounds (metabolites) in the pathway
    For pathway_id, use the KEGG map ID (e.g. 'map00010' for Glycolysis, 'map00020' for TCA cycle).
    To find the correct map ID, first call tool_search_kegg_pathway with just_list=True.
    
    CRITICAL: Do NOT use this tool if the user is asking to trace a specific compound or filter reactions 
    by a compound (e.g., 'reactions where pyruvate is a product'). In that case, use 
    tool_pathway_trace_compound INSTEAD. Calling both is a waste of API resources.
    """
    # pathway_id: KEGG map ID (e.g. map00010). Use tool_search_kegg_pathway to find the right ID.
    # link_type: 'rn' for reactions, 'ec' for enzymes, 'cpd' for compounds.
    print(f"\n[Tool: get_pathway_linked_data] pathway='{pathway_id}', type='{link_type}'")
    data = get_pathway_linked_data(pathway_id, link_type)
    
    # Safety truncation to prevent Groq API 8000 token limit crashes
    if "results" in data and len(data["results"]) > 50:
        data["results"] = data["results"][:50]
        data["warning"] = "Results truncated to 50 items to prevent rate limits. The pathway has more items."
        
    return data


# ---- PATHWAY ANALYTICS TOOLS (Zero LLM token processing) ---- #

@tool
def tool_pathway_summarize(pathway_id: str):
    """
    Returns a full statistical summary of a KEGG pathway: total reactions,
    reversible vs irreversible count, coverage (how many reactions have an EC),
    orphan (enzyme-unknown) reactions, EC class breakdown, and top 5 hub metabolites.
    Use this when the user asks for a general overview or statistics of a pathway.
    """
    # pathway_id: KEGG map ID such as map00010 (Glycolysis) or map00020 (TCA cycle)
    print(f"\n[Tool: pathway_summarize] pathway='{pathway_id}'")
    return pa.summarize_pathway(pathway_id)


@tool
def tool_pathway_isozymes(pathway_id: str):
    """
    Finds reactions in a pathway that are catalyzed by MORE than one EC number.
    These are called isozymes or redundant enzymes.
    Use this when the user asks: 'Which ECs catalyze the same reaction?' or
    'Which reactions have multiple enzymes?' in a given pathway.
    """
    # pathway_id: KEGG map ID such as map00010 (Glycolysis)
    print(f"\n[Tool: pathway_isozymes] pathway='{pathway_id}'")
    return pa.find_isozymes(pathway_id)


@tool
def tool_pathway_irreversible(pathway_id: str):
    """
    Returns only the irreversible reactions (one-directional, using =>) in a pathway.
    Irreversible reactions are often rate-limiting or regulatory steps.
    Use this when the user asks: 'Which reactions are irreversible?' or
    'What are the rate-limiting steps in this pathway?'
    """
    # pathway_id: KEGG map ID such as map00010 (Glycolysis)
    print(f"\n[Tool: pathway_irreversible] pathway='{pathway_id}'")
    return pa.find_irreversible_reactions(pathway_id)


@tool
def tool_pathway_orphan_reactions(pathway_id: str):
    """
    Finds reactions in a pathway that have NO EC number assigned.
    These are 'orphan' reactions where the chemistry is known but the enzyme
    is not yet characterised.
    Use this when the user asks: 'Which reactions lack an enzyme?' or
    'Are there any uncharacterised steps in this pathway?'
    """
    # pathway_id: KEGG map ID such as map00010 (Glycolysis)
    print(f"\n[Tool: pathway_orphan_reactions] pathway='{pathway_id}'")
    return pa.find_orphan_reactions(pathway_id)


@tool
def tool_pathway_ec_class_breakdown(pathway_id: str):
    """
    Groups all EC numbers in the pathway by their top-level enzyme class
    (Oxidoreductases, Transferases, Hydrolases, Lyases, Isomerases, Ligases, Translocases).
    Use this when the user asks: 'What types of reactions dominate this pathway?' or
    'Break down the enzymes in this pathway by class.'
    """
    # pathway_id: KEGG map ID such as map00010 (Glycolysis)
    print(f"\n[Tool: pathway_ec_class_breakdown] pathway='{pathway_id}'")
    return pa.get_ec_class_breakdown(pathway_id)


@tool
def tool_pathway_hub_metabolites(pathway_id: str, top_n: int = 10):
    """
    Finds the most frequently appearing compounds (hub metabolites) across all
    reactions in the pathway. Hub metabolites like ATP, NAD+, or CoA appear
    in many reactions and act as metabolic currencies.
    Use this when the user asks: 'Which compounds are most central in this pathway?'
    or 'What are the key metabolites here?'
    """
    # pathway_id: KEGG map ID such as map00010 (Glycolysis)
    # top_n: Number of top metabolites to return (default 10)
    print(f"\n[Tool: pathway_hub_metabolites] pathway='{pathway_id}', top={top_n}")
    return pa.find_hub_metabolites(pathway_id, top_n)


@tool
def tool_pathway_compound_intersection(pathway_ids: list[str]):
    """
    Finds compounds shared between 2 or more KEGG pathways. Useful for finding
    metabolic crossroads between pathways.
    Use this when the user asks: 'Which compounds are shared between Glycolysis
    and TCA cycle?' or 'What do these three pathways have in common?'
    """
    # pathway_ids: List of KEGG map IDs (e.g. ['map00010', 'map00020'])
    print(f"\n[Tool: pathway_compound_intersection] {pathway_ids}")
    return pa.find_pathway_compound_intersection(pathway_ids)


@tool
def tool_pathway_ec_intersection(pathway_ids: list[str]):
    """
    Finds enzymes (EC numbers) shared between 2 or more KEGG pathways. Shared enzymes
    can be dual-function enzymes or regulatory branch points.
    Use this when the user asks: 'Which enzymes are shared between these pathways?'
    """
    # pathway_ids: List of KEGG map IDs (e.g. ['map00010', 'map00020'])
    print(f"\n[Tool: pathway_ec_intersection] {pathway_ids}")
    return pa.find_pathway_ec_intersection(pathway_ids)


@tool
def tool_pathway_trace_compound(pathway_id: str, compound_name: str):
    """
    Traces a specific compound through a pathway. Shows exactly which reactions
    consume it (as substrate) and which produce it (as product).
    Use this when the user asks: 'Where is pyruvate consumed in Glycolysis?' or
    'Which reactions produce ATP in this pathway?'
    
    CRITICAL: This tool automatically fetches all reactions internally. Do NOT simultaneously 
    call tool_get_pathway_linked_data or tool_search_reaction when using this tool.
    """
    # pathway_id: KEGG map ID such as map00010 (Glycolysis)
    # compound_name: Partial or full name of the compound (e.g. 'pyruvate', 'ATP')
    print(f"\n[Tool: pathway_trace_compound] pathway='{pathway_id}', compound='{compound_name}'")
    return pa.trace_compound_in_pathway(pathway_id, compound_name)


# ---- REACTOME ANALYTICS TOOLS ---- #

@tool
def tool_reactome_summarize(pathway_name: str):
    """
    Returns a full statistical summary of a Reactome pathway: total reactions,
    disease vs normal reactions, compartments, unique genes, and hub metabolites.
    """
    print(f"\n[Tool: reactome_summarize] pathway='{pathway_name}'")
    return ra.summarize_pathway(pathway_name)

@tool
def tool_reactome_isozymes(pathway_name: str):
    """
    Finds reactions in a Reactome pathway that are catalyzed by multiple genes.
    """
    print(f"\n[Tool: reactome_isozymes] pathway='{pathway_name}'")
    return ra.find_reactions_with_multiple_genes(pathway_name)

@tool
def tool_reactome_irreversible(pathway_name: str):
    """
    Returns only the irreversible reactions in a Reactome pathway.
    """
    print(f"\n[Tool: reactome_irreversible] pathway='{pathway_name}'")
    return ra.find_irreversible_reactions(pathway_name)

@tool
def tool_reactome_orphan_reactions(pathway_name: str):
    """
    Finds reactions in a Reactome pathway that have NO genes or catalysts assigned.
    """
    print(f"\n[Tool: reactome_orphan_reactions] pathway='{pathway_name}'")
    return ra.find_orphan_reactions(pathway_name)

@tool
def tool_reactome_hub_metabolites(pathway_name: str, top_n: int = 10):
    """
    Finds the most frequently appearing compounds (hub metabolites) in a Reactome pathway.
    """
    print(f"\n[Tool: reactome_hub_metabolites] pathway='{pathway_name}', top={top_n}")
    return ra.find_hub_metabolites(pathway_name, top_n)

@tool
def tool_reactome_compound_intersection(pathway_names: list[str]):
    """
    Finds compounds shared between 2 or more Reactome pathways.
    """
    print(f"\n[Tool: reactome_compound_intersection] {pathway_names}")
    return ra.find_pathway_compound_intersection(pathway_names)

@tool
def tool_reactome_gene_intersection(pathway_names: list[str]):
    """
    Finds genes shared between 2 or more Reactome pathways.
    """
    print(f"\n[Tool: reactome_gene_intersection] {pathway_names}")
    return ra.find_pathway_gene_intersection(pathway_names)

@tool
def tool_reactome_trace_compound(pathway_name: str, compound_name: str):
    """
    Traces a specific compound through a Reactome pathway. Shows exactly which 
    reactions consume it and which produce it.
    """
    print(f"\n[Tool: reactome_trace_compound] pathway='{pathway_name}', compound='{compound_name}'")
    return ra.trace_compound(pathway_name, compound_name)

@tool
def tool_reactome_disease_reactions(pathway_name: str):
    """
    Finds reactions in a Reactome pathway that are specifically related to a disease or mutation.
    """
    print(f"\n[Tool: reactome_disease_reactions] pathway='{pathway_name}'")
    return ra.find_disease_reactions(pathway_name)

@tool
def tool_reactome_get_reactions(pathway_name: str, basic_only: bool = True):
    """
    Returns the complete list of ALL reactions for a given Reactome pathway.
    If basic_only is True, returns only ID, Name, and Equation.
    If the user explicitly asks for genes, compartments, or descriptions, set basic_only to False.
    """
    print(f"\n[Tool: reactome_get_reactions] pathway='{pathway_name}', basic={basic_only}")
    reactions = ra._get_pathway_reactions(pathway_name)
    
    # Return a clean list of dictionaries with important user-facing fields
    full_rxns = []
    for r in reactions:
        rxn = {
            "Reaction_ID": r.get("Reaction_ID", ""),
            "Reaction_Name": r.get("Reaction_Name", ""),
            "Equation": r.get("Equation", "")
        }
        if not basic_only:
            rxn["Genes"] = r.get("Genes", "")
            rxn["Compartment"] = r.get("Compartment", "")
            rxn["Summation"] = r.get("Summation", "")
        full_rxns.append(rxn)
        
    # Truncate if massive to prevent token limits
    if len(full_rxns) > 50:
        full_rxns = full_rxns[:50]
        return {"warning": "Truncated to 50 reactions.", "reactions": full_rxns}
        
    return {"total_reactions": len(full_rxns), "reactions": full_rxns}


@tool
def tool_compare_kegg_reactome_reactions(kegg_pathway_id: str, reactome_pathway_name: str):
    """
    Finds common reactions between a KEGG pathway and a Reactome pathway using ChEBI IDs as a universal bridge.
    This is the ONLY correct tool to use when the user asks for 'common reactions', 'shared reactions',
    or wants to compare reactions between KEGG and Reactome for the same biological pathway.
    IMPORTANT: Requires the KEGG map ID (e.g. 'map00010') and the exact Reactome pathway name (e.g. 'Glycolysis').
    Always search for the pathway first if you do not have these IDs.
    """
    # kegg_pathway_id: The KEGG map ID of the pathway (e.g. 'map00010' for Glycolysis)
    # reactome_pathway_name: The exact Reactome pathway name (e.g. 'Glycolysis')
    print(f"\n[Tool: compare_kegg_reactome_reactions] KEGG='{kegg_pathway_id}', Reactome='{reactome_pathway_name}'")
    return ra.find_common_reactions_cross_db(kegg_pathway_id, reactome_pathway_name)


@tool
def tool_resolve_chemical(name: str, target_db: str = "all"):
    """
    Resolves a chemical name to its exact universal structure (InChIKey) and maps it 
    across all major biological databases (KEGG, ChEBI, HMDB, DrugBank, PubChem).
    Also returns a perfectly cleaned list of known synonyms.
    target_db can be 'all', 'KEGG', 'ChEBI', 'HMDB', 'DrugBank', or 'PubChem'.
    """
    print(f"\n[Tool: resolve_chemical] Resolving '{name}' (target: {target_db})")
    return cr.resolve_chemical(name, target_db=target_db)


@tool
def tool_query_rhea(rhea_id: str):
    """
    Directly query a Rhea Master ID to get its chemical equation, ChEBI identifiers, 
    and cross-database mappings (KEGG, Reactome, UniProt, EC).
    Input should be the numeric ID or 'RHEA:10000' format.
    """
    print(f"\n[Tool: query_rhea] Fetching details for Rhea ID: '{rhea_id}'")
    return ra.query_rhea_master_id(rhea_id)


@tool
def tool_pathway_rhea_mapping(kegg_pathway_id: str):
    """
    Gets ALL reactions in a KEGG pathway and instantly maps them to their Rhea Master IDs and UniProt IDs.
    IMPORTANT: Requires the KEGG map ID (e.g. 'map00010').
    Use this when asked for "all Rhea IDs in a pathway" or "Rhea mapping for a pathway".
    """
    print(f"\n[Tool: pathway_rhea_mapping] KEGG='{kegg_pathway_id}'")
    return ra.get_pathway_rhea_mapping(kegg_pathway_id)


# ---- AGENT BUILD FUNCTION ---- #

def build_ncbi_kegg_graph():
    # Builds and returns the dynamic ReAct tool-calling agent.
    # The LLM will decide which tools to call and in what order based on the user query.

    from core.ai_model_wrapper import get_llm
    llm = get_llm()

    tools = [
        tool_search_kegg_pathway,
        tool_search_kegg_disease,
        tool_get_ncbi_gene_info,
        tool_fetch_clinvar_mutations,
        tool_search_compound,
        tool_search_ko,
        tool_search_enzyme,
        tool_search_module,
        tool_search_reaction,
        tool_get_pathway_linked_data,
        # Deterministic analytics tools (zero LLM token processing)
        tool_pathway_summarize,
        tool_pathway_isozymes,
        tool_pathway_irreversible,
        tool_pathway_orphan_reactions,
        tool_pathway_ec_class_breakdown,
        tool_pathway_hub_metabolites,
        tool_pathway_compound_intersection,
        tool_pathway_ec_intersection,
        tool_pathway_trace_compound,
        # Reactome analytics tools
        tool_reactome_summarize,
        tool_reactome_isozymes,
        tool_reactome_irreversible,
        tool_reactome_orphan_reactions,
        tool_reactome_hub_metabolites,
        tool_reactome_compound_intersection,
        tool_reactome_gene_intersection,
        tool_reactome_trace_compound,
        tool_reactome_disease_reactions,
        tool_reactome_get_reactions,
        tool_compare_kegg_reactome_reactions,
        tool_resolve_chemical,
        tool_query_rhea,
        tool_pathway_rhea_mapping,
    ]

    system_prompt = (
        "You are a biomedical research assistant. Use tools to search NCBI, KEGG, and ClinVar.\n\n"
        "### CRITICAL RULES:\n"
        "1. **OFF-TOPIC**: You are a specialized agent. If a query is COMPLETELY unrelated to biology, chemistry, medicine, genetics, or diseases (e.g., asking for recipes or general chatter), refuse it. However, if the query contains ANY biochemical terms (like 'glycolysis', 'rhea', 'reactions', 'kegg'), you MUST process it and NEVER refuse it.\n"
        "2. **FORMATTING**: No LaTeX math. Use plain text.\n"
        "3. **KNOWLEDGE**: For database queries, answer ONLY from tool output and never invent fields. **CRITICAL: NEVER summarize, paraphrase, or rewrite descriptions retrieved from KEGG or Reactome (e.g., pathway or disease overviews). You MUST output the exact word-for-word text provided by the tool.** For general conceptual questions (e.g. 'what are modules?', 'explain glycolysis' without asking for database info), you may use your internal knowledge and answer concisely.\n"
        "4. **TOOL SELECTION & PAGINATION**:\n"
        "   - Disease/Pathway genes: call `tool_search_kegg_disease` / `tool_search_kegg_pathway` ONCE. Do NOT loop `tool_get_ncbi_gene_info`.\n"
        "   - `tool_get_ncbi_gene_info`: use only for a single specific gene the user explicitly names.\n"
        "   - `tool_fetch_clinvar_mutations`: use only if explicitly asked for mutations. (Note: Many Reactome pathways have long names like 'APC truncation mutants have impaired AXIN binding'. If a query is a full descriptive phrase about mutants, search it as a Reactome pathway FIRST, rather than assuming it is a ClinVar mutation query!)\n"
        "   - Do NOT auto-fetch additional pages. Return the first page and stop.\n"
        "5. **TRUNCATION**: Never truncate lists. Output every item returned by the tool.\n"
        "6. **EFFICIENT FETCHING**: Do NOT auto-fetch genes/compounds/drugs unless explicitly asked. For general pathway info, only fetch 'description'. For generic 'list' queries (e.g., 'list all pathways', 'list a few pathways from Reactome'), you MUST strictly set `query='all'`, `just_list=True`, and `fields_needed=[]`. NEVER pick a random specific pathway (like apoptosis) to deep-dive. NEVER output descriptions when the user just wants a list of pathways.\n"
        "7. **WARNINGS**: If a tool returns '(Fallback Match)', explicitly warn the user.\n"
        "8. **RELEVANCE**: Discard irrelevant tool results silently. If data is a bad match, reply: 'I could not find relevant information.'\n"
        "9. **STRUCTURED DATA**: For KO/Enzyme/Module/Reaction queries, read explicit fields (substrates, products, equation, genes, pathways) and output ONLY what the user asked for. If the user just provides an ID or asks for general/basic info, provide only a CONCISE summary (e.g. ID, Name, Class, and Equation/Reaction). Do NOT dump the full list of genes, orthology, pathways, substrates, or products unless explicitly requested.\n"
        "10. **SMART FILTERING**: To find reactions/enzymes/modules involving a specific compound (e.g. 'pyruvate'), ALWAYS provide the compound name in the 'query' parameter, and use 'filter_field' (e.g. 'products') and 'filter_value' (e.g. 'pyruvate'). Never omit 'query'.\n"
        "11. **LOCAL DATA PRIORITY**: Reactions, Enzymes, Modules, and KOs are parsed perfectly from local offline files. NEVER manually loop through a huge list of reaction/enzyme IDs returned by a compound search. ALWAYS use the compound NAME directly in the specialized search tool (e.g., `tool_search_reaction`).\n"
        "12. **PATHWAY DATA**: To get ALL reactions/enzymes/compounds in a KEGG pathway, use `tool_get_pathway_linked_data`. For Reactome data, use `tool_search_kegg_pathway` with `database_preference='reactome'`.\n"
        "13. **CRITICAL ANALYTICS OVERRIDE**: NEVER manually download raw pathway data to trace compounds, count enzymes, or check reversibility. YOU MUST STRICTLY use the dedicated analytics tools:\n"
        "    - 'same reaction / multiple ECs / isozymes' -> `tool_pathway_isozymes`\n"
        "    - 'irreversible / rate-limiting steps' -> `tool_pathway_irreversible`\n"
        "    - 'no enzyme / orphan reactions' -> `tool_pathway_orphan_reactions`\n"
        "    - 'overview / statistics / summary of pathway' -> `tool_pathway_summarize`\n"
        "    - 'enzyme type breakdown' -> `tool_pathway_ec_class_breakdown`\n"
        "    - 'hub metabolites / central metabolites' -> `tool_pathway_hub_metabolites`\n"
        "    - 'shared compounds between two pathways' -> `tool_pathway_compound_intersection`\n"
        "    - 'shared enzymes between two pathways' -> `tool_pathway_ec_intersection`\n"
        "    - 'where is compound X / gimme reactions for compound X' -> `tool_pathway_trace_compound`\n"
        "    *IMPORTANT*: All KEGG analytics tools require a pathway ID (e.g., 'map00010'). If you only have the pathway name, call `tool_search_kegg_pathway` FIRST to get the ID. For Reactome analytics, you can use EITHER the exact pathway name or the Reactome ID (e.g., 'R-HSA-9663891') in the tool_reactome_* tools.\n"
        "    - If asked for disease-specific reactions in a Reactome pathway, use `tool_reactome_disease_reactions`.\n"
        "    - If asked to list reactions in a Reactome pathway, or if the user simply provides a Reactome ID (e.g., 'R-HSA-9663891'), they want to see its basic reactions. Use `tool_reactome_get_reactions`. By default, only output basic info (ID and Equation). NEVER output descriptions, full genes, compartments, or summations unless explicitly asked.\n"
        "    - If asked for 'common reactions' / 'shared reactions' / 'reactions in both KEGG and Reactome' for any pathway -> ALWAYS use `tool_compare_kegg_reactome_reactions`. NEVER manually compare reaction lists yourself.\n"
        "    - If asked to 'resolve', 'identify', find 'synonyms', or map a chemical across databases (e.g. 'what is the KEGG ID for pyruvate?'), ALWAYS use `tool_resolve_chemical`.\n"
        "    - If asked about a Rhea ID (e.g., 'What is Rhea 10000?'), ALWAYS use `tool_query_rhea` to get its equation and cross-database mappings.\n"
        "    - If asked to get Rhea IDs for an entire pathway (e.g. 'get all reactions from glycolysis and their rhea ids'), ALWAYS use `tool_pathway_rhea_mapping`.\n"
        "14. **SOURCES, REFERENCES & FALLBACKS**: You MUST explicitly and accurately state the exact source of ALL data you provide. If a tool uses a fallback source (like fetching a Wikipedia summary because a KEGG disease description was missing), you MUST explicitly state that the data is from Wikipedia (or the actual fallback source) and briefly explain WHY (e.g., 'Source: Wikipedia (Fallback used because KEGG lacked a description)'). Do NOT falsely claim data is from KEGG or Reactome if the tool pulled it from Wikipedia, HGNC, or another fallback. This strict citation rule applies to EVERYTHING: diseases, pathways, genes, compounds, etc. If explaining concepts purely from your own internal training data, explicitly state 'Source: Internal LLM Knowledge'.\n\n"
        "Format answers as clean Markdown. No filler sentences. No 'Further Reading' unless asked."
    )


    agent = create_react_agent(llm, tools, prompt=system_prompt)
    # Hard cap: max 10 tool calls per query to prevent infinite looping
    return agent.with_config({"recursion_limit": 10})


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
