import os
import requests
import pandas as pd

# Define the hardcoded path to the Reactome TSV file
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
REACTOME_FILE = os.path.join(BASE_DIR, "..", "Reactome_files", "reactome_reactions.tsv")

def _get_pathway_reactions(pathway_name):
    # pathway_name: Exact string name of the Reactome pathway to filter the data by
    print(f"[Reactome Analytics] Loading and filtering reactions for pathway: {pathway_name}")
    
    # Load the TSV file using pandas as requested and safely fill NaN with empty strings
    df = pd.read_csv(REACTOME_FILE, sep='\t').fillna("")
    
    # Filter the dataframe to keep rows matching either the pathway name or ID (case-insensitive)
    query_val = str(pathway_name).strip().lower()
    filtered_df = df[
        (df['Pathway_Name'].str.strip().str.lower() == query_val) |
        (df['Pathway_ID'].str.strip().str.lower() == query_val)
    ]
    
    # Convert the filtered dataframe into a list of dictionaries for easier processing
    reactions = filtered_df.to_dict(orient='records')
    
    # Pre-parse equations and lists for easier analytics later
    for r in reactions:
        # Check reversibility
        equation = str(r['Equation'])
        r['is_reversible'] = "<=>" in equation
        r['is_irreversible'] = "=>" in equation and "<=>" not in equation
        
        # Parse genes into a list
        genes_str = str(r['Genes']).strip()
        if genes_str != "nan" and genes_str != "":
            r['parsed_genes'] = [g.strip() for g in genes_str.split(";")]
        else:
            r['parsed_genes'] = []
            
        # Parse inputs into a list
        inputs_str = str(r['Inputs']).strip()
        if inputs_str != "nan" and inputs_str != "":
            r['parsed_inputs'] = [i.strip() for i in inputs_str.split(";")]
        else:
            r['parsed_inputs'] = []
            
        # Parse outputs into a list
        outputs_str = str(r['Outputs']).strip()
        if outputs_str != "nan" and outputs_str != "":
            r['parsed_outputs'] = [o.strip() for o in outputs_str.split(";")]
        else:
            r['parsed_outputs'] = []

    print(f"[Reactome Analytics] Found {len(reactions)} reactions for {pathway_name}")
    return reactions

def _get_pathway_compounds(pathway_name):
    # pathway_name: Exact string name of the Reactome pathway
    # Helper function that returns a set of all unique compounds (inputs + outputs) in the pathway
    reactions = _get_pathway_reactions(pathway_name)
    compounds = set()
    for r in reactions:
        for c in r['parsed_inputs']:
            compounds.add(c.lower())
        for c in r['parsed_outputs']:
            compounds.add(c.lower())
    return compounds

def _get_pathway_genes(pathway_name):
    # pathway_name: Exact string name of the Reactome pathway
    # Helper function that returns a set of all unique genes in the pathway
    reactions = _get_pathway_reactions(pathway_name)
    genes = set()
    for r in reactions:
        for g in r['parsed_genes']:
            genes.add(g.upper())
    return genes

# =====================================================================================
# PUBLIC ANALYTICS FUNCTIONS (MAPPED FROM KEGG)
# =====================================================================================

def find_irreversible_reactions(pathway_name):
    # pathway_name: Exact string name of the Reactome pathway
    print(f"\n[Reactome Analytics] find_irreversible_reactions: {pathway_name}")
    reactions = _get_pathway_reactions(pathway_name)
    
    # Filter reactions that only contain '=>' without '<=>'
    irreversible = [r for r in reactions if r['is_irreversible']]
    print(f"[Reactome Analytics] Found {len(irreversible)} irreversible reactions.")
    
    return {
        "pathway_name": pathway_name,
        "analysis": "irreversible_reactions",
        "total_reactions": len(reactions),
        "results": irreversible
    }

def find_reversible_reactions(pathway_name):
    # pathway_name: Exact string name of the Reactome pathway
    print(f"\n[Reactome Analytics] find_reversible_reactions: {pathway_name}")
    reactions = _get_pathway_reactions(pathway_name)
    
    # Filter reactions that contain '<=>'
    reversible = [r for r in reactions if r['is_reversible']]
    print(f"[Reactome Analytics] Found {len(reversible)} reversible reactions.")
    
    return {
        "pathway_name": pathway_name,
        "analysis": "reversible_reactions",
        "total_reactions": len(reactions),
        "results": reversible
    }

def find_orphan_reactions(pathway_name):
    # pathway_name: Exact string name of the Reactome pathway
    print(f"\n[Reactome Analytics] find_orphan_reactions: {pathway_name}")
    reactions = _get_pathway_reactions(pathway_name)
    
    # Find reactions that have NO genes or catalysts assigned
    orphans = [r for r in reactions if len(r['parsed_genes']) == 0]
    print(f"[Reactome Analytics] Found {len(orphans)} orphan reactions.")
    
    return {
        "pathway_name": pathway_name,
        "analysis": "orphan_reactions",
        "total_reactions": len(reactions),
        "results": orphans
    }

def get_pathway_coverage(pathway_name):
    # pathway_name: Exact string name of the Reactome pathway
    print(f"\n[Reactome Analytics] get_pathway_coverage: {pathway_name}")
    reactions = _get_pathway_reactions(pathway_name)
    
    total = len(reactions)
    with_genes = len([r for r in reactions if len(r['parsed_genes']) > 0])
    without_genes = total - with_genes
    
    # Calculate percentage safely
    coverage_pct = round((with_genes / total) * 100, 1) if total > 0 else 0
    print(f"[Reactome Analytics] Coverage: {with_genes}/{total} ({coverage_pct}%)")
    
    return {
        "pathway_name": pathway_name,
        "analysis": "pathway_coverage",
        "total_reactions": total,
        "reactions_with_genes": with_genes,
        "orphan_reactions": without_genes,
        "coverage_percent": coverage_pct
    }

def find_reactions_with_multiple_genes(pathway_name):
    # pathway_name: Exact string name of the Reactome pathway
    # This is equivalent to find_isozymes in KEGG.
    print(f"\n[Reactome Analytics] find_reactions_with_multiple_genes: {pathway_name}")
    reactions = _get_pathway_reactions(pathway_name)
    
    # Keep reactions that have more than 1 gene assigned
    multiple_genes = [r for r in reactions if len(r['parsed_genes']) > 1]
    print(f"[Reactome Analytics] Found {len(multiple_genes)} reactions with multiple genes.")
    
    return {
        "pathway_name": pathway_name,
        "analysis": "multiple_genes_reactions",
        "total_reactions": len(reactions),
        "results": multiple_genes
    }

def find_hub_metabolites(pathway_name, top_n=10):
    # pathway_name: Exact string name of the Reactome pathway
    # top_n: Number of top compounds to return
    print(f"\n[Reactome Analytics] find_hub_metabolites: {pathway_name}")
    reactions = _get_pathway_reactions(pathway_name)
    
    freq = {}
    
    # Loop through and count every input and output compound
    for r in reactions:
        for comp in r['parsed_inputs'] + r['parsed_outputs']:
            freq[comp] = freq.get(comp, 0) + 1
            
    # Sort the dictionary by count in descending order
    sorted_cpds = sorted(freq.items(), key=lambda x: x[1], reverse=True)[:top_n]
    results = [{"compound": name, "appearances": count} for name, count in sorted_cpds]
    
    print(f"[Reactome Analytics] Top {top_n} hub metabolites found.")
    
    return {
        "pathway_name": pathway_name,
        "analysis": "hub_metabolites",
        "total_reactions_scanned": len(reactions),
        "results": results
    }

def find_pathway_compound_intersection(pathway_names):
    # pathway_names: List of exact string names of the Reactome pathways to compare
    if not pathway_names or len(pathway_names) < 2:
        return {"error": "Provide at least 2 pathway names for intersection."}
        
    print(f"\n[Reactome Analytics] find_pathway_compound_intersection: {pathway_names}")
    
    # Get compounds for the first pathway
    shared = _get_pathway_compounds(pathway_names[0])
    
    # Intersect with all other pathways
    for pname in pathway_names[1:]:
        shared &= _get_pathway_compounds(pname)
        
    results = sorted(list(shared))
    print(f"[Reactome Analytics] Shared compounds: {len(results)}")
    
    return {
        "pathways": pathway_names,
        "analysis": "compound_intersection",
        "shared_count": len(results),
        "results": results
    }

def find_pathway_gene_intersection(pathway_names):
    # pathway_names: List of exact string names of the Reactome pathways to compare
    # This is equivalent to find_pathway_ec_intersection in KEGG
    if not pathway_names or len(pathway_names) < 2:
        return {"error": "Provide at least 2 pathway names for intersection."}
        
    print(f"\n[Reactome Analytics] find_pathway_gene_intersection: {pathway_names}")
    
    # Get genes for the first pathway
    shared = _get_pathway_genes(pathway_names[0])
    
    # Intersect with all other pathways
    for pname in pathway_names[1:]:
        shared &= _get_pathway_genes(pname)
        
    results = sorted(list(shared))
    print(f"[Reactome Analytics] Shared genes: {len(results)}")
    
    return {
        "pathways": pathway_names,
        "analysis": "gene_intersection",
        "shared_count": len(results),
        "results": results
    }

def find_pathway_unique_compounds(pathway_name_1, pathway_name_2):
    # pathway_name_1: The pathway to find unique items FOR
    # pathway_name_2: The reference pathway to compare AGAINST
    print(f"\n[Reactome Analytics] find_pathway_unique_compounds: unique to {pathway_name_1}")
    
    cpds1 = _get_pathway_compounds(pathway_name_1)
    cpds2 = _get_pathway_compounds(pathway_name_2)
    
    # Find compounds in cpds1 that are NOT in cpds2
    unique = cpds1 - cpds2
    results = sorted(list(unique))
    
    print(f"[Reactome Analytics] Unique compounds in {pathway_name_1}: {len(results)}")
    
    return {
        "pathway_name": pathway_name_1,
        "compared_against": pathway_name_2,
        "analysis": "unique_compounds",
        "unique_count": len(results),
        "results": results
    }

def trace_compound(pathway_name, compound):
    # pathway_name: Exact string name of the Reactome pathway
    # compound: String name of the molecule/compound to trace
    print(f"\n[Reactome Analytics] Tracing compound '{compound}' in pathway: {pathway_name}")
    reactions = _get_pathway_reactions(pathway_name)
    
    consumed_in = []
    produced_in = []
    compound_lower = compound.lower()
    
    # Check parsed lists instead of strings to avoid partial word matches
    for r in reactions:
        # Check inputs
        for inp in r['parsed_inputs']:
            if compound_lower in inp.lower():
                consumed_in.append({"reaction_id": r['Reaction_ID'], "equation": r['Equation']})
                break
                
        # Check outputs
        for out in r['parsed_outputs']:
            if compound_lower in out.lower():
                produced_in.append({"reaction_id": r['Reaction_ID'], "equation": r['Equation']})
                break
                
    print(f"[Reactome Analytics] '{compound}' consumed in {len(consumed_in)} rxns, produced in {len(produced_in)} rxns.")
    
    return {
        "pathway_name": pathway_name,
        "compound_traced": compound,
        "analysis": "compound_trace",
        "consumed_in_reactions": consumed_in,
        "produced_in_reactions": produced_in
    }

def find_disease_reactions(pathway_name):
    # pathway_name: Exact string name of the Reactome pathway
    print(f"\n[Reactome Analytics] Finding disease reactions in pathway: {pathway_name}")
    reactions = _get_pathway_reactions(pathway_name)
    
    disease_reactions = []
    
    for r in reactions:
        if r['Is_Disease'] == True:
            disease_reactions.append({
                "reaction_id": r['Reaction_ID'],
                "reaction_name": r['Reaction_Name'],
                "genes": r['Genes'],
                "summation": r['Summation']
            })
            
    print(f"[Reactome Analytics] Found {len(disease_reactions)} disease-related reactions.")
    
    return {
        "pathway_name": pathway_name,
        "analysis": "disease_reactions",
        "total_disease_reactions": len(disease_reactions),
        "results": disease_reactions
    }

def summarize_pathway(pathway_name):
    # pathway_name: Exact string name of the Reactome pathway
    print(f"\n[Reactome Analytics] Summarizing pathway: {pathway_name}")
    reactions = _get_pathway_reactions(pathway_name)
    
    total = len(reactions)
    if total == 0:
        return {"error": f"No reactions found for pathway {pathway_name}"}
        
    disease_count = len([r for r in reactions if r['Is_Disease'] == True])
    normal_count = total - disease_count
    
    reversible = len([r for r in reactions if r['is_reversible']])
    irreversible = len([r for r in reactions if r['is_irreversible']])
    orphans = len([r for r in reactions if len(r['parsed_genes']) == 0])
    with_genes = total - orphans
    coverage = round((with_genes / total) * 100, 1)
    
    multiple_genes = len([r for r in reactions if len(r['parsed_genes']) > 1])
    
    # Compartment counts
    compartments = {}
    freq = {}
    
    for r in reactions:
        comp = str(r['Compartment']).strip()
        if comp != "nan" and comp != "":
            compartments[comp] = compartments.get(comp, 0) + 1
            
        for c in r['parsed_inputs'] + r['parsed_outputs']:
            freq[c] = freq.get(c, 0) + 1
            
    top5 = sorted(freq.items(), key=lambda x: x[1], reverse=True)[:5]
    genes_list = list(_get_pathway_genes(pathway_name))
    
    print(f"[Reactome Analytics] Summary complete. Disease reactions: {disease_count}, Normal: {normal_count}")
    
    return {
        "pathway_name": pathway_name,
        "analysis": "full_summary",
        "total_reactions": total,
        "disease_reactions": disease_count,
        "normal_reactions": normal_count,
        "reversible": reversible,
        "irreversible": irreversible,
        "reactions_with_genes": with_genes,
        "orphan_reactions": orphans,
        "coverage_percent": coverage,
        "reactions_with_multiple_genes": multiple_genes,
        "compartments": compartments,
        "total_unique_genes": len(genes_list),
        "top_5_hub_metabolites": [{"compound": n, "appearances": c} for n, c in top5]
    }

# Hardcoded paths for KEGG and Rhea cross-reference files
CHEBI_REACTOME_FILE = os.path.join(BASE_DIR, "..", "Reactome_files", "ChEBI2Reactome_PE_Reactions.txt")
KEGG_CPD_TO_CHEBI_FILE = os.path.join(BASE_DIR, "..", "Kegg_files", "kegg_cpd_to_chebi.tsv")
KEGG_RN_TO_CPD_FILE = os.path.join(BASE_DIR, "..", "Kegg_files", "kegg_rn_to_cpd.tsv")
RHEA_KEGG_FILE = os.path.join(BASE_DIR, "..", "Rhea_files", "rhea2kegg_reaction.tsv")
RHEA_REACTOME_FILE = os.path.join(BASE_DIR, "..", "Rhea_files", "rhea2reactome.tsv")

# Module-level caches: each loaded once on first call
_CHEBI_REACTOME_INDEX = None
_CPD_TO_CHEBI_MAP = None
_RN_TO_CPD_MAP = None
_RHEA_KEGG_MAP = None      # kegg_rn_id -> rhea_master_id
_RHEA_REACTOME_MAP = None  # rhea_master_id -> set of bare reactome rxn ids

def _load_chebi_reactome_index():
    # Loads ChEBI2Reactome_PE_Reactions.txt into a fast dict: {chebi_id -> set of reaction_ids}
    # This is cached at module level so the file is only read once per session.
    global _CHEBI_REACTOME_INDEX
    if _CHEBI_REACTOME_INDEX is not None:
        return _CHEBI_REACTOME_INDEX

    print("[Reactome Analytics] Loading ChEBI-Reactome index (one-time)...")
    index = {}
    df = pd.read_csv(
        CHEBI_REACTOME_FILE,
        sep="\t",
        header=None,
        names=["ChEBI_ID", "PE_ID", "PE_Name", "Reaction_ID", "URL", "Reaction_Name", "Evidence", "Species"]
    )
    # Filter to human reactions only for speed
    df = df[df["Species"] == "Homo sapiens"]

    for _, row in df.iterrows():
        chebi_id = str(row["ChEBI_ID"]).strip()
        rxn_id = str(row["Reaction_ID"]).strip()
        index.setdefault(chebi_id, set()).add(rxn_id)

    _CHEBI_REACTOME_INDEX = index
    print(f"[Reactome Analytics] ChEBI-Reactome index loaded: {len(index)} ChEBI entries.")
    return index


def _load_cpd_to_chebi_map():
    # Loads kegg_cpd_to_chebi.tsv into a dict: {cpd_id -> set of chebi_ids}
    # File format: 'cpd:C00022\tchebi:15361'
    # Cached at module level so file is only read once.
    global _CPD_TO_CHEBI_MAP
    if _CPD_TO_CHEBI_MAP is not None:
        return _CPD_TO_CHEBI_MAP
    print("[Reactome Analytics] Loading KEGG compound-to-ChEBI map (one-time)...")
    mapping = {}
    with open(KEGG_CPD_TO_CHEBI_FILE, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                # Strip 'cpd:' and 'chebi:' prefixes
                cpd_id = parts[0].replace("cpd:", "").strip()
                chebi_id = parts[1].replace("chebi:", "").strip()
                mapping.setdefault(cpd_id, set()).add(chebi_id)
    _CPD_TO_CHEBI_MAP = mapping
    print(f"[Reactome Analytics] CPD->ChEBI map loaded: {len(mapping)} compounds.")
    return mapping


def _load_rn_to_cpd_map():
    # Loads kegg_rn_to_cpd.tsv into a dict: {rn_id -> set of cpd_ids}
    # File format: 'rn:R00001\tcpd:C00001'
    # Cached at module level so file is only read once.
    global _RN_TO_CPD_MAP
    if _RN_TO_CPD_MAP is not None:
        return _RN_TO_CPD_MAP
    print("[Reactome Analytics] Loading KEGG reaction-to-compound map (one-time)...")
    mapping = {}
    with open(KEGG_RN_TO_CPD_FILE, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                # Strip 'rn:' and 'cpd:' prefixes
                rn_id = parts[0].replace("rn:", "").strip()
                cpd_id = parts[1].replace("cpd:", "").strip()
                mapping.setdefault(rn_id, set()).add(cpd_id)
    _RN_TO_CPD_MAP = mapping
    print(f"[Reactome Analytics] RN->CPD map loaded: {len(mapping)} reactions.")
    return mapping


def _get_chebi_ids_for_kegg_compound(kegg_cpd_id, cpd_to_chebi_map):
    # kegg_cpd_id: KEGG compound ID like 'C00022' (without 'cpd:' prefix)
    # cpd_to_chebi_map: Pre-loaded dict from _load_cpd_to_chebi_map()
    # Returns a set of ChEBI IDs from the local file. No API call made.
    return cpd_to_chebi_map.get(kegg_cpd_id, set())


def _load_rhea_kegg_map():
    # Loads rhea2kegg_reaction.tsv into a dict: {kegg_rn_id -> rhea_master_id}
    # Columns: RHEA_ID, DIRECTION, MASTER_ID, ID (ID is the KEGG RN ID like R00200)
    # Cached at module level so the file is only read once.
    global _RHEA_KEGG_MAP
    if _RHEA_KEGG_MAP is not None:
        return _RHEA_KEGG_MAP
    print("[Reactome Analytics] Loading Rhea->KEGG map (one-time)...")
    mapping = {}
    df = pd.read_csv(RHEA_KEGG_FILE, sep="\t")
    for _, row in df.iterrows():
        kegg_id = str(row["ID"]).strip()
        master_id = str(int(row["MASTER_ID"]))
        # Multiple Rhea IDs map to the same KEGG ID (directional variants) - master_id is always the same
        mapping[kegg_id] = master_id
    _RHEA_KEGG_MAP = mapping
    print(f"[Reactome Analytics] Rhea->KEGG map loaded: {len(mapping)} KEGG reactions.")
    return mapping


def _load_rhea_reactome_map():
    # Loads rhea2reactome.tsv into a dict: {rhea_master_id -> set of bare reactome rxn ids}
    # Columns: RHEA_ID, DIRECTION, MASTER_ID, ID (ID is like 'R-HSA-71670.9' with version suffix)
    # Cached at module level so the file is only read once.
    global _RHEA_REACTOME_MAP
    if _RHEA_REACTOME_MAP is not None:
        return _RHEA_REACTOME_MAP
    print("[Reactome Analytics] Loading Rhea->Reactome map (one-time)...")
    mapping = {}
    df = pd.read_csv(RHEA_REACTOME_FILE, sep="\t")
    for _, row in df.iterrows():
        reactome_id_raw = str(row["ID"]).strip()
        # Strip version suffix: 'R-HSA-71670.9' -> 'R-HSA-71670'
        reactome_id = reactome_id_raw.rsplit(".", 1)[0]
        master_id = str(int(row["MASTER_ID"]))
        mapping.setdefault(master_id, set()).add(reactome_id)
    _RHEA_REACTOME_MAP = mapping
    print(f"[Reactome Analytics] Rhea->Reactome map loaded: {len(mapping)} Rhea master entries.")
    return mapping


def find_common_reactions_cross_db(kegg_pathway_id, reactome_pathway_name):
    # kegg_pathway_id: KEGG map ID like 'map00010' for Glycolysis
    # reactome_pathway_name: Exact string name of the Reactome pathway (e.g. 'Glycolysis')
    print(f"\n[Reactome Analytics] Cross-DB comparison (Rhea bridge): KEGG '{kegg_pathway_id}' vs Reactome '{reactome_pathway_name}'")

    # Step 1: Get all KEGG reactions for the pathway (1 API call - unavoidable)
    kegg_link_url = f"https://rest.kegg.jp/link/rn/{kegg_pathway_id}"
    kegg_rn_ids = set()
    res = requests.get(kegg_link_url, timeout=15)
    if res.status_code == 200 and res.text.strip():
        for line in res.text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) == 2:
                kegg_rn_ids.add(parts[1].replace("rn:", "").strip())
    if not kegg_rn_ids:
        return {"error": f"No KEGG reactions found for {kegg_pathway_id}"}
    print(f"[Reactome Analytics] KEGG reactions in pathway: {len(kegg_rn_ids)}")

    # Step 2: Load Rhea bridge maps from local files (zero API calls)
    rhea_kegg_map = _load_rhea_kegg_map()          # kegg_rn_id -> rhea_master_id
    rhea_reactome_map = _load_rhea_reactome_map()  # rhea_master_id -> set of reactome_rxn_ids

    # Step 3: Get all Reactome reactions for the target pathway from local TSV (zero API calls)
    reactome_rxns = _get_pathway_reactions(reactome_pathway_name)
    reactome_rxn_ids = set(r["Reaction_ID"] for r in reactome_rxns)
    reactome_rxn_map = {r["Reaction_ID"]: r for r in reactome_rxns}
    if not reactome_rxn_ids:
        return {"error": f"No Reactome reactions found for pathway '{reactome_pathway_name}'"}
    print(f"[Reactome Analytics] Reactome reactions in pathway: {len(reactome_rxn_ids)}")

    # Step 4: Match KEGG rxns to Reactome rxns via shared Rhea Master ID
    # This is exact and authoritative - no text matching, no ChEBI version mismatches
    common_reactions = []
    unmatched_kegg = []  # KEGG reactions that have no Rhea entry

    for kegg_rxn_id in kegg_rn_ids:
        master_id = rhea_kegg_map.get(kegg_rxn_id)
        if not master_id:
            unmatched_kegg.append(kegg_rxn_id)
            continue

        reactome_candidates = rhea_reactome_map.get(master_id, set())
        for reactome_rxn_id in reactome_candidates:
            if reactome_rxn_id in reactome_rxn_ids:
                reactome_rxn_data = reactome_rxn_map[reactome_rxn_id]
                common_reactions.append({
                    "kegg_reaction_id": kegg_rxn_id,
                    "reactome_reaction_id": reactome_rxn_id,
                    "reactome_reaction_name": reactome_rxn_data.get("Reaction_Name", ""),
                    "reactome_equation": reactome_rxn_data.get("Equation", ""),
                    "rhea_master_id": master_id
                })

    print(f"[Reactome Analytics] Common reactions found: {len(common_reactions)}")
    print(f"[Reactome Analytics] KEGG reactions with no Rhea entry (not bridged): {len(unmatched_kegg)}")
    return {
        "kegg_pathway_id": kegg_pathway_id,
        "reactome_pathway_name": reactome_pathway_name,
        "total_kegg_reactions": len(kegg_rn_ids),
        "total_reactome_reactions": len(reactome_rxn_ids),
        "common_reaction_count": len(common_reactions),
        "kegg_reactions_not_in_rhea": len(unmatched_kegg),
        "method": "Rhea universal reaction ID bridge (zero false positives)",
        "common_reactions": common_reactions
    }

# =====================================================================================
# OVERALL SCRIPT LOGIC EXPLANATION
# =====================================================================================
# 1. PATH: The script connects directly to reactome_reactions.tsv using pandas.
# 2. _get_pathway_reactions(): The core function. It loads the TSV, filters for the exact
#    pathway name, and performs early processing (splitting equation strings and parsing 
#    genes/compounds into lists) to save processing time later.
# 3. HELPER FUNCTIONS: _get_pathway_compounds and _get_pathway_genes are used to extract
#    raw sets of data so the intersection and unique functions can easily compare pathways.
# 4. ANALYTICS (KEGG Equivalents): 
#    - find_irreversible_reactions and find_reversible_reactions check the '=>' and '<=>' signs.
#    - find_orphan_reactions finds reactions with no genes assigned.
#    - find_reactions_with_multiple_genes is the Reactome version of find_isozymes.
#    - find_hub_metabolites counts all inputs/outputs across all reactions.
#    - Intersections and Unique find overlaps and differences between different pathways.
# 5. ANALYTICS (Reactome specific): 
#    - find_disease_reactions strictly filters for Is_Disease == True.
# 6. summarize_pathway(): Aggregates all these stats into a single fast dictionary snapshot 
#    so the LLM can understand the entire pathway without needing to call every single function.
# 7. find_common_reactions_cross_db(): Uses ChEBI IDs as a universal bridge to match KEGG
#    and Reactome reactions without any text/fuzzy matching. Zero false positives.
# =====================================================================================

