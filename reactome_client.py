import sys
sys.dont_write_bytecode = True
import requests
import re
import urllib.parse
from hgnc_client import expand_gene_symbol
from functools import lru_cache

sys.dont_write_bytecode = True

# Common cofactors, ions, and inorganic molecules to exclude from the metabolite list
COFACTORS_STOPLIST = {
    "ATP(4-)", "ADP(3-)", "AMP(1-)", "adenosine 5'-monophosphate(2-)",
    "NAD(1-)", "NADH(2-)", "NADP", "NADPH", "FAD", "FADH2",
    "water", "hydron", "hydrogenphosphate", "phosphate(2-)", "phosphate(4-)",
    "magnesium(2+)", "potassium(1+)", "ammonium", "calcium(2+)",
    "3',5'-cyclic AMP(1-)", "acetyl-CoA(4-)", "coenzyme A(4-)"
}


@lru_cache(maxsize=32)
def search_reactome_pathway_id(pathway_name):
    # pathway_name: The name of the pathway to search for (e.g. "glycolysis")
    print(f"Searching Reactome for pathway: '{pathway_name}'")
    try:
        import os
        file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Reactome_files", "ReactomePathways.txt")
        candidates = []
        
        if os.path.exists(file_path):
            with open(file_path, 'r', encoding='utf-8') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        r_id, r_name, species = parts[0], parts[1], parts[2]
                        if species == "Homo sapiens":
                            candidates.append({"id": r_id, "name": r_name})
                            
        if not candidates:
            print(f"No local Homo sapiens Reactome pathways found for: {pathway_name}")
            return None

        # Robust core-word scoring algorithm
        # Removes punctuation and generic biological filler words to focus on the unique nouns
        q_clean = re.sub(r'[^a-z0-9\s]', ' ', pathway_name.lower())
        q_words = set(q_clean.split())
        stop_words = {'pathway', 'signaling', 'cycle', 'metabolism', 'of', 'by', 'and', 'the', 'in', 'to', 'a', 'disease', 'diseases'}
        q_core = q_words - stop_words
        if not q_core: 
            q_core = q_words
            
        def score_candidate(c):
            c_clean = re.sub(r'[^a-z0-9\s]', ' ', c["name"].lower())
            c_words = set(c_clean.split())
            c_core = c_words - stop_words
            if not c_core: 
                c_core = c_words
                
            # Exact match is always 0 penalty
            if c_clean.strip() == q_clean.strip():
                return 0
                
            # Penalty based on missing query words (heavy) + extra candidate words (light)
            missing = len(q_core - c_core)
            extra = len(c_core - q_core)
            return (missing * 100) + (extra * 10)

        candidates.sort(key=score_candidate)
        
        best = candidates[0]
        best_score = score_candidate(best)
        
        # Prevent hallucinated matches: if the score is too high (terrible match), safely reject it
        if best_score > 90:
            print(f"No sufficiently close Reactome pathway found for: '{pathway_name}'. (Closest was '{best['name']}' but it failed the similarity threshold).")
            return None
            
        print(f"Found Reactome pathway: {best['name']} ({best['id']})")
        return best['id']
    except Exception as e:
        print(f"Error searching Reactome: {e}")
        return None


@lru_cache(maxsize=32)
def get_reactome_pathway_data(pathway_name):
    # pathway_name: The name of the pathway to search for (e.g. "glycolysis")
    # Returns a dict with 'reactome_id', 'genes', and 'compounds'.
    reactome_id = search_reactome_pathway_id(pathway_name)
    if not reactome_id:
        print(f"No Reactome ID found for '{pathway_name}', skipping Reactome data.")
        return {"reactome_id": None, "genes": [], "compounds": []}

    print(f"Fetching Reactome participants for: {reactome_id}")
    url = f"https://reactome.org/ContentService/data/participants/{reactome_id}"

    try:
        response = requests.get(url, timeout=15)
        data = response.json()

        genes = set()
        compounds = set()

        for entity in data:
            if "refEntities" not in entity:
                continue
            for ref in entity["refEntities"]:
                schema = ref.get("schemaClass", "")

                if schema == "ReferenceGeneProduct":
                    # displayName is like "UniProt:P30153 PPP2R1A"
                    display = ref.get("displayName", "")
                    parts = display.split()
                    gene_name = parts[-1] if len(parts) > 1 else display
                    genes.add(gene_name)

                elif schema == "ReferenceMolecule":
                    # displayName is like "pyruvate [ChEBI:15361]"
                    display = ref.get("displayName", "")
                    clean_name = display.split(" [")[0].strip()
                    # Filter out cofactors and inorganic molecules
                    if clean_name not in COFACTORS_STOPLIST:
                        compounds.add(display)

        # Expand gene names instantly using local HGNC file
        expanded_genes = []
        for g in genes:
            expanded_genes.append(expand_gene_symbol(g))

        genes = sorted(expanded_genes)
        compounds = sorted(list(compounds))

        print(f"Reactome: Found {len(genes)} genes and {len(compounds)} metabolites.")
        return {
            "reactome_id": reactome_id,
            "genes": genes,
            "compounds": compounds
        }

    except Exception as e:
        print(f"Error fetching Reactome participants for {reactome_id}: {e}")
        return {"reactome_id": reactome_id, "genes": [], "compounds": []}


# How this script works:
# 1. search_reactome_pathway_id: Hits the Reactome ContentService search API to find the
#    correct Reactome stable ID (e.g. R-HSA-70171) for a given pathway name.
#    It filters strictly for Homo sapiens Pathways only.
#
# 2. get_reactome_pathway_data: Uses the stable ID from step 1 to call the
#    participants API. It separates the results into genes (ReferenceGeneProduct)
#    and compounds (ReferenceMolecule). Cofactors and inorganic ions are filtered
#    out using COFACTORS_STOPLIST so only true metabolites are returned.
#
# 3. The returned dict with 'reactome_id', 'genes', and 'compounds' is designed
#    to be directly used by the agent_graph.py pathway tool.
