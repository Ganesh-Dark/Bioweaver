import sys
sys.dont_write_bytecode = True
import os
import pickle
import requests
import re
import urllib.parse
from hgnc_client import expand_gene_symbol
from functools import lru_cache
import pandas as pd
import numpy as np
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity

# Global cache for the ML model and data
ML_MODEL = None
TFIDF_MATRIX = None
REACTOME_DATA = None

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
# Path where the trained model is saved so we don't retrain every time
PKL_PATH = os.path.join(BASE_DIR, "Reactome_files", "reactome_tfidf_model.pkl")

# Common cofactors, ions, and inorganic molecules to exclude from the metabolite list
COFACTORS_STOPLIST = {
    "ATP(4-)", "ADP(3-)", "AMP(1-)", "adenosine 5'-monophosphate(2-)",
    "NAD(1-)", "NADH(2-)", "NADP", "NADPH", "FAD", "FADH2",
    "water", "hydron", "hydrogenphosphate", "phosphate(2-)", "phosphate(4-)",
    "magnesium(2+)", "potassium(1+)", "ammonium", "calcium(2+)",
    "3',5'-cyclic AMP(1-)", "acetyl-CoA(4-)", "coenzyme A(4-)"
}

def _initialize_ml_model():
    global ML_MODEL, TFIDF_MATRIX, REACTOME_DATA
    if ML_MODEL is not None:
        return

    # The .pkl file must be built first using build_reactome_model.py
    if not os.path.exists(PKL_PATH):
        print(f"Reactome model .pkl not found at: {PKL_PATH}")
        print("Please run build_reactome_model.py first to generate the model file.")
        return

    print("Loading Reactome TF-IDF model from .pkl file...")
    with open(PKL_PATH, "rb") as f:
        cached = pickle.load(f)
    ML_MODEL = cached["model"]
    TFIDF_MATRIX = cached["matrix"]
    REACTOME_DATA = cached["data"]
    print(f"Loaded {len(REACTOME_DATA)} Reactome pathways from cache.")

@lru_cache(maxsize=32)
def search_reactome_pathway_id(pathway_name, top_k=3):
    """
    Searches Reactome using a TF-IDF ML model.
    If the top match is very confident, returns its ID.
    If there are multiple close matches (ambiguity), it returns a dict with options
    so the AI can ask the user to clarify.
    """
    print(f"Searching Reactome for pathway: '{pathway_name}'")
    if ML_MODEL is None:
        _initialize_ml_model()
        
    if ML_MODEL is None or REACTOME_DATA is None or REACTOME_DATA.empty:
        return None

    # Check if the query is a direct Reactome ID
    query_upper = str(pathway_name).strip().upper()
    
    if query_upper.startswith("R-HSA-"):
        if query_upper in REACTOME_DATA["id"].values:
            print(f"Direct ID match found: {query_upper}")
            return query_upper
            
    if re.fullmatch(r'\d+', query_upper):
        possible_id = f"R-HSA-{query_upper}"
        if possible_id in REACTOME_DATA["id"].values:
            print(f"Direct ID match found: {possible_id}")
            return possible_id

    try:
        # Preprocess query: remove generic words that KEGG adds (like " / Gluconeogenesis")
        clean_query = pathway_name.lower().split(" / ")[0].strip()
        
        query_vec = ML_MODEL.transform([clean_query])
        sim_scores = cosine_similarity(query_vec, TFIDF_MATRIX).flatten()
        top_indices = np.argsort(sim_scores)[::-1][:top_k+2] # Grab a bit more in case we filter
        
        results = []
        
        GENERIC_STOP_PATHWAYS = {"metabolism", "disease", "diseases", "signaling by rtks", "transport of small molecules", "metabolism of proteins"}
        
        for idx in top_indices:
            score = float(sim_scores[idx])
            name = REACTOME_DATA.iloc[idx]["name"]
            
            # Penalize ultra-generic pathways unless the user typed exactly that word
            if name.lower() in GENERIC_STOP_PATHWAYS and clean_query.lower() != name.lower():
                score -= 0.4
                
            if score > 0.15 and len(results) < top_k:
                results.append({
                    "id": REACTOME_DATA.iloc[idx]["id"],
                    "name": name,
                    "score": round(score, 3)
                })
        
        if not results:
            print(f"No Reactome pathway found for: '{pathway_name}'")
            return None
            
        # If the first match is extremely confident or it's the only match
        if len(results) == 1 or results[0]["score"] > 0.85:
            print(f"Found confident Reactome pathway: {results[0]['name']} ({results[0]['id']})")
            return results[0]['id']
            
        # If there is a big drop off between 1st and 2nd, take the 1st
        if len(results) > 1 and (results[0]["score"] - results[1]["score"] > 0.2):
            print(f"Found confident Reactome pathway: {results[0]['name']} ({results[0]['id']})")
            return results[0]['id']
            
        # AMBIGUITY DETECTED: Return options so the AI asks the user
        print(f"Ambiguity detected for Reactome query '{pathway_name}'. Returning options.")
        return {"ambiguous_options": results}

    except Exception as e:
        print(f"Error searching Reactome ML model: {e}")
        return None


def get_reactome_pathways_list(pathway_name, top_k=10):
    # pathway_name: search query (e.g. "amino acid metabolism")
    # top_k: max number of results to return
    # Returns a list of dicts with 'id' and 'name', bypassing the ambiguity check.
    if ML_MODEL is None:
        _initialize_ml_model()
    if ML_MODEL is None or REACTOME_DATA is None or REACTOME_DATA.empty:
        return []

    # Direct ID check
    query_upper = str(pathway_name).strip().upper()
    if query_upper.startswith("R-HSA-"):
        match = REACTOME_DATA[REACTOME_DATA["id"] == query_upper]
        if not match.empty:
            return [{"id": query_upper, "name": match.iloc[0]["name"]}]
    if re.fullmatch(r'\d+', query_upper):
        possible_id = f"R-HSA-{query_upper}"
        match = REACTOME_DATA[REACTOME_DATA["id"] == possible_id]
        if not match.empty:
            return [{"id": possible_id, "name": match.iloc[0]["name"]}]

    clean_query = str(pathway_name).lower().split(" / ")[0].strip()
    query_vec = ML_MODEL.transform([clean_query])
    sim_scores = cosine_similarity(query_vec, TFIDF_MATRIX).flatten()
    top_indices = np.argsort(sim_scores)[::-1][:top_k + 2]

    results = []
    GENERIC_STOP = {"metabolism", "disease", "diseases", "signaling by rtks", "transport of small molecules", "metabolism of proteins"}

    for idx in top_indices:
        score = float(sim_scores[idx])
        name = REACTOME_DATA.iloc[idx]["name"]
        if name.lower() in GENERIC_STOP and clean_query.lower() != name.lower():
            score -= 0.4
        if score > 0.15 and len(results) < top_k:
            results.append({"id": REACTOME_DATA.iloc[idx]["id"], "name": name})

    print(f"Reactome list search for '{pathway_name}': found {len(results)} results.")
    return results


@lru_cache(maxsize=32)
def get_reactome_pathway_data(pathway_name):
    # pathway_name: The name of the pathway to search for (e.g. "glycolysis")
    # Returns a dict with 'reactome_id', 'genes', and 'compounds'.
    reactome_id = search_reactome_pathway_id(pathway_name)
    if isinstance(reactome_id, dict) and "ambiguous_options" in reactome_id:
        return reactome_id
        
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
