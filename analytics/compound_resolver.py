import requests
import json
import re

def is_valid_synonym(name):
    """Filters out known database IDs, registry numbers, and vendor codes from chemical synonyms."""
    bad_prefixes = (
        "CHEBI:", "CHEMBL", "UNII-", "DTXSID", "DTXCID", "RefChem:", 
        "SCHEMBL", "BDBM", "NSC", "EINECS", "Tox21_", "NCGC", "MFCD", 
        "BRN", "AI3-", "FEMA No.", "HY-", "DB-", "orb", "AKOS", "EN300-", 
        "InChI=", "CAS-", "GTPL", "HMS", "CS-W", "FP", "Tox21"
    )
    if name.startswith(bad_prefixes):
        return False
        
    if re.match(r'^\d{2,7}-\d{2}-\d$', name):
        return False
        
    if re.match(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$', name):
        return False
        
    if re.match(r'^Q\d+$', name):
        return False
        
    if re.match(r'^C\d{5}$', name):
        return False
        
    if re.match(r'^[A-Z0-9]+$', name) and any(c.isdigit() for c in name):
        return False
        
    if re.match(r'^[a-z0-9]+$', name) and any(c.isdigit() for c in name) and len(name) < 10:
        return False

    return True

def _get_pubchem_synonyms(inchikey):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/synonyms/JSON"
    res = requests.get(url)
    if res.status_code == 200:
        try:
            return res.json()["InformationList"]["Information"][0]["Synonym"]
        except (KeyError, IndexError):
            pass
    return []

def _get_kegg_synonyms(kegg_id):
    url = f"https://rest.kegg.jp/get/cpd:{kegg_id}"
    res = requests.get(url)
    if res.status_code == 200:
        synonyms = []
        for line in res.text.split("\n"):
            if line.startswith("NAME"):
                names_str = line.replace("NAME", "").strip()
                synonyms.extend([n.strip().rstrip(';') for n in names_str.split(";") if n.strip()])
                break
        return synonyms
    return []

def _get_unichem_xrefs(inchikey):
    url = f"https://www.ebi.ac.uk/unichem/rest/inchikey/{inchikey}"
    res = requests.get(url)
    source_map = {"1": "ChEMBL", "2": "DrugBank", "6": "KEGG", "7": "ChEBI", "18": "HMDB", "22": "PubChem"}
    results = {}
    if res.status_code == 200:
        for entry in res.json():
            src_id = str(entry["src_id"])
            if src_id in source_map:
                db_name = source_map[src_id]
                if db_name not in results:
                    results[db_name] = []
                results[db_name].append(entry["src_compound_id"])
    return results

def resolve_chemical(name, target_db="all"):
    """
    Resolves a chemical name using PubChem and UniChem to cross-reference identifiers and clean synonyms.
    
    :param name: The name of the chemical to search (e.g., 'pyruvate')
    :param target_db: Specific database to return ('KEGG', 'ChEBI', 'all')
    :return: A dictionary of results.
    """
    pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/InChIKey/JSON"
    res = requests.get(pubchem_url)
    
    if res.status_code != 200:
        return {"error": f"Could not find '{name}' in PubChem. Server returned {res.status_code}"}
        
    inchikey = res.json()["PropertyTable"]["Properties"][0]["InChIKey"]
    
    parts = inchikey.split('-')
    if len(parts) < 3:
        return {"error": f"Invalid InChIKey format returned: {inchikey}"}
        
    skeleton = parts[0]
    stereo = parts[1]
    
    # We forcefully check the Neutral form (Acid) to catch KEGG and other specific subsets
    inchikey_neutral = f"{skeleton}-{stereo}-N"
    
    results_original = _get_unichem_xrefs(inchikey)
    results_neutral = _get_unichem_xrefs(inchikey_neutral)
    
    merged_results = {}
    for db in ["ChEMBL", "DrugBank", "KEGG", "ChEBI", "HMDB", "PubChem"]:
        ids = set(results_original.get(db, []) + results_neutral.get(db, []))
        if ids:
            merged_results[db] = list(ids)
            
    # Handle database targeting
    if target_db.lower() != "all":
        actual_db_name = next((db for db in merged_results.keys() if db.lower() == target_db.lower()), None)
        
        if not actual_db_name:
            return {
                "query": name,
                "target_db": target_db,
                "error": f"No IDs found for target database: {target_db}",
                "available_dbs": list(merged_results.keys())
            }
        
        return {
            "query": name,
            "target_db": actual_db_name,
            "ids": merged_results[actual_db_name]
        }
        
    # If 'all', we fetch synonyms
    pubchem_syns_orig = _get_pubchem_synonyms(inchikey)
    pubchem_syns_neut = _get_pubchem_synonyms(inchikey_neutral)
    
    raw_syns = list(dict.fromkeys(pubchem_syns_orig + pubchem_syns_neut))
    clean_syns = [syn for syn in raw_syns if is_valid_synonym(syn)]
    
    kegg_syns = []
    if "KEGG" in merged_results:
        for k_id in merged_results["KEGG"]:
            kegg_syns.extend(_get_kegg_synonyms(k_id))
    kegg_syns = list(dict.fromkeys(kegg_syns))
    
    return {
        "query": name,
        "inchikey_base": inchikey,
        "inchikey_skeleton": skeleton,
        "cross_references": merged_results,
        "kegg_synonyms": kegg_syns,
        "pubchem_synonyms": clean_syns
    }
