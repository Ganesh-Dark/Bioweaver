import sys
sys.dont_write_bytecode = True
import os
import requests
import json
import urllib.parse
import re
from functools import lru_cache

sys.dont_write_bytecode = True

KEGG_REST_BASE = "https://rest.kegg.jp"

# We will cache the pathways and diseases lists so we don't fetch them multiple times
# if the agent searches repeatedly.

def convert_ncbi_to_kegg_id(ncbi_gene_id):
    # ncbi_gene_id: The numeric NCBI Gene ID string (e.g. 7157).

    if ncbi_gene_id is None:
        print("Cannot convert to KEGG ID because NCBI Gene ID is None.")
        return None

    print(f"Converting NCBI Gene ID {ncbi_gene_id} to KEGG Gene ID")

    url = f"https://rest.kegg.jp/conv/genes/ncbi-geneid:{ncbi_gene_id}"
    response = requests.get(url, timeout=10)

    text = response.text.strip()
    if len(text) == 0:
        print(f"No KEGG Gene ID mapping found for NCBI Gene ID: {ncbi_gene_id}")
        return None

    lines = text.split("\n")
    first_line = lines[0]
    parts = first_line.split("\t")

    if len(parts) >= 2:
        kegg_gene_id = parts[1]
        print(f"Found KEGG Gene ID: {kegg_gene_id}")
        return kegg_gene_id

    print(f"Could not parse KEGG mapping response for NCBI Gene ID: {ncbi_gene_id}")
    return None

def get_kegg_pathways(kegg_gene_id):
    # kegg_gene_id: The KEGG Gene ID string (e.g. hsa:7157).

    if kegg_gene_id is None:
        print("Cannot fetch pathways because KEGG Gene ID is None.")
        return []

    print(f"Fetching pathways for KEGG Gene ID: {kegg_gene_id}")

    url = f"https://rest.kegg.jp/link/pathway/{kegg_gene_id}"
    response = requests.get(url, timeout=10)

    text = response.text.strip()
    if len(text) == 0:
        print(f"No pathways found for KEGG Gene ID: {kegg_gene_id}")
        return []

    pathway_ids = []
    lines = text.split("\n")
    for line in lines:
        parts = line.split("\t")
        if len(parts) >= 2:
            pathway_id = parts[1].replace("path:", "")
            pathway_ids.append(pathway_id)

    print(f"Found {len(pathway_ids)} pathway IDs for {kegg_gene_id}")

    org_code = kegg_gene_id.split(":")[0]
    
    text = ""
    base_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(base_dir, "Kegg_files", "kegg_pathways.tsv")
    if os.path.exists(file_path):
        print("Reading KEGG pathways from local file...")
        with open(file_path, "r", encoding="utf-8") as f:
            text = f.read().strip()
    else:
        print("Local file not found, fetching from KEGG API...")
        list_url = "https://rest.kegg.jp/list/pathway"
        list_res = requests.get(list_url, timeout=10)
        text = list_res.text.strip()

    pathway_name_map = {}
    if len(text) > 0:
        for line in text.split("\n"):
            parts = line.split("\t")
            if len(parts) >= 2:
                p_code = parts[0].replace("path:", "").replace("map", "")
                p_name = parts[1]
                pathway_name_map[p_code] = p_name

    pathway_details = []
    for p_id in pathway_ids:
        clean_p_id = p_id.replace(org_code, "")
        name = pathway_name_map.get(clean_p_id, p_id)
        pathway_details.append({
            "id": p_id,
            "name": name
        })

    return pathway_details

def get_kegg_diseases(kegg_gene_id):
    # kegg_gene_id: The KEGG Gene ID string (e.g. hsa:7157).

    if kegg_gene_id is None:
        print("Cannot fetch diseases because KEGG Gene ID is None.")
        return []

    print(f"Fetching diseases directly linked to KEGG Gene ID: {kegg_gene_id}")

    url = f"https://rest.kegg.jp/link/disease/{kegg_gene_id}"
    response = requests.get(url, timeout=10)

    text = response.text.strip()
    if len(text) == 0:
        print(f"No direct diseases found for KEGG Gene ID: {kegg_gene_id}")
        return []

    disease_ids = []
    lines = text.split("\n")
    for line in lines:
        parts = line.split("\t")
        if len(parts) >= 2:
            disease_id = parts[1].replace("ds:", "")
            disease_ids.append(disease_id)

    print(f"Found {len(disease_ids)} disease IDs for {kegg_gene_id}")

    text = ""
    base_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(base_dir, "Kegg_files", "kegg_diseases.tsv")
    if os.path.exists(file_path):
        print("Reading KEGG diseases from local file...")
        with open(file_path, "r", encoding="utf-8") as f:
            text = f.read().strip()
    else:
        print("Local file not found, fetching from KEGG API...")
        list_url = "https://rest.kegg.jp/list/disease"
        list_res = requests.get(list_url, timeout=10)
        text = list_res.text.strip()

    disease_name_map = {}
    if len(text) > 0:
        for line in text.split("\n"):
            parts = line.split("\t")
            if len(parts) >= 2:
                d_code = parts[0].replace("ds:", "")
                d_name = parts[1]
                disease_name_map[d_code] = d_name

    disease_details = []
    for d_id in disease_ids:
        name = disease_name_map.get(d_id, d_id)
        disease_details.append({
            "id": d_id,
            "name": name
        })

    return disease_details

@lru_cache(maxsize=32)
def search_kegg_compound(compound_name):
    """
    Searches the KEGG COMPOUND database for the given query.
    Returns a list of dictionaries with 'id' and 'name'.
    """
    print(f"Searching KEGG Compound database for: '{compound_name}'")
    encoded_query = urllib.parse.quote(compound_name)
    url = f"{KEGG_REST_BASE}/find/compound/{encoded_query}"
    
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200 and response.text.strip():
            matches = []
            for line in response.text.strip().split('\n'):
                parts = line.split('\t')
                if len(parts) >= 2:
                    matches.append({"id": parts[0], "name": parts[1]})
                    
            # Sort so exact matches or compounds starting with the query appear first (e.g. D-Glucose over UDP-glucose)
            query_lower = compound_name.lower()
            matches.sort(key=lambda x: (
                query_lower not in [n.strip().lower() for n in x["name"].split(';')],
                not x["name"].lower().startswith(query_lower),
                len(x["name"])
            ))
            
            print(f"Found {len(matches)} compound matches. Top match: {matches[0]['name']}")
            return matches
        else:
            print(f"No compound matches found for '{compound_name}' in KEGG API.")
            return []
    except Exception as e:
        print(f"Error searching KEGG compounds: {e}")
        return []

@lru_cache(maxsize=32)
def get_kegg_compound_details(compound_id):
    """
    Fetches details of a KEGG compound.
    """
    if not compound_id.startswith("cpd:"):
        compound_id = f"cpd:{compound_id}"
        
    print(f"Fetching full details for KEGG Compound: {compound_id}")
    url = f"{KEGG_REST_BASE}/get/{compound_id}"
    
    details = {
        "id": compound_id,
        "names": [],
        "formula": "",
        "exact_mass": "",
        "mol_weight": "",
        "pathways": [],
        "reactions": [],
        "modules": [],
        "enzymes": []
    }
    
    try:
        response = requests.get(url, timeout=10)
        if response.status_code != 200:
            return details
            
        current_section = None
        for line in response.text.split('\n'):
            if line.startswith(' '):
                # Continuation of previous section
                content = line.strip()
            else:
                # New section
                parts = line.split(maxsplit=1)
                if not parts: continue
                current_section = parts[0]
                content = parts[1] if len(parts) > 1 else ""
                
            if not content: continue
            
            if current_section == "NAME":
                # Remove trailing semicolons
                details["names"].append(content.strip(';'))
            elif current_section == "FORMULA":
                details["formula"] = content
            elif current_section == "EXACT_MASS":
                details["exact_mass"] = content
            elif current_section == "MOL_WEIGHT":
                details["mol_weight"] = content
            elif current_section == "PATHWAY":
                details["pathways"].append(content)
            elif current_section == "MODULE":
                details["modules"].append(content)
            elif current_section == "ENZYME":
                details["enzymes"].extend(content.split())
            elif current_section == "REACTION":
                # Multiple reactions can be space-separated
                details["reactions"].extend(content.split())
                
        # Deduplicate list elements
        details["names"] = list(set(details["names"]))
        details["reactions"] = list(set(details["reactions"]))
        details["enzymes"] = list(set(details["enzymes"]))
        
        return details
    except Exception as e:
        print(f"Error fetching KEGG compound details: {e}")
        return details


@lru_cache(maxsize=32)
def search_kegg_disease(disease_name):
    """
    Searches the KEGG DISEASE database for the given query.
    Returns a list of dictionaries with 'id' and 'name'.
    """
    # disease_name: The disease query string to search in KEGG (e.g. breast cancer).

    print(f"Searching KEGG Disease database for: '{disease_name}'")
    disease_matches = []
    base_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(base_dir, "Kegg_files", "kegg_diseases.tsv")
    if os.path.exists(file_path):
        print("Searching KEGG diseases from local file...")
        with open(file_path, "r", encoding="utf-8") as f:
            lines = f.readlines()
        
        for line in lines:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                d_id = parts[0].replace("ds:", "")
                d_name = parts[1]
                if disease_name.lower() in d_name.lower():
                    disease_matches.append({
                        "id": d_id,
                        "name": d_name
                    })
                    
        if len(disease_matches) == 0:
            clean_name = re.sub(r"'s\b", "", disease_name, flags=re.IGNORECASE).replace("'", "").replace('"', '')
            words = [w.lower() for w in clean_name.split() if w.lower() not in ["related", "diseases", "disease", "things", "problems", "in", "the", "and", "human", "humans"]]
            for word in words:
                if len(word) > 2:
                    print(f"Fallback searching KEGG Disease from local file for keyword: '{word}'")
                    for line in lines:
                        parts = line.strip().split("\t")
                        if len(parts) >= 2:
                            d_id = parts[0].replace("ds:", "")
                            d_name = parts[1]
                            if word in d_name.lower():
                                if not any(m["id"] == d_id for m in disease_matches):
                                    disease_matches.append({
                                        "id": d_id,
                                        "name": d_name
                                    })
    else:
        print("Local file not found, fetching from KEGG API...")
        url = f"https://rest.kegg.jp/find/disease/{disease_name.replace(' ', '+')}"
        response = requests.get(url, timeout=10)
        text = response.text.strip()
    
        if len(text) > 0:
            lines = text.split("\n")
            for line in lines:
                parts = line.split("\t")
                if len(parts) >= 2:
                    d_id = parts[0].replace("ds:", "")
                    d_name = parts[1]
                    disease_matches.append({
                        "id": d_id,
                        "name": d_name
                    })
    
        if len(disease_matches) == 0:
            clean_name = re.sub(r"'s\b", "", disease_name, flags=re.IGNORECASE).replace("'", "").replace('"', '')
            words = [w.lower() for w in clean_name.split() if w.lower() not in ["related", "diseases", "disease", "things", "problems", "in", "the", "and", "human", "humans"]]
            for word in words:
                if len(word) > 2:
                    print(f"Fallback searching KEGG Disease for keyword: '{word}'")
                    fb_url = f"https://rest.kegg.jp/find/disease/{word}"
                    fb_res = requests.get(fb_url, timeout=10)
                    fb_text = fb_res.text.strip()
                    if len(fb_text) > 0:
                        for line in fb_text.split("\n"):
                            parts = line.split("\t")
                            if len(parts) >= 2:
                                d_id = parts[0].replace("ds:", "")
                                d_name = parts[1]
                                if not any(m["id"] == d_id for m in disease_matches):
                                    disease_matches.append({
                                        "id": d_id,
                                        "name": d_name
                                    })

    print(f"Found {len(disease_matches)} matching KEGG Disease entries")
    return disease_matches

KEGG_ORGANISM_CACHE = {}

def get_kegg_organism_code(organism_name):
    # organism_name: The organism name or common name string (e.g. mouse, human, rat, dog).

    if not organism_name:
        return "hsa"

    clean_org = organism_name.lower().strip()

    if clean_org in KEGG_ORGANISM_CACHE:
        return KEGG_ORGANISM_CACHE[clean_org]

    print(f"Resolving KEGG organism code automatically for: '{organism_name}'")

    text = ""
    base_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(base_dir, "Kegg_files", "kegg_organisms.tsv")
    
    if os.path.exists(file_path):
        print("Reading KEGG organisms from local file...")
        with open(file_path, "r", encoding="utf-8") as f:
            text = f.read().strip()
    else:
        print("Local file not found, fetching from KEGG API...")
        url = "https://rest.kegg.jp/list/genome"
        res = requests.get(url, timeout=10)
        text = res.text.strip()

    if len(text) > 0:
        lines = text.split("\n")
        for line in lines:
            parts = line.split("\t")
            if len(parts) >= 2:
                info = parts[1]
                if ";" in info:
                    code, full_name = info.split(";", 1)
                    code = code.strip()
                    full_name = full_name.lower().strip()
    
                    if clean_org in full_name or full_name in clean_org:
                        print(f"Automatically resolved KEGG organism code: '{code}' for '{organism_name}'")
                        KEGG_ORGANISM_CACHE[clean_org] = code
                        return code

    if "mouse" in clean_org or "musculus" in clean_org:
        KEGG_ORGANISM_CACHE[clean_org] = "mmu"
        return "mmu"
    elif "rat" in clean_org:
        KEGG_ORGANISM_CACHE[clean_org] = "rno"
        return "rno"
    elif "human" in clean_org:
        KEGG_ORGANISM_CACHE[clean_org] = "hsa"
        return "hsa"

    print(f"Could not resolve KEGG organism code for: '{organism_name}'")
    KEGG_ORGANISM_CACHE[clean_org] = None
    return None

def get_disease_linked_genes_and_pathways(disease_id, organism="human"):
    # disease_id: The KEGG Disease entry ID (e.g. H00031).
    # organism: The organism name string to map KEGG organism code (e.g. mouse or human).

    print(f"Fetching linked genes and pathways for KEGG Disease ID: {disease_id}, organism: {organism}")

    org_code = get_kegg_organism_code(organism)
    if org_code is None:
        print(f"Cannot fetch linked genes because organism code is None for: '{organism}'")
        return {"gene_ids": [], "pathway_ids": []}

    target_urls = [
        f"https://rest.kegg.jp/link/{org_code}/{disease_id}",
        f"https://rest.kegg.jp/link/{org_code}/ds:{disease_id}",
        f"https://rest.kegg.jp/link/hsa/{disease_id}",
        f"https://rest.kegg.jp/link/hsa/ds:{disease_id}",
        f"https://rest.kegg.jp/link/genes/{disease_id}",
        f"https://rest.kegg.jp/link/genes/ds:{disease_id}"
    ]

    gene_ids = []
    for g_url in target_urls:
        g_res = requests.get(g_url, timeout=10)
        text = g_res.text.strip()
        if len(text) > 0:
            for line in text.split("\n"):
                parts = line.split("\t")
                if len(parts) >= 2:
                    target_id = parts[1]
                    if target_id not in gene_ids:
                        gene_ids.append(target_id)
            if len(gene_ids) > 0:
                break

    path_urls = [
        f"https://rest.kegg.jp/link/pathway/{disease_id}",
        f"https://rest.kegg.jp/link/pathway/ds:{disease_id}"
    ]

    pathway_ids = []
    for p_url in path_urls:
        p_res = requests.get(p_url, timeout=10)
        text = p_res.text.strip()
        if len(text) > 0:
            for line in text.split("\n"):
                parts = line.split("\t")
                if len(parts) >= 2:
                    pid = parts[1].replace("path:", "")
                    if pid not in pathway_ids:
                        pathway_ids.append(pid)
            if len(pathway_ids) > 0:
                break

    return {
        "gene_ids": gene_ids,
        "pathway_ids": pathway_ids
    }

def translate_kegg_genes_bulk(gene_ids):
    """
    Takes a list of KEGG gene IDs (e.g. ['hsa:351', 'hsa:348'])
    and bulk translates them to readable names (e.g. ['APP (hsa:351)', ...]).
    """
    if not gene_ids:
        return []
        
    query = "+".join(gene_ids)
    url = f"https://rest.kegg.jp/list/{query}"
    
    translated = []
    try:
        res = requests.get(url, timeout=10)
        text = res.text.strip()
        if text:
            for line in text.split("\n"):
                parts = line.split("\t")
                if len(parts) >= 2:
                    k_id = parts[0]
                    # Extract the primary gene symbol
                    desc = parts[1]
                    primary_name = desc.replace(";", ",").split(",")[0].strip()
                    translated.append(f"{primary_name} ({k_id})")
                    
        return translated if translated else gene_ids
    except Exception as e:
        print(f"Error bulk translating KEGG genes: {e}")
        return gene_ids

# Detailed Script Workings:
# 1. convert_ncbi_to_kegg_id sends a GET request to KEGG conv API endpoint to translate numeric NCBI Gene ID into KEGG Gene ID format (such as hsa:7157).
# 2. get_kegg_pathways queries the KEGG link pathway endpoint to get pathway IDs and maps their names in bulk using list/pathway/<org>.
# 3. get_kegg_diseases queries the KEGG link disease endpoint to extract disease IDs and maps their names in bulk using list/disease.
# 4. search_kegg_disease searches KEGG for matching disease entries by name.
# 5. get_disease_linked_genes_and_pathways fetches genes and pathway IDs associated with a specific KEGG disease ID.

@lru_cache(maxsize=32)
def search_kegg_pathway(pathway_name):
    """
    Searches the KEGG PATHWAY database for the given query.
    Returns a list of dictionaries with 'id' and 'name'.
    """
    # pathway_name: The name of the pathway to search for (e.g. mTOR signaling pathway)
    print(f"Searching KEGG Pathway database for: '{pathway_name}'")
    
    def _do_search(query):
        matches = []
        query_lower = query.lower()
        base_dir = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(base_dir, "Kegg_files", "Kegg_Pathways_All_Levels.tsv")
        
        if os.path.exists(file_path):
            try:
                with open(file_path, "r", encoding="utf-8") as f:
                    lines = f.readlines()
                # Skip header
                for line in lines[1:]:
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        pid = parts[0]
                        level1 = parts[1]
                        level2 = parts[2]
                        pname = parts[3]
                        if query_lower in pname.lower() or query_lower in level1.lower() or query_lower in level2.lower():
                            matches.append({
                                "id": pid, 
                                "name": pname,
                                "category": f"{level1} > {level2}"
                            })
                if matches:
                    print("Found pathway matches in local Kegg_Pathways_All_Levels.tsv file")
                    return matches
            except Exception as e:
                print(f"Error reading local KEGG pathways file: {e}")
                
        # Fallback to API if file doesn't exist or no matches found in file
        print("Local match not found or file missing, fetching from KEGG API...")
        try:
            url = f"https://rest.kegg.jp/find/pathway/{urllib.parse.quote(query)}"
            response = requests.get(url, timeout=30)
            text = response.text.strip()
            if text:
                for line in text.split("\n"):
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        pid = parts[0]
                        pname = parts[1]
                        if not any(exist["id"] == pid for exist in matches):
                            matches.append({"id": pid, "name": pname})
        except Exception as e:
            print(f"Error fetching from KEGG API: {e}")
            
        return matches

    try:
        matches = _do_search(pathway_name)
        
        if len(matches) == 0 and " " in pathway_name:
            words = [w.lower() for w in pathway_name.split() if w.lower() not in ["pathway", "pathways", "signaling", "metabolism", "degradation", "biosynthesis", "related", "and", "the", "in", "of"]]
            for word in words:
                if len(word) > 2:
                    print(f"Fallback searching KEGG Pathway for keyword: '{word}'")
                    fb_matches = _do_search(word)
                    for m in fb_matches:
                        if not any(exist["id"] == m["id"] for exist in matches):
                            matches.append(m)
                    if len(matches) > 0:
                        break
                        
        return matches
    except Exception as e:
        print(f"Error searching KEGG pathway: {e}")
        return []

def get_reactome_pathway_description(pathway_name):
    # pathway_name: The name of the pathway to search in Reactome
    print(f"Fetching description for Reactome Pathway: {pathway_name}")
    url = f"https://reactome.org/ContentService/search/query?query={urllib.parse.quote(pathway_name)}&species=Homo%20sapiens"
    
    try:
        response = requests.get(url, timeout=15)
        data = response.json()
        
        if "results" in data and len(data["results"]) > 0:
            entries = data["results"][0].get("entries", [])
            if entries:
                # Robust core-word scoring algorithm
                q_clean = re.sub(r'[^a-z0-9\s]', ' ', pathway_name.lower())
                q_words = set(q_clean.split())
                stop_words = {'pathway', 'signaling', 'cycle', 'metabolism', 'of', 'by', 'and', 'the', 'in', 'to', 'a'}
                q_core = q_words - stop_words
                if not q_core: 
                    q_core = q_words
                    
                best_entry = None
                best_score = float('inf')
                
                for entry in entries:
                    if entry.get("exactType") == "Pathway" and entry.get("summation"):
                        c_name = re.sub(r'<[^>]+>', '', entry.get("name", ""))
                        c_clean = re.sub(r'[^a-z0-9\s]', ' ', c_name.lower())
                        c_words = set(c_clean.split())
                        c_core = c_words - stop_words
                        if not c_core:
                            c_core = c_words
                            
                        # Calculate Jaccard-like distance for core words
                        intersection = q_core.intersection(c_core)
                        if not intersection:
                            continue # No core words match
                            
                        score = len(c_core) - len(intersection)
                        if score < best_score:
                            best_score = score
                            best_entry = entry
                            
                if best_entry:
                    summary = best_entry.get("summation")
                    # Strip Reactome's weird HTML highlighting tags
                    summary = re.sub(r'<[^>]+>', '', summary)
                    return summary.strip()
        return ""
    except Exception as e:
        print(f"Error fetching Reactome description: {e}")
        return ""

def get_kegg_pathway_description(pathway_id):
    # pathway_id: The KEGG pathway ID (e.g. map04150)
    print(f"Fetching description for KEGG Pathway: {pathway_id}")
    url = f"https://rest.kegg.jp/get/{pathway_id}"
    
    try:
        response = requests.get(url, timeout=30)
        text = response.text
        
        description = ""
        in_description = False
        
        for line in text.split("\n"):
            if line.startswith("DESCRIPTION"):
                in_description = True
                description += line.replace("DESCRIPTION", "").strip() + " "
            elif in_description:
                # If we hit a new all-caps keyword block at the start of a line (e.g., CLASS, PATHWAY_MAP), stop collecting
                if re.match(r'^[A-Z_]+\s', line):
                    break
                description += line.strip() + " "
                
        return description.strip()
    except Exception as e:
        print(f"Error fetching KEGG pathway description: {e}")
        return ""

@lru_cache(maxsize=32)
def get_kegg_pathway_details(pathway_id, organism="hsa"):
    """
    Fetches the full details of a KEGG pathway.
    """
    # pathway_id: The KEGG pathway ID (e.g. map04340)
    
    if pathway_id.startswith("map"):
        pathway_id = pathway_id.replace("map", organism)
    elif pathway_id.startswith("path:map"):
        pathway_id = pathway_id.replace("path:map", organism)
        
    print(f"Fetching full details for KEGG Pathway: {pathway_id}")
    url = f"https://rest.kegg.jp/get/{pathway_id}"
    
    details = {
        "description": "",
        "genes": [],
        "compounds": [],
        "drugs": [],
        "modules": [],
        "rel_pathways": []
    }
    
    try:
        response = requests.get(url, timeout=15)
        text = response.text
        
        current_section = None
        
        for line in text.split("\n"):
            if not line.startswith(" "):
                if line.strip() == "":
                    continue
                parts = line.split(maxsplit=1)
                current_section = parts[0]
                content = parts[1] if len(parts) > 1 else ""
            else:
                content = line.strip()
                
            if not content:
                continue
                
            if current_section == "DESCRIPTION":
                details["description"] += content + " "
                
            elif current_section == "GENE":
                parts = content.split(";")
                gene_symbols = parts[0].split()
                if len(gene_symbols) > 1:
                    details["genes"].append(f"{gene_symbols[1].strip(',')} (ID: {gene_symbols[0]})")
            
            elif current_section == "COMPOUND":
                parts = content.split(maxsplit=1)
                if len(parts) > 1:
                    details["compounds"].append(f"{parts[1].strip()} (ID: {parts[0]})")
                    
            elif current_section == "DRUG":
                parts = content.split(maxsplit=1)
                if len(parts) > 1:
                    details["drugs"].append(f"{parts[1].strip()} (ID: {parts[0]})")
                    
            elif current_section == "MODULE":
                parts = content.split(maxsplit=1)
                if len(parts) > 1:
                    details["modules"].append(f"{parts[1].strip()} (ID: {parts[0]})")
                    
            elif current_section == "REL_PATHWAY":
                parts = content.split(maxsplit=1)
                if len(parts) > 1:
                    details["rel_pathways"].append(f"{parts[1].strip()} (ID: {parts[0]})")
                    
        #for key in details:
        #    details[key] = details[key][:10]
            
        return details
    except Exception as e:
        print(f"Error fetching pathway details: {e}")
        return details

