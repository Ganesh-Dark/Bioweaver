import sys
sys.dont_write_bytecode = True
import os
import requests
import json
from functools import lru_cache
from hgnc_client import HGNC_CACHE, _load_hgnc_data
import urllib.parse
import re

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
                
        # Deduplicate and limit list elements to 10 to save tokens
        details["names"] = list(set(details["names"]))[:10]
        details["reactions"] = list(set(details["reactions"]))[:10]
        details["enzymes"] = list(set(details["enzymes"]))[:10]
        details["pathways"] = list(set(details["pathways"]))[:10]
        details["modules"] = list(set(details["modules"]))[:10]
        
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
                    # KEGG descriptions format: "SYMBOL, ALIAS; FULL NAME"
                    raw_desc = parts[1].strip()
                    if ";" in raw_desc:
                        symbols_part, fullname_part = raw_desc.split(";", 1)
                        primary_symbol = symbols_part.split(",")[0].strip()
                        fallback_full_name = fullname_part.strip()
                    else:
                        primary_symbol = raw_desc.split(",")[0].strip()
                        fallback_full_name = None
                        
                    # Try to get official name from HGNC
                    if not HGNC_CACHE:
                        _load_hgnc_data()
                        
                    official_name = HGNC_CACHE.get(primary_symbol)
                    
                    if official_name:
                        translated.append(f"{primary_symbol} (Full Name: {official_name}) (KEGG ID: {k_id})")
                    elif fallback_full_name:
                        # Fallback to KEGG name if not found in HGNC
                        translated.append(f"{primary_symbol} (Full Name: {fallback_full_name}) (KEGG ID: {k_id})")
                    else:
                        translated.append(f"{raw_desc} (KEGG ID: {k_id})")
                        
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
                            m_copy = m.copy()
                            if "(Fallback Match)" not in m_copy.get("name", ""):
                                m_copy["name"] = m_copy.get("name", "") + " (Fallback Match)"
                            matches.append(m_copy)
                    if len(matches) > 0:
                        break
                        
        return matches
    except Exception as e:
        print(f"Error searching KEGG pathway: {e}")
        return []

def get_reactome_pathway_description(pathway_name):
    print(f"Fetching description for Reactome Pathway: {pathway_name}")
    import reactome_client
    reactome_id = reactome_client.search_reactome_pathway_id(pathway_name)
    
    if not reactome_id:
        return ""
        
    # If ML model found ambiguity, return the options dict so agent_graph can catch it
    if isinstance(reactome_id, dict) and "ambiguous_options" in reactome_id:
        return reactome_id
        
    url = f"https://reactome.org/ContentService/data/query/{reactome_id}"
    try:
        response = requests.get(url, timeout=15)
        data = response.json()
        if "summation" in data and len(data["summation"]) > 0:
            summary = data["summation"][0].get("text", "")
            summary = re.sub(r'<[^>]+>', '', summary)
            return summary.strip()
        return ""
    except Exception as e:
        print(f"Error fetching Reactome description for {reactome_id}: {e}")
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
                parts = content.split(";", 1)
                gene_symbols = parts[0].split()
                if len(gene_symbols) > 1:
                    gene_id = gene_symbols[0]
                    primary_symbol = gene_symbols[1].strip(',')
                    full_name = parts[1].strip() if len(parts) > 1 else ""
                    if full_name:
                        details["genes"].append(f"{primary_symbol} (Full Name: {full_name}) (ID: {gene_id})")
                    else:
                        details["genes"].append(f"{primary_symbol} (ID: {gene_id})")
            
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

def get_compound_synonyms(query: str) -> list:
    """
    Looks up a query in kegg_compounds_synonyms.tsv to find all its aliases.
    Returns a list of lowercase synonyms (including the original query).
    """
    query_lower = query.lower()
    synonyms = [query_lower]
    base_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(base_dir, "Kegg_files", "kegg_compounds_synonyms.tsv")
    
    if os.path.exists(file_path):
        try:
            with open(file_path, "r", encoding="utf-8") as f:
                for line in f:
                    if query_lower in line.lower():
                        parts = line.strip().split("\t", 1)
                        if len(parts) >= 2:
                            syns = [s.strip().lower() for s in parts[1].split(";")]
                            if query_lower in syns:
                                synonyms.extend(syns)
                                synonyms.append(parts[0].strip().lower())
                                break
        except Exception as e:
            print(f"Error reading synonyms file: {e}")
            
    return list(set(synonyms))

def search_local_kegg_tsv(query, filename):
    """
    Generic function to search a local TSV file.
    """
    # query: Search term
    # filename: The name of the TSV file (e.g., "ko.tsv")
    print(f"Searching {filename} for: '{query}'")
    
    # Smart synonym expansion
    synonyms_to_check = get_compound_synonyms(query)
    if len(synonyms_to_check) > 1:
        print(f"Expanded search to include synonyms: {synonyms_to_check}")
        
    matches = []
    base_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(base_dir, "Kegg_files", filename)
    
    if os.path.exists(file_path):
        try:
            with open(file_path, "r", encoding="utf-8") as f:
                for line in f:
                    # Pad with spaces so 'glucose' won't match 'UDP-glucose' (dash, not space, precedes glucose)
                    padded = " " + line.lower().replace("\t", " ") + " "
                    if any(f" {syn} " in padded for syn in synonyms_to_check):
                        parts = line.strip().split("\t", 1)
                        if len(parts) >= 2:
                            matches.append({"id": parts[0], "name": parts[1]})
            if matches:
                # Sort matches to prioritize exact starts with the original query
                matches.sort(key=lambda x: (
                    not x["name"].lower().startswith(query.lower()),
                    len(x["name"])
                ))
                print(f"Found {len(matches)} matches in local {filename}")
                return matches
        except Exception as e:
            print(f"Error reading local {filename}: {e}")
    else:
        print(f"File {filename} not found.")
    return []

@lru_cache(maxsize=32)
def get_generic_kegg_details(kegg_id):
    """
    Generic function to fetch details from KEGG API.
    """
    # kegg_id: The full KEGG ID including prefix (e.g. "ko:K00001")
    print(f"Fetching full details for: {kegg_id}")
    url = f"{KEGG_REST_BASE}/get/{kegg_id}"
    
    details = {}
    
    try:
        response = requests.get(url, timeout=10)
        if response.status_code != 200:
            return {"error": f"Failed to fetch details for {kegg_id}"}
            
        current_section = None
        for line in response.text.split('\n'):
            if not line:
                continue
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
            
            if current_section not in details:
                details[current_section] = []
            
            # Simple deduplication or accumulation
            details[current_section].append(content.strip(';'))
            
        # Clean up and filter
        keys_to_remove = ["DBLINKS", "COMMENT", "HISTORY", "REFERENCE"]
        for k in keys_to_remove:
            if k in details:
                del details[k]
                
        for key in details:
            details[key] = [item for item in details[key] if item][:10]
            
        return details
    except Exception as e:
        print(f"Error fetching details: {e}")
        return {"error": str(e)}

def parse_equation(equation_str):
    """
    Parses an equation string into substrates and products.
    """
    # equation_str: The string containing the reaction equation
    substrates = []
    products = []
    
    if "<=>" in equation_str:
        left, right = equation_str.split("<=>", 1)
    elif "=>" in equation_str:
        left, right = equation_str.split("=>", 1)
    else:
        left, right = equation_str, ""
        
    if left:
        substrates = [s.strip() for s in left.split("+") if s.strip()]
    if right:
        products = [p.strip() for p in right.split("+") if p.strip()]
        
    return substrates, products

def format_kegg_details(raw, default_id):
    """
    Formats the raw KEGG dictionary into a complete, lowercase dictionary.
    Automatically parses equations into substrates and products.
    """
    # raw: The raw dictionary directly from KEGG API
    # default_id: The ID to use as a fallback if ENTRY is missing
    if "error" in raw: return raw
    
    result = {k.lower(): v for k, v in raw.items()}
    result["id"] = raw.get("ENTRY", [""])[0].split()[0] if raw.get("ENTRY") else default_id
    
    if "EQUATION" in raw and raw["EQUATION"]:
        equation_str = raw["EQUATION"][0]
        result["equation"] = equation_str
        subs, prods = parse_equation(equation_str)
        result["substrates"] = subs
        result["products"] = prods
        
    # Some entries have REACTION but no EQUATION
    if "REACTION" in raw and "EQUATION" not in raw:
        rxn_str = " ".join(raw["REACTION"])
        if "<=>" in rxn_str or "=>" in rxn_str:
            subs, prods = parse_equation(rxn_str)
            if "substrates" not in result:
                result["substrates"] = []
            if "products" not in result:
                result["products"] = []
            result["substrates"].extend(subs)
            result["products"].extend(prods)
            
    return result

def filter_kegg_matches(matches, detail_fetcher, filter_field, filter_value, max_to_check=10, max_to_return=10):
    """
    Filters initial matches by fetching details and strictly verifying the filter.
    Stops after finding max_to_return matches or checking max_to_check items.
    """
    # matches: The basic matches list from TSV search
    # detail_fetcher: The function to fetch full details for a match
    # filter_field: The specific dictionary field to inspect
    # filter_value: The value that must exist in that field
    if not filter_field or not filter_value:
        return matches
        
    filtered = []
    synonyms = get_compound_synonyms(filter_value)
    
    for m in matches[:max_to_check]:
        details = detail_fetcher(m["id"])
        if "error" in details:
            continue
            
        field_data = details.get(filter_field.lower())
        if not field_data:
            continue
            
        found = False
        if isinstance(field_data, list):
            for item in field_data:
                padded = f" {str(item).lower()} "
                if any(f" {syn} " in padded for syn in synonyms):
                    found = True
                    break
        elif isinstance(field_data, str):
            padded = f" {field_data.lower()} "
            if any(f" {syn} " in padded for syn in synonyms):
                found = True
                
        if found:
            filtered.append(m)
            if len(filtered) >= max_to_return:
                break
            
    return filtered

def remove_stoichiometry(compound_str):
    compound_str = compound_str.strip()
    cleaned = re.sub(r'^(\d+|[a-zA-Z]|\([^\)]+\))\s+', '', compound_str)
    return cleaned.strip()

def apply_local_filter(data, detail_fetcher, filter_field, filter_value):
    if not filter_field or not filter_value:
        return True
        
    field_data = data.get(filter_field.lower())
    if field_data is None:
        details = detail_fetcher(data["id"])
        if "error" in details: return False
        field_data = details.get(filter_field.lower())
        if field_data is None: return False
        
    filter_values = [v.strip() for v in filter_value.split(",")]
    for fv in filter_values:
        fv_syns = get_compound_synonyms(fv)
        found_fv = False
        if isinstance(field_data, list):
            for item in field_data:
                if any(syn in str(item).lower() for syn in fv_syns):
                    found_fv = True
                    break
        elif isinstance(field_data, str):
            if any(syn in field_data.lower() for syn in fv_syns):
                found_fv = True
                
        if not found_fv:
            return False
            
    return True

# KO wrappers
def search_kegg_ko(query, filter_field=None, filter_value=None):
    """
    Searches for a KO locally and optionally filters.
    """
    matches = search_local_kegg_tsv(query, "ko.tsv")
    parsed_matches = []
    for m in matches:
        names = [n.strip() for n in m["name"].split("; ")]
        data = {"id": m["id"], "names": names}
        if apply_local_filter(data, get_kegg_ko_details, filter_field, filter_value):
            parsed_matches.append(data)
            if len(parsed_matches) >= 10: break
    return parsed_matches

def get_kegg_ko_details(ko_id):
    if not ko_id.startswith("ko:"): ko_id = f"ko:{ko_id}"
    raw = get_generic_kegg_details(ko_id)
    return format_kegg_details(raw, ko_id)

# Enzyme wrappers
def search_kegg_enzyme(query, filter_field=None, filter_value=None):
    """
    Searches for an enzyme locally and optionally filters.
    """
    matches = search_local_kegg_tsv(query, "enzyme.tsv")
    parsed_matches = []
    for m in matches:
        names = [n.strip() for n in m["name"].split("; ")]
        data = {"id": m["id"], "names": names}
        if apply_local_filter(data, get_kegg_enzyme_details, filter_field, filter_value):
            parsed_matches.append(data)
            if len(parsed_matches) >= 10: break
    return parsed_matches

def get_kegg_enzyme_details(ec_id):
    if not ec_id.startswith("ec:"): ec_id = f"ec:{ec_id}"
    raw = get_generic_kegg_details(ec_id)
    return format_kegg_details(raw, ec_id)

# Module wrappers
def search_kegg_module(query, filter_field=None, filter_value=None):
    """
    Searches for a Module locally and optionally filters.
    """
    matches = search_local_kegg_tsv(query, "module.tsv")
    parsed_matches = []
    for m in matches:
        data = {"id": m["id"], "name": m["name"]}
        if apply_local_filter(data, get_kegg_module_details, filter_field, filter_value):
            parsed_matches.append(data)
            if len(parsed_matches) >= 10: break
    return parsed_matches

def get_kegg_module_details(module_id):
    if not module_id.startswith("md:"): module_id = f"md:{module_id}"
    raw = get_generic_kegg_details(module_id)
    return format_kegg_details(raw, module_id)

# Reaction wrappers
def search_kegg_reaction(query, filter_field=None, filter_value=None):
    """
    Searches for a Reaction locally and optionally filters.
    """
    matches = search_local_kegg_tsv(query, "Reactions.tsv")
    parsed_matches = []
    for m in matches:
        rest = m["name"]
        equation = ""
        name_str = rest
        if "<=>" in rest or "=>" in rest:
            chunks = rest.split("; ")
            for i in range(len(chunks)-1, -1, -1):
                if "<=>" in chunks[i] or "=>" in chunks[i]:
                    equation = chunks[i]
                    name_str = "; ".join(chunks[:i])
                    break
                    
        substrates = []
        products = []
        is_reversible = None
        
        if equation:
            if "<=>" in equation:
                left, right = equation.split("<=>", 1)
                is_reversible = True
            elif "=>" in equation:
                left, right = equation.split("=>", 1)
                is_reversible = False
            else:
                left, right = equation, ""
                
            subs = [s.strip() for s in left.split(" + ") if s.strip()]
            prods = [p.strip() for p in right.split(" + ") if p.strip()]
            
            substrates = [remove_stoichiometry(s) for s in subs]
            products = [remove_stoichiometry(p) for p in prods]
            
        data = {
            "id": m["id"],
            "names": name_str.split("; "),
            "equation": equation,
            "substrates": substrates,
            "products": products,
            "reversible": is_reversible
        }
        
        if apply_local_filter(data, get_kegg_reaction_details, filter_field, filter_value):
            parsed_matches.append(data)
            if len(parsed_matches) >= 10: break
            
    return parsed_matches

def get_kegg_reaction_details(reaction_id):
    if not reaction_id.startswith("rn:"): reaction_id = f"rn:{reaction_id}"
    raw = get_generic_kegg_details(reaction_id)
    return format_kegg_details(raw, reaction_id)


def get_pathway_linked_data(pathway_id, link_type):
    """
    Calls the KEGG Link API to get IDs linked to a pathway, then cross-references
    them with local TSV files to return full details.
    """
    # pathway_id: A KEGG pathway ID such as map00010 or hsa00010.
    # link_type: One of "rn" (reactions), "ec" (enzymes), or "cpd" (compounds).

    # Normalise pathway_id: strip any "path:" prefix
    clean_id = pathway_id.replace("path:", "").strip()

    valid_types = ("rn", "ec", "cpd")
    if link_type not in valid_types:
        print(f"Invalid link_type '{link_type}'. Must be one of {valid_types}.")
        return {"error": f"link_type must be one of {valid_types}"}

    print(f"Calling KEGG Link API: /link/{link_type}/{clean_id}")
    url = f"https://rest.kegg.jp/link/{link_type}/{clean_id}"
    response = requests.get(url, timeout=15)

    if response.status_code != 200 or not response.text.strip():
        print(f"KEGG Link API returned no data for pathway: {clean_id}, type: {link_type}")
        return {"pathway_id": clean_id, "link_type": link_type, "results": []}

    # Parse the tab-delimited response to extract all IDs
    # Each line is: path:mapXXXXX  rn:RXXXXX  (or ec:X.X.X.X or cpd:CXXXXX)
    linked_ids = []
    for line in response.text.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) >= 2:
            # The second column is the linked ID (e.g. rn:R00299)
            raw_id = parts[1].strip()
            # Strip the prefix (rn:, ec:, cpd:) to get the bare ID
            bare_id = raw_id.split(":")[-1]
            linked_ids.append(bare_id)

    print(f"Found {len(linked_ids)} linked IDs for {clean_id} (type={link_type}). Looking up local data...")

    results = []
    base_dir = os.path.dirname(os.path.abspath(__file__))

    if link_type == "rn":
        # Cross-reference with local Reactions.tsv
        file_path = os.path.join(base_dir, "Kegg_files", "Reactions.tsv")
        if os.path.exists(file_path):
            # Build a set for fast lookup
            id_set = set(linked_ids)
            
            # Fetch EC mappings in bulk to prevent AI from having to loop!
            ec_mapping = {}
            try:
                print("Fetching EC mappings for reactions...")
                ec_rn_res = requests.get("https://rest.kegg.jp/link/ec/rn", timeout=10)
                if ec_rn_res.status_code == 200:
                    for line in ec_rn_res.text.strip().split("\n"):
                        if "\t" in line:
                            rn_part, ec_part = line.split("\t")
                            rn_id = rn_part.replace("rn:", "")
                            ec_id = ec_part.replace("ec:", "")
                            if rn_id in id_set:
                                ec_mapping.setdefault(rn_id, []).append(ec_id)
            except Exception as e:
                print("Failed to fetch EC mappings:", e)

            with open(file_path, "r", encoding="utf-8") as f:
                for line in f:
                    parts = line.strip().split("\t", 1)
                    if len(parts) == 2 and parts[0] in id_set:
                        rxn_id = parts[0]
                        rest = parts[1]

                        # Split name from equation
                        equation = ""
                        name_str = rest
                        if "<=>" in rest or "=>" in rest:
                            chunks = rest.split("; ")
                            for i in range(len(chunks) - 1, -1, -1):
                                if "<=>" in chunks[i] or "=>" in chunks[i]:
                                    equation = chunks[i]
                                    name_str = "; ".join(chunks[:i])
                                    break

                        results.append({
                            "id": rxn_id,
                            "name": name_str.split("; ")[0] if name_str else rxn_id,
                            "equation": equation,
                            "ec_numbers": ec_mapping.get(rxn_id, [])
                        })
        else:
            print("Reactions.tsv not found locally. Cannot cross-reference.")

    elif link_type == "ec":
        # Cross-reference with local enzyme.tsv
        file_path = os.path.join(base_dir, "Kegg_files", "enzyme.tsv")
        if os.path.exists(file_path):
            id_set = set(linked_ids)
            with open(file_path, "r", encoding="utf-8") as f:
                for line in f:
                    parts = line.strip().split("\t", 1)
                    if len(parts) == 2 and parts[0] in id_set:
                        # Take the first common name only (before the first semicolon)
                        common_name = parts[1].split(";")[0].strip()
                        results.append({
                            "ec": parts[0],
                            "name": common_name
                        })
        else:
            print("enzyme.tsv not found locally. Cannot cross-reference.")

    elif link_type == "cpd":
        # Cross-reference with local kegg_compounds_synonyms.tsv
        file_path = os.path.join(base_dir, "Kegg_files", "kegg_compounds_synonyms.tsv")
        if os.path.exists(file_path):
            id_set = set(linked_ids)
            with open(file_path, "r", encoding="utf-8") as f:
                for line in f:
                    parts = line.strip().split("\t", 1)
                    if len(parts) == 2 and parts[0] in id_set:
                        # Take the first listed synonym as the primary name
                        primary_name = parts[1].split(";")[0].strip()
                        results.append({
                            "id": parts[0],
                            "name": primary_name
                        })
        else:
            print("kegg_compounds_synonyms.tsv not found locally. Cannot cross-reference.")

    print(f"Returning {len(results)} locally matched results for pathway {clean_id}.")
    return {
        "pathway_id": clean_id,
        "link_type": link_type,
        "total_found": len(results),
        "results": results
    }

# Detailed Script Workings (new addition):
# get_pathway_linked_data fetches cross-reference IDs from the KEGG Link API
# (/link/rn, /link/ec, /link/cpd) and then looks up the bare IDs in the local
# Reactions.tsv, enzyme.tsv, or kegg_compounds_synonyms.tsv files respectively.
# This avoids multiple slow API calls for individual entries.
