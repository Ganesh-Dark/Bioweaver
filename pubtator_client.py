import sys
sys.dont_write_bytecode = True
import requests
import urllib.parse

def search_pubtator_papers(query):
    # query: A search string (e.g. "TP53 cancer")
    
    if not query:
        return []
        
    print(f"Fetching PubTator papers for: '{query}'")
    
    encoded_query = urllib.parse.quote(query)
    url = f"https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/?text={encoded_query}"
    
    try:
        response = requests.get(url, timeout=15)
        if response.status_code != 200:
            return []
            
        data = response.json()
        results = data.get("results", [])
        
        papers = []
        for r in results[:5]:
            pmid = str(r.get("pmid", "N/A"))
            title = r.get("title", "No Title")
            authors_list = r.get("authors", [])
            authors = ", ".join(authors_list) if isinstance(authors_list, list) and len(authors_list) > 0 else "Unknown Authors"
            journal = r.get("journal", "Unknown Journal")
            
            date_str = r.get("date", "")
            year = date_str[:4] if date_str and len(date_str) >= 4 else "Unknown Year"
            
            papers.append({
                "pmid": pmid,
                "title": title,
                "authors": authors,
                "journal": journal,
                "year": year
            })
            
        print(f"Found {len(papers)} papers for query.")
        return papers
        
    except Exception as e:
        print(f"Error searching PubTator3: {e}")
        return []

def fetch_pubtator_annotations(pmid):
    # pmid: The PubMed ID of the paper
    
    # Return empty lists if there's no valid PMID
    if not pmid or str(pmid).strip() == "N/A":
        return {"Genes": [], "Diseases": [], "Chemicals": []}

    url = f"https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson?pmids={pmid}"
    
    try:
        response = requests.get(url, timeout=10)
        
        if response.status_code != 200:
            return {"Genes": [], "Diseases": [], "Chemicals": []}
            
        data = response.json()
        
        genes = set()
        diseases = set()
        chemicals = set()
        species = set()
        cell_lines = set()
        variants = set()
        relations = set()
        
        # Parse the BioC JSON format from PubTator 3
        pubtator_data = data.get("PubTator3", [])
        if not pubtator_data:
            return {"Genes": [], "Diseases": [], "Chemicals": [], "Species": [], "CellLines": [], "Variants": [], "Relations": []}
            
        for doc in pubtator_data:
            # 1. Parse annotations
            for passage in doc.get("passages", []):
                for annotation in passage.get("annotations", []):
                    infons = annotation.get("infons", {})
                    ann_type = infons.get("type", "")
                    
                    # Try to get the normalized scientific name first, if not fall back to the raw text matched
                    name = infons.get("name") or annotation.get("text", "")
                    
                    if not name:
                        continue
                        
                    # Capitalize for neatness
                    name = name.capitalize()
                    
                    if ann_type == "Gene":
                        genes.add(name)
                    elif ann_type == "Disease":
                        diseases.add(name)
                    elif ann_type == "Chemical":
                        chemicals.add(name)
                    elif ann_type == "Species":
                        species.add(name)
                    elif ann_type == "CellLine":
                        cell_lines.add(name)
                    elif ann_type == "Variant" or ann_type == "Mutation":
                        variants.add(name)
                        
            # 2. Parse relations
            for relation in doc.get("relations", []):
                infons = relation.get("infons", {})
                rel_type = infons.get("type", "Associates_with")
                role1 = infons.get("role1", {})
                role2 = infons.get("role2", {})
                
                name1 = role1.get("name", "").capitalize()
                name2 = role2.get("name", "").capitalize()
                
                if name1 and name2:
                    relations.add(f"{name1} -> [{rel_type}] -> {name2}")
                        
        return {
            "Genes": list(genes)[:5], # limit to 5 per category to keep terminal output clean
            "Diseases": list(diseases)[:5],
            "Chemicals": list(chemicals)[:5],
            "Species": list(species)[:5],
            "CellLines": list(cell_lines)[:5],
            "Variants": list(variants)[:5],
            "Relations": list(relations)[:5]
        }
        
    except Exception as e:
        print(f"Error fetching PubTator3 data for {pmid}: {e}")
        return {"Genes": [], "Diseases": [], "Chemicals": []}
