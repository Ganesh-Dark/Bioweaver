import os
import time
import requests
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

# --- CONFIGURATION ---
MAX_WORKERS = 10 # Number of parallel downloads (keeps it fast without getting banned)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
REACTOME_DIR = os.path.join(BASE_DIR, "Reactome_files")
PATHWAYS_FILE = os.path.join(REACTOME_DIR, "ReactomePathways.txt")
OUTPUT_FILE = os.path.join(REACTOME_DIR, "reactome_reactions.tsv")

def fetch_json(url, retries=3):
    """
    Fetches JSON from a URL with automatic retries if the server is busy.
    """
    for attempt in range(retries):
        try:
            res = requests.get(url, timeout=15)
            if res.status_code == 200:
                return res.json()
            elif res.status_code == 429: # Rate limited by Reactome
                time.sleep(2)
        except Exception:
            time.sleep(1)
    return None

def get_reaction_details(reaction_id):
    """
    Fetches the deep details (11 fields) for a specific reaction.
    """
    query_url = f"https://reactome.org/ContentService/data/query/{reaction_id}"
    parts_url = f"https://reactome.org/ContentService/data/participants/{reaction_id}"
    
    # 1. Fetch metadata (Summary, Literature, Disease, Equation, Compartment)
    meta = fetch_json(query_url)
    if not meta:
        return None
        
    # 2. Fetch specific physical parts (Genes, Enzymes, Substrates)
    parts = fetch_json(parts_url)
    
    # --- PARSE THE 11 FIELDS ---
    
    # 1. Reaction ID
    r_id = meta.get("stId", reaction_id)
    
    # 2 & 3. Name and Equation
    names = meta.get("name", [])
    r_name = names[0] if len(names) > 0 else meta.get("displayName", "")
    r_eq = names[1] if len(names) > 1 else ""
    
    # 4 & 5 & 6 & 7. Inputs, Outputs, Catalysts, Genes
    inputs = []
    outputs = []
    catalysts = []
    genes = set()
    
    # Simple strings from the query endpoint
    for item in meta.get("input", []): 
        if isinstance(item, dict): inputs.append(item.get("displayName", ""))
    for item in meta.get("output", []): 
        if isinstance(item, dict): outputs.append(item.get("displayName", ""))
    for item in meta.get("catalystActivity", []): 
        if isinstance(item, dict) and "physicalEntity" in item:
            catalysts.append(item["physicalEntity"].get("displayName", ""))
            
    # Dig deep into participants to get pure gene symbols
    if parts:
        for p in parts:
            if "refEntities" in p:
                for ref in p["refEntities"]:
                    if ref.get("schemaClass") in ["ReferenceGeneProduct", "ReferenceIsoform", "ReferenceSequence"]:
                        gene_symbol = ref.get("displayName", "").split()[-1]
                        genes.add(gene_symbol)
                        
    # 8. Compartment
    compartments = [c.get("displayName", "") for c in meta.get("compartment", [])]
    
    # 9. Summation
    summation = ""
    summations = meta.get("summation", [])
    if len(summations) > 0:
        summation = summations[0].get("text", "")
        
    # 10. Literature References
    lit_refs = []
    for lit in meta.get("literatureReference", []):
        if "pubMedIdentifier" in lit:
            lit_refs.append(str(lit["pubMedIdentifier"]))
            
    # 11. Disease Flag
    is_disease = meta.get("isInDisease", False)
    
    return {
        "Reaction_ID": r_id,
        "Reaction_Name": r_name,
        "Equation": r_eq,
        "Inputs": "; ".join(inputs),
        "Outputs": "; ".join(outputs),
        "Catalysts": "; ".join(catalysts),
        "Genes": "; ".join(sorted(list(genes))),
        "Compartment": "; ".join(compartments),
        "Summation": summation,
        "PubMed_IDs": "; ".join(lit_refs),
        "Is_Disease": is_disease
    }

def process_pathway(pathway_id, pathway_name):
    """
    Downloads all reactions for a single pathway.
    """
    events_url = f"https://reactome.org/ContentService/data/pathway/{pathway_id}/containedEvents"
    events = fetch_json(events_url)
    if not events:
        return []
        
    results = []
    for event in events:
        # Event can be an integer (dbId) or a dictionary
        is_reaction = False
        event_id = None
        
        if isinstance(event, dict):
            if event.get("className") == "Reaction":
                is_reaction = True
                event_id = event.get("stId", event.get("dbId"))
        elif isinstance(event, int):
            # If it's just an int, we query it to see if it's a reaction
            meta = fetch_json(f"https://reactome.org/ContentService/data/query/{event}")
            if meta and meta.get("className") == "Reaction":
                is_reaction = True
                event_id = meta.get("stId", event)
                
        if is_reaction and event_id:
            details = get_reaction_details(event_id)
            if details:
                # Add the pathway info to the 11 fields
                details["Pathway_ID"] = pathway_id
                details["Pathway_Name"] = pathway_name
                results.append(details)
                
    return results

def main():
    print(f"Reading human pathways from {PATHWAYS_FILE}...")
    if not os.path.exists(PATHWAYS_FILE):
        print("Error: ReactomePathways.txt not found. Cannot proceed.")
        return
        
    pathways = []
    with open(PATHWAYS_FILE, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3 and parts[2] == "Homo sapiens":
                pathways.append((parts[0], parts[1]))
                
    # --- RESUME CAPABILITY ---
    processed_pathways = set()
    write_header = True
    if os.path.exists(OUTPUT_FILE):
        try:
            existing_df = pd.read_csv(OUTPUT_FILE, sep="\t", usecols=["Pathway_ID"])
            processed_pathways = set(existing_df["Pathway_ID"].unique())
            write_header = False
            print(f"Found {len(processed_pathways)} already processed pathways. Resuming...")
        except:
            pass

    pathways_to_process = [p for p in pathways if p[0] not in processed_pathways]
    total_pathways = len(pathways_to_process)
    
    if total_pathways == 0:
        print("All pathways are already downloaded!")
        return
        
    print(f"Starting parallel download for {total_pathways} pathways...")
    
    cols = ["Pathway_ID", "Pathway_Name", "Reaction_ID", "Reaction_Name", "Equation", 
            "Inputs", "Outputs", "Catalysts", "Genes", "Compartment", "Summation", 
            "PubMed_IDs", "Is_Disease"]
            
    processed = 0
    
    try:
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            future_to_pathway = {executor.submit(process_pathway, pid, pname): (pid, pname) for pid, pname in pathways_to_process}
            
            for future in as_completed(future_to_pathway):
                pid, pname = future_to_pathway[future]
                processed += 1
                try:
                    reactions = future.result()
                    if reactions:
                        df = pd.DataFrame(reactions)
                        # Ensure columns exist even if empty
                        for c in cols:
                            if c not in df.columns:
                                df[c] = ""
                        df = df[cols]
                        
                        # Append to file instantly!
                        df.to_csv(OUTPUT_FILE, sep="\t", index=False, header=write_header, mode="a", encoding="utf-8")
                        write_header = False
                        
                    print(f"[{processed}/{total_pathways}] Completed: {pname} ({len(reactions)} reactions found)")
                except Exception as exc:
                    print(f"[{processed}/{total_pathways}] Error processing {pname}: {exc}")
                    
    except KeyboardInterrupt:
        print("\nDownload safely halted by user! All progress has been saved directly to the file.")
        
    print(f"\nFinished! All data is securely saved in: {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
