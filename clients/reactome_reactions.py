import requests
from functools import lru_cache

@lru_cache(maxsize=32)
def get_reactome_pathway_reactions(pathway_id):
    """
    Fetches all biochemical reactions contained within a specific Reactome pathway.
    """
    print(f"Fetching reactions for Reactome Pathway: {pathway_id}")
    url = f"https://reactome.org/ContentService/data/pathway/{pathway_id}/containedEvents"
    
    try:
        response = requests.get(url, timeout=15)
        if response.status_code != 200:
            print(f"Failed to fetch reactions for {pathway_id} (Status: {response.status_code})")
            return []
            
        data = response.json()
        reactions = []
        
        for event_item in data:
            if isinstance(event_item, int):
                # Reactome sometimes returns raw integer IDs instead of full objects
                try:
                    q_res = requests.get(f"https://reactome.org/ContentService/data/query/{event_item}", timeout=10)
                    if q_res.status_code == 200:
                        event = q_res.json()
                    else:
                        continue
                except:
                    continue
            elif isinstance(event_item, dict):
                event = event_item
            else:
                continue
                
            # We only want actual biochemical reactions, not sub-pathways or black boxes
            if event.get("className") == "Reaction":
                names = event.get("name", [])
                
                # The first name is usually the human-readable description
                primary_name = names[0] if len(names) > 0 else "Unknown Reaction"
                
                # The second name (if it exists) is usually the chemical equation
                equation = names[1] if len(names) > 1 else ""
                
                reactions.append({
                    "id": event.get("stId"),
                    "name": primary_name,
                    "equation": equation
                })
                
        print(f"Found {len(reactions)} reactions for pathway {pathway_id}")
        return reactions
        
    except Exception as e:
        print(f"Error fetching Reactome reactions for {pathway_id}: {e}")
        return []

@lru_cache(maxsize=32)
def get_reactome_reaction_details(reaction_id):
    """
    Fetches the specific participants (enzymes, substrates, products) for a single Reactome reaction.
    """
    print(f"Fetching details for Reactome Reaction: {reaction_id}")
    url = f"https://reactome.org/ContentService/data/event/{reaction_id}/participants"
    
    try:
        response = requests.get(url, timeout=15)
        if response.status_code != 200:
            return {"error": f"Failed to fetch details for {reaction_id}"}
            
        data = response.json()
        
        enzymes = set()
        compounds = set()
        
        for entity in data:
            if "refEntities" not in entity:
                continue
                
            for ref in entity["refEntities"]:
                schema = ref.get("schemaClass", "")
                display = ref.get("displayName", "")
                
                # Proteins / Enzymes / Genes
                if schema == "ReferenceGeneProduct":
                    parts = display.split()
                    gene_name = parts[-1] if len(parts) > 1 else display
                    enzymes.add(gene_name)
                    
                # Small Molecules / Compounds
                elif schema == "ReferenceMolecule":
                    clean_name = display.split(" [")[0].strip()
                    compounds.add(clean_name)
                    
        return {
            "id": reaction_id,
            "enzymes": sorted(list(enzymes)),
            "compounds": sorted(list(compounds))
        }
        
    except Exception as e:
        print(f"Error fetching Reactome reaction details for {reaction_id}: {e}")
        return {"error": str(e)}

# --- HOW TO TEST ---
# Run this from the terminal to see it work instantly:
# uv run python -c "import reactome_reactions as rr; print(rr.get_reactome_pathway_reactions('R-HSA-70171'))"
