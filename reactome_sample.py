import sys
sys.dont_write_bytecode = True
import requests

def fetch_reactome_participants(pathway_id):
    """
    Fetches all physical entities (genes, proteins, compounds) involved in a Reactome pathway.
    """
    print(f"Fetching Reactome participants for pathway: {pathway_id}")
    url = f"https://reactome.org/ContentService/data/participants/{pathway_id}"
    
    try:
        response = requests.get(url, timeout=15)
        response.raise_for_status()
        data = response.json()
        
        genes = set()
        compounds = set()
        
        # Reactome returns a list of PhysicalEntities (complexes, sets, single molecules).
        # We need to iterate through them and pull out the specific "refEntities" (the actual molecules/genes).
        for entity in data:
            if "refEntities" in entity:
                for ref in entity["refEntities"]:
                    # Is it a Gene/Protein?
                    if ref.get("schemaClass") == "ReferenceGeneProduct":
                        name = ref.get("displayName", "Unknown Gene")
                        # e.g., "UniProt:P30153 PPP2R1A" -> Extract just the gene name "PPP2R1A" if possible
                        parts = name.split()
                        if len(parts) > 1:
                            genes.add(parts[-1])
                        else:
                            genes.add(name)
                            
                    # Is it a Chemical Compound?
                    elif ref.get("schemaClass") == "ReferenceMolecule":
                        name = ref.get("displayName", "Unknown Compound")
                        # Reactome already formats it perfectly as "Compound Name [ChEBI:ID]"
                        compounds.add(name)
                            
        # Sort for clean output
        genes = sorted(list(genes))
        compounds = sorted(list(compounds))
        
        print("\n" + "="*50)
        print(f"REACTOME PATHWAY: {pathway_id}")
        print("="*50)
        
        print(f"\nCORE GENES/PROTEINS ({len(genes)} found):")
        for g in genes:
            print(f"  - {g}")
            
        print(f"\ncORE COMPOUNDS ({len(compounds)} found):")
        for c in compounds:
            print(f"  - {c}")
                   
    except Exception as e:
        print(f"Error fetching Reactome data: {e}")

if __name__ == "__main__":
    # Test it on the strict Human Glycolysis pathway (R-HSA-70171)
    fetch_reactome_participants("R-HSA-70171")
