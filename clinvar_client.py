import sys
sys.dont_write_bytecode = True
import duckdb

#wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

def fetch_top_mutations(gene_symbols, limit=5):
    """
    Fetches the top pathogenic mutations for a list of genes using DuckDB.
    LIMIT restricts the number of mutations RETURNED PER GENE to keep the LLM context clean.
    """
    if isinstance(gene_symbols, str):
        gene_symbols = [gene_symbols]
        
    print(f"Fetching ClinVar mutations for genes: {gene_symbols} using DuckDB")
    
    import os
    # Dynamically build the path based on the script's location so it works on any computer!
    base_dir = os.path.dirname(os.path.abspath(__file__))
    db_path = os.path.join(base_dir, "Clinvar_files", "variant_summary.txt.gz")
    
    formatted_genes = ", ".join([f"'{g}'" for g in gene_symbols])
    
    # DuckDB can query directly from the gzipped CSV file without loading it into memory!
    # We use all_varchar=true to prevent any strict data type errors while scanning.
    # The ROW_NUMBER() window function ensures we only grab {limit} mutations per gene!
    query = f"""
    WITH UniqueMutations AS (
        SELECT DISTINCT
            GeneSymbol,
            Name as VariantName, 
            "RS# (dbSNP)" as rs_id, 
            Type, 
            PhenotypeList, 
            Origin, 
            ClinicalSignificance
        FROM read_csv_auto('{db_path}', header=True, sep='\\t', ignore_errors=true, all_varchar=true)
        WHERE GeneSymbol IN ({formatted_genes})
          AND ClinicalSignificance ILIKE '%Pathogenic%'
    ),
    RankedMutations AS (
        SELECT *, ROW_NUMBER() OVER(PARTITION BY GeneSymbol ORDER BY VariantName) as rn
        FROM UniqueMutations
    )
    SELECT * EXCLUDE (rn)
    FROM RankedMutations
    WHERE rn <= {limit}
    """
    
    try:
        # Give each concurrent gene search its own private DuckDB connection to prevent thread crashes
        with duckdb.connect() as con:
            result_relation = con.sql(query)
            rows = result_relation.fetchall()
        
        if not rows:
            print("No pathogenic mutations found.")
            return []
            
        mutations_list = []
        for row in rows:
            mut_dict = {
                "gene": str(row[0]),
                "variant_name": str(row[1]),
                "rs_id": str(row[2]),
                "type": str(row[3]),
                "phenotypes": str(row[4]),
                "origin": str(row[5]),
                "significance": str(row[6])
            }
            mutations_list.append(mut_dict)
            
        print(f"Found {len(mutations_list)} pathogenic mutations.")
        #print(mutations_list)
        return mutations_list
        
    except Exception as e:
        print(f"Failed to query ClinVar database with DuckDB: {e}")
        return []


if __name__ == "__main__":
    fetch_top_mutations("BRCA1")

# Script Logic Summary:
# 1. This script uses the incredibly fast DuckDB library to query the massive `variant_summary.txt.gz` file.
# 2. It does NOT load the 443MB file into memory. Instead, it runs a raw SQL query directly against the compressed file on your hard drive!
# 3. It filters for the specific GeneSymbol and for ClinicalSignificance containing 'Pathogenic' using ILIKE for case-insensitivity.
# 4. It returns the top N mutations as a cleanly structured dictionary for the AI agent.
