import sys
sys.dont_write_bytecode = True
import duckdb

#wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

def fetch_top_mutations(targets: list, target_type="gene", mutation_type="single nucleotide variant", fallback_allowed=True, limit=5):
    """
    Fetches the top N most clinically significant pathogenic mutations for given targets.
    Supports smart routing: searching by GeneSymbol OR RS ID, and filtering by mutation type.
    
    Parameters:
    targets: A list of gene names or RS IDs to search for simultaneously.
    target_type: Specifies if the targets are "gene" or "rs_id".
    mutation_type: The exact type of mutation to filter by (e.g., "single nucleotide variant").
    fallback_allowed: If True, it retries the search without the mutation_type filter if 0 results are found.
    limit: The maximum number of results to return PER target.
    """
    
    print(f"Fetching ClinVar mutations for targets: {targets} (Type: {target_type}) using DuckDB")
    if mutation_type:
        print(f"Applying strict mutation filter: '{mutation_type}'")
    
    import os
    base_dir = os.path.dirname(os.path.abspath(__file__))
    db_path = os.path.join(base_dir, "Clinvar_files", "variant_summary.txt.gz")
    
    # Format the targets for SQL IN clause (e.g., 'LRRK2', 'SNCA')
    formatted_targets = ", ".join([f"'{t}'" for t in targets])
    
    # Determine the WHERE clause dynamically based on what the HF Router extracted
    if target_type == "rs_id":
        target_clause = f"\"RS# (dbSNP)\" IN ({formatted_targets})"
    else:
        # Use a list of ILIKE statements or an IN clause (IN is exact, but GeneSymbol is usually exact)
        # However, to maintain case insensitivity across a list in DuckDB, we can just upper() it
        formatted_upper_targets = ", ".join([f"'{str(t).upper()}'" for t in targets])
        target_clause = f"UPPER(GeneSymbol) IN ({formatted_upper_targets})"
        
    type_clause = f"AND Type ILIKE '%{mutation_type}%'" if mutation_type else ""
    
    # Helper function to run the query
    def run_query(current_type_clause):
        return f"""
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
            WHERE {target_clause}
              AND ClinicalSignificance ILIKE '%Pathogenic%'
              {current_type_clause}
        ),
        RankedMutations AS (
            SELECT *, ROW_NUMBER() OVER(PARTITION BY GeneSymbol ORDER BY VariantName) as rn
            FROM UniqueMutations
        )
        SELECT * EXCLUDE (rn)
        FROM RankedMutations
        WHERE rn <= {limit}
        """
        
    query = run_query(type_clause)
    
    try:
        # Give each concurrent gene search its own private DuckDB connection to prevent thread crashes
        with duckdb.connect() as con:
            result_relation = con.sql(query)
            rows = result_relation.fetchall()
        
        if not rows and fallback_allowed and mutation_type:
            print(f"0 rows found for '{mutation_type}'. Smart Fallback Activated: Retrying without type filter...")
            fallback_query = run_query("")
            with duckdb.connect() as con:
                rows = con.sql(fallback_query).fetchall()

        if not rows:
            print("No pathogenic mutations found.")
            return []
            
        mutations_list = []
        for row in rows:
            raw_rs_id = str(row[2])
            clean_rs_id = "N/A" if raw_rs_id in ["-1", "-", "na"] else raw_rs_id
            
            mut_dict = {
                "gene": str(row[0]),
                "variant_name": str(row[1]),
                "rs_id": clean_rs_id,
                "type": str(row[3]),
                "phenotypes": str(row[4]),
                "origin": str(row[5]),
                "significance": str(row[6])
            }
            mutations_list.append(mut_dict)
            
        print(f"Found {len(mutations_list)} pathogenic mutations.")
        print(mutations_list)
        return mutations_list
        
    except Exception as e:
        print(f"Failed to query ClinVar database with DuckDB: {e}")
        return []


if __name__ == "__main__":
    # Test fallback logic
    fetch_top_mutations("TP53", target_type="gene", mutation_type="copy number loss", fallback_allowed=True)
    # Test RS ID routing
    fetch_top_mutations("563571689", target_type="rs_id", mutation_type="")

# Script Logic Summary:
# 1. This script uses the incredibly fast DuckDB library to query the massive `variant_summary.txt.gz` file.
# 2. It does NOT load the 443MB file into memory. Instead, it runs a raw SQL query directly against the compressed file on your hard drive!
# 3. It filters for the specific GeneSymbol and for ClinicalSignificance containing 'Pathogenic' using ILIKE for case-insensitivity.
# 4. It returns the top N mutations as a cleanly structured dictionary for the AI agent.
