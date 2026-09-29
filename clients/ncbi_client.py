import sys
sys.dont_write_bytecode = True
import requests



def get_ncbi_gene_id(gene_symbol, organism):
    # gene_symbol: The gene symbol or name string to search for in NCBI (e.g. TP53 or tumor protein p53).
    # organism: The organism name string to filter by (e.g. human or Homo sapiens).

    print(f"Searching NCBI Gene database for gene name/symbol: '{gene_symbol}', organism: '{organism}'")

    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    term = f"{gene_symbol}[gene] AND {organism}[orgn]"
    params = {
        "db": "gene",
        "term": term,
        "retmode": "json"
    }

    response = requests.get(base_url, params=params, timeout=10)
    data = response.json()

    id_list = data.get("esearchresult", {}).get("idlist", [])

    if len(id_list) == 0:
        # Try searching the Title field (good for full names like "Sonic Hedgehog")
        title_term = f"{gene_symbol}[Title] AND {organism}[orgn]"
        params["term"] = title_term
        title_res = requests.get(base_url, params=params, timeout=10)
        id_list = title_res.json().get("esearchresult", {}).get("idlist", [])

    if len(id_list) == 0:
        # Last resort: broad free-text search
        fallback_term = f"{gene_symbol} AND {organism}[orgn]"
        params["term"] = fallback_term
        fallback_res = requests.get(base_url, params=params, timeout=10)
        id_list = fallback_res.json().get("esearchresult", {}).get("idlist", [])

    if len(id_list) > 0:
        gene_id = id_list[0]
        print(f"Found NCBI Gene ID: {gene_id}")
        return gene_id

    print(f"No NCBI Gene ID found for symbol: {gene_symbol}")
    return None

def get_ncbi_gene_summary(ncbi_gene_id):
    # ncbi_gene_id: The numeric NCBI Gene ID string (e.g. 7157).

    if ncbi_gene_id is None:
        print("Cannot fetch summary because NCBI Gene ID is None.")
        return None

    print(f"Fetching NCBI gene summary for Gene ID: {ncbi_gene_id}")

    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": "gene",
        "id": ncbi_gene_id,
        "retmode": "json"
    }

    response = requests.get(base_url, params=params, timeout=10)
    data = response.json()

    result = data.get("result", {}).get(str(ncbi_gene_id), {})

    name = result.get("name", "")
    description = result.get("description", "")
    summary = result.get("summary", "")
    organism_official = result.get("organism", {}).get("scientificname", "")

    print(f"Retrieved NCBI summary for {name} ({organism_official})")

    summary_info = {
        "gene_id": ncbi_gene_id,
        "symbol": name,
        "description": description,
        "summary": summary,
        "organism": organism_official
    }

    return summary_info

# Detailed Script Workings:
# 1. get_ncbi_gene_id takes a gene symbol and organism name, construct an NCBI esearch query, and sends an HTTP GET request to NCBI E-utilities.
# 2. It parses the JSON response to extract the first matching NCBI Gene ID.
# 3. get_ncbi_gene_summary takes the numeric NCBI Gene ID, queries NCBI esummary, and extracts the gene symbol, official name, description, summary text, and scientific organism name.
