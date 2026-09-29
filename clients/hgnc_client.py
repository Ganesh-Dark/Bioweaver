import sys
sys.dont_write_bytecode = True
import os
#wget https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt

# Global cache to hold the loaded HGNC data
HGNC_CACHE = {}
HGNC_FILE = os.path.join(os.path.dirname(__file__), "hgnc_files", "hgnc_complete_set.txt")

def _load_hgnc_data():
    """
    Loads the HGNC database from the local text file into memory.
    Only runs once. Creates a mapping of {symbol: full_name}.
    """
    if HGNC_CACHE:
        return
        
    print(f"Loading HGNC database from {HGNC_FILE}...")
    try:
        with open(HGNC_FILE, "r", encoding="utf-8") as f:
            headers = f.readline().strip().split("\t")
            # Find index of symbol and name columns
            try:
                symbol_idx = headers.index("symbol")
                name_idx = headers.index("name")
            except ValueError:
                print("Error: Could not find 'symbol' or 'name' columns in HGNC file.")
                return
                
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) > max(symbol_idx, name_idx):
                    symbol = parts[symbol_idx].strip()
                    name = parts[name_idx].strip()
                    if symbol:
                        HGNC_CACHE[symbol] = name
        print(f"Successfully loaded {len(HGNC_CACHE)} genes from HGNC.")
    except Exception as e:
        print(f"Error loading HGNC file: {e}")

def expand_gene_symbol(symbol):
    """
    Takes a short gene symbol (e.g. 'PFKFB4') and returns the formatted
    expansion 'PFKFB4 (6-phosphofructo-2-kinase...)' if found in HGNC.
    """
    # symbol: The short gene symbol to expand.
    if not HGNC_CACHE:
        _load_hgnc_data()
        
    full_name = HGNC_CACHE.get(symbol)
    if full_name:
        return f"{symbol} ({full_name})"
    return symbol
