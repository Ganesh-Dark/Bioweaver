import os
import urllib.request


def download_clinvar_database():
    """
    Downloads the official ClinVar variant summary file directly from the NCBI FTP server.
    """
    
    clinvar_url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
    
    base_dir = os.path.dirname(os.path.abspath(__file__))
    clinvar_dir = os.path.join(base_dir, "Clinvar_files")

    if not os.path.exists(clinvar_dir):
        print(f"Creating directory: {clinvar_dir}")
        os.makedirs(clinvar_dir)
        
    destination_file = os.path.join(clinvar_dir, "variant_summary.txt.gz")
    
    if os.path.exists(destination_file):
        print(f"The file {destination_file} already exists. Skipping download.")
        return

    print(f"Starting download of ClinVar database from NCBI...")
    print(f"URL: {clinvar_url}")
    print(f"Please be patient, this is a large file (approx 150MB+)...")
    
    try:
        urllib.request.urlretrieve(clinvar_url, destination_file)
        print("\nDownload complete!")
        print(f"Saved securely to: {destination_file}")
    except Exception as e:
        print(f"\nAn error occurred during download: {e}")

def download_kegg_databases():
    """
    Downloads KEGG databases using the KEGG REST API and formats them.
    """
    base_dir = os.path.dirname(os.path.abspath(__file__))
    kegg_dir = os.path.join(base_dir, "Kegg_files")
    
    if not os.path.exists(kegg_dir):
        print(f"Creating directory: {kegg_dir}")
        os.makedirs(kegg_dir)
        
    kegg_files = {
        "enzyme.tsv": "enzyme",
        "kegg_compounds_synonyms.tsv": "compound",
        "kegg_diseases.tsv": "disease",
        "kegg_organisms.tsv": "organism",
        "ko.tsv": "ko",
        "kegg_pathways.tsv": "pathway",
        "module.tsv": "module",
        "Reactions.tsv": "rn",
        # Cross-reference files used for KEGG<->Reactome reaction matching
        "kegg_cpd_to_chebi.tsv": "conv/chebi/cpd",
        "kegg_rn_to_cpd.tsv": "link/cpd/rn"
    }
    
    for filename, endpoint in kegg_files.items():
        destination_file = os.path.join(kegg_dir, filename)
        if os.path.exists(destination_file):
            print(f"The file {destination_file} already exists. Skipping download.")
            continue
            
        # kegg_cpd_to_chebi and kegg_rn_to_cpd use different base URL patterns
        if endpoint.startswith("conv/") or endpoint.startswith("link/"):
            url = f"https://rest.kegg.jp/{endpoint}"
        else:
            url = f"https://rest.kegg.jp/list/{endpoint}"
        print(f"Downloading KEGG {filename} from {url}...")
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=60) as response:
                content = response.read().decode('utf-8')
            
            processed_lines = []
            for line in content.split('\n'):
                if not line.strip():
                    continue
                cols = line.split('\t')
                #strip KEGG prefixes from the first column (e.g., 'path:map01100' -> 'map01100')
                if len(cols) > 0 and ':' in cols[0]:
                    cols[0] = cols[0].split(':', 1)[1]
                processed_lines.append('\t'.join(cols))
            
            with open(destination_file, 'w', encoding='utf-8') as f:
                f.write('\n'.join(processed_lines) + '\n')
            
            print(f"Saved: {destination_file}")
        except Exception as e:
            print(f"Error downloading {filename}: {e}")

def download_reactome_chebi():
    """
    Downloads the official Reactome ChEBI-to-Reactions cross-reference file.
    Used to match KEGG reactions with Reactome reactions via shared ChEBI compound IDs.
    """
    # URL: Reactome's official download server
    url = "https://reactome.org/download/current/ChEBI2Reactome_PE_Reactions.txt"
    base_dir = os.path.dirname(os.path.abspath(__file__))
    reactome_dir = os.path.join(base_dir, "Reactome_files")

    if not os.path.exists(reactome_dir):
        print(f"Creating directory: {reactome_dir}")
        os.makedirs(reactome_dir)

    destination_file = os.path.join(reactome_dir, "ChEBI2Reactome_PE_Reactions.txt")

    if os.path.exists(destination_file):
        print(f"The file {destination_file} already exists. Skipping download.")
        return

    print(f"Downloading Reactome ChEBI cross-reference file (~41MB)...")
    print(f"URL: {url}")
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=120) as response:
            content = response.read()
        with open(destination_file, 'wb') as f:
            f.write(content)
        print(f"Saved: {destination_file}")
    except Exception as e:
        print(f"Error downloading ChEBI2Reactome file: {e}")


def download_hgnc_database():
    """
    Downloads the complete HGNC dataset for gene symbol resolution.
    """
    hgnc_url = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
    base_dir = os.path.dirname(os.path.abspath(__file__))
    hgnc_dir = os.path.join(base_dir, "hgnc_files")
    
    if not os.path.exists(hgnc_dir):
        print(f"Creating directory: {hgnc_dir}")
        os.makedirs(hgnc_dir)
        
    destination_file = os.path.join(hgnc_dir, "hgnc_complete_set.txt")
    
    if os.path.exists(destination_file):
        print(f"The file {destination_file} already exists. Skipping download.")
        return

    print(f"Starting download of HGNC database...")
    print(f"URL: {hgnc_url}")
    try:
        urllib.request.urlretrieve(hgnc_url, destination_file)
        print("Download completed")
        print(f"Saved securely to: {destination_file}")
    except Exception as e:
        print(f"An error occurred during HGNC download: {e}")


def download_rhea_files():
    """
    Downloads the Rhea reaction cross-reference files.
    These are used to map KEGG reaction IDs to Reactome reaction IDs via a shared Rhea master ID.
    """
    # Base URL for Rhea FTP server
    rhea_files = {
        "rhea2kegg_reaction.tsv": "https://ftp.expasy.org/databases/rhea/tsv/rhea2kegg_reaction.tsv",
        "rhea2reactome.tsv": "https://ftp.expasy.org/databases/rhea/tsv/rhea2reactome.tsv"
    }
    base_dir = os.path.dirname(os.path.abspath(__file__))
    rhea_dir = os.path.join(base_dir, "Rhea_files")

    if not os.path.exists(rhea_dir):
        print(f"Creating directory: {rhea_dir}")
        os.makedirs(rhea_dir)

    for filename, url in rhea_files.items():
        destination_file = os.path.join(rhea_dir, filename)
        if os.path.exists(destination_file):
            print(f"The file {destination_file} already exists. Skipping download.")
            continue
        print(f"Downloading Rhea file: {filename} from {url}...")
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=60) as response:
            content = response.read()
        with open(destination_file, 'wb') as f:
            f.write(content)
        print(f"Saved: {destination_file}")


if __name__ == "__main__":
    download_clinvar_database()
    download_kegg_databases()
    download_hgnc_database()
    download_reactome_chebi()
    download_rhea_files()
