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

if __name__ == "__main__":
    download_clinvar_database()
