import sys
sys.dont_write_bytecode = True
import requests
import urllib.parse
import re



def fetch_wikipedia_summary(topic_name):
    # topic_name: The name of the disease or pathway (e.g. "Glycolysis" or "Duchenne muscular dystrophy")
    
    if not topic_name:
        return ""
        
    # Clean the name: split on ';' or '/' and take the first part to avoid URL encoding issues
    clean_name = re.split(r'[;/]', topic_name)[0].strip()
        
    print(f"Fetching Wikipedia summary for: '{clean_name}'")
    
    # Format the name for Wikipedia URL (replace spaces with underscores)
    formatted_name = clean_name.replace(" ", "_")
    encoded_name = urllib.parse.quote(formatted_name)
    
    url = f"https://en.wikipedia.org/api/rest_v1/page/summary/{encoded_name}"
    
    try:
        # Wikipedia requires a User-Agent header, otherwise it blocks default Python requests with a 403 error
        headers = {'User-Agent': 'NCBIAgentBot/1.0'}
        response = requests.get(url, headers=headers, timeout=10)
        
        if response.status_code != 200:
            print(f"No Wikipedia page found for '{topic_name}'.")
            return ""
            
        data = response.json()
        
        # Extract the plain text summary
        extract = data.get("extract", "")
        
        if extract:
            print(f"Successfully retrieved Wikipedia summary for '{topic_name}'.")
            
        return extract
        
    except Exception as e:
        print(f"Error connecting to Wikipedia API: {e}")
        return ""
