# Bioweaver

Bioweaver is a ReAct-powered biomedical AI agent. It acts as a central hub, dynamically querying and weaving together data from NCBI, KEGG, Reactome, and ClinVar. It allows users to ask natural language questions about genes, pathways, diseases, compounds, drugs, mutations and retrieve factual answers directly from official biological databases.

> **⚠️ IMPORTANT OVERVIEW DISCLAIMER ⚠️**
> **Bioweaver is designed as an "Overview Agent." Because biological databases are incredibly vast (pathways can have thousands of genes and metabolites), the script is hardcoded to return TRIMMED and PAGINATED results to prevent the AI from crashing or running out of context memory.**
> **If you need comprehensive, exhaustive data (e.g., all 400 genes in a pathway), you MUST manually tweak the script parameters and use your own heavy-duty API keys.**

## Features & Capabilities

Bioweaver is equipped with tools to query specific biological domains:
*   **Pathways & Metabolism (KEGG & Reactome):** Fetches comprehensive lists of genes, compounds, and drugs involved in specific biological pathways. It dynamically merges data from both databases and deduplicates the results.
*   **Diseases & Modules (KEGG & Wikipedia):** Searches the KEGG Disease database to find genes linked to specific diseases, and automatically falls back to Wikipedia to provide human-readable overviews of rare diseases.
*   **Clinical Variants (ClinVar):** Uses a local DuckDB engine to rapidly search the massive ClinVar database for known pathogenic mutations associated with specific genes.
*   **Gene Deep-Dives (NCBI HGNC):** Retrieves official gene summaries, descriptions, and organism data directly from NCBI.
*   **Reactions & Enzymes (KEGG Offline TSVs):** Local parsing of KEGG Reactions, Enzymes, Modules, and Orthology (KO) using offline `.tsv` databases to isolate specific substrates and products without API rate limits.
*   **Local Pathway Analytics:** Custom Python modules (`kegg_analytics.py` & `reactome_analytics.py`) that perform set operations and pathway analysis locally (e.g., finding isozymes, shared compounds, orphan reactions). This reduces token usage by preprocessing datasets.
*   **Chemical Resolution:** A robust internal resolver using PubChem to map generic chemical names to universal InChIKeys and extract official KEGG, ChEBI, HMDB, and DrugBank identifiers.
*   **Rhea Universal Reaction Mapping:** Seamless, offline cross-referencing between KEGG, Reactome, and Rhea. Instantly maps reactions to universal Rhea Master IDs, ChEBI participants, UniProt proteins, and EC numbers.
*   **Query Routing:** The ReAct agent differentiates between database queries and general biological concepts. Standard conceptual questions (e.g., "What are the 10 steps of glycolysis?") are answered using the LLM's internal knowledge, saving API calls.

## Setup Instructions

**Important Note on Data Files:**
Almost all required database files (HGNC, KEGG Pathways, KEGG Diseases, KEGG Organisms, and Reactome) are fully included in this repository in their respective folders. 

The **only** file that is not included is the ClinVar database (`variant_summary.txt.gz`). This file is over 150MB, and GitHub strictly blocks the uploading of any files larger than 100MB.

Before running the application for the first time, you must download the ClinVar database by running this simple script in your terminal:
```bash
python setup_dbs.py
```
*This will safely download the ClinVar file directly from the NCBI FTP servers and place it in your local directory.*

> **Want the latest KEGG & HGNC data?** 
> Even though KEGG and HGNC databases are already included, they update frequently. If you ever want to fetch the absolute latest, most recent versions directly from the official KEGG/HGNC APIs, simply delete the existing files in your data folder and run `python setup_dbs.py` again. The script will automatically detect they are missing and download the newest versions for you!

**API Key Setup:**
I strongly recommend to use Google free tier gemini API as this tool sometimes consumes too many tokens depending on context. You must provide your API key for that. There is also option to use Hugging face models or Groq if you prefer.
1. Create a file named `.env` in the root directory.
2. Add your API key exactly like this:
```env
LLM_API=your_actual_api_key_here
```

## Installation and Execution

This project uses the modern `uv` package manager for lightning-fast dependency resolution.

To install all dependencies and launch the application instantly, you can just run:
```bash
uv run python app.py
```
*(This will start the backend server on **http://localhost:8000**)*

Alternatively, to explicitly create a virtual environment and install from `requirements.txt` using `uv`:
```bash
uv venv
source .venv/bin/activate
uv pip install -r requirements.txt
python app.py
```
*(This will start the backend server on **http://localhost:8000**)*

### Alternative Installation (using pip)
If you do not have `uv` installed, you can use standard Python tools:
```bash
python -m venv venv
source venv/bin/activate  # On Windows use: venv\Scripts\activate
pip install -r requirements.txt
python app.py
```

## Example Queries (Test the AI!)

Copy and paste these queries into the Bioweaver chat interface to see its full power:

### Pathway & Metabolism Queries
1. *"List all the genes, compounds, and drugs involved in carbohydrate metabolism."*
2. *"I want to know the genes in the TCA cycle. I only want data from Reactome."*
3. *"What compounds are created during glycolysis?"*
4. *"Show me the difference between KEGG and Reactome data for the pentose phosphate pathway."*

### Disease & Clinical Queries
5. *"Search for Parkinson's disease. Get me the overview and find all linked genes from KEGG."*
6. *"What genes are associated with Alzheimer's disease in the KEGG database?"*
7. *"Tell me about Glioblastoma and list any related metabolic pathways."*
8. *"What does the KEGG database say about Type 2 Diabetes?"*

### Gene Deep-Dives
9. *"Tell me everything you know about the BRCA1 gene. What is its official NCBI summary?"*
10. *"Look up the TP53 gene. What diseases is it linked to in KEGG?"*
11. *"What is the official description and organism for the EGFR gene?"*
12. *"Compare the official NCBI summaries of BRCA1 and BRCA2."*

### ClinVar Mutation Queries
13. *"What are the most dangerous pathogenic mutations for the TP53 gene?"*
14. *"Find the top pathogenic variants for BRCA1 in ClinVar."*
15. *"Are there any known mutations for the APOE gene?"*
16. *"List the pathogenic mutations for EGFR and tell me what phenotypes they cause."*

### Reaction, Chemical & Cross-Database Queries
17. *"Find all reactions where pyruvate is a product."*
18. *"What are the substrates and products of reaction R00224?"*
19. *"List enzymes that involve ATP as a substrate."*
20. *"What is a KEGG module and how is it related to glycolysis?"*
21. *"Resolve 'pyruvate' and give me all its database IDs."*
22. *"Get all reactions from glycolysis and get their rhea ids."*
23. *"What exactly is Rhea 10012?"*
24. *"What are the common reactions between KEGG Glycolysis and Reactome Glycolysis?"*

### Advanced Pathway Analytics
25. *"Which ECs catalyze the same reaction in the TCA cycle?"*
26. *"What are the irreversible reactions in the Pentose Phosphate Pathway?"*
27. *"Which reactions in the Alzheimer's disease pathway do not have a known enzyme yet?"*
28. *"What compounds are common between Glycolysis, the TCA cycle, and the Pentose Phosphate Pathway?"*
29. *"Exactly which reactions produce ATP, and which reactions consume ATP in Glycolysis?"*
30. *"Break down the enzymes in Glycolysis by their EC class. How many are oxidoreductases?"*

### General & Conceptual Biology
*Bioweaver routes classic textbook questions directly to the LLM's internal knowledge rather than querying databases.*
31. *"What are the basic steps of glycolysis in order?"*
32. *"Explain the central dogma of molecular biology."*
33. *"What is the difference between a kinase and a phosphatase?"*
34. *"How do competitive and non-competitive inhibitors differ?"*

## Data Acknowledgments & Credits

This project queries and processes data from several public biological databases. Please refer to their respective licenses and terms of use if you plan to use Bioweaver for commercial purposes:

*   **KEGG (Kyoto Encyclopedia of Genes and Genomes):** Pathway, disease, reaction, and enzyme data are sourced from KEGG. KEGG data is freely available for academic use, but commercial use requires a license. Please visit the [KEGG Legal & Copyright](https://www.kegg.jp/kegg/legal.html) page for full details.
*   **NCBI / NLM:** Gene summaries and ClinVar variant data are sourced from the National Center for Biotechnology Information.
*   **Reactome:** Additional pathway data is integrated from the open-source Reactome knowledgebase. Note: Bioweaver uses a custom local Machine Learning model (TF-IDF vectorizer) to perform intelligent, fuzzy semantic matching of Reactome pathway names, preventing exact-character mismatch errors.
