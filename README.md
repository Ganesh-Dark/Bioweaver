# Bioweaver

Bioweaver is a ReAct-powered biomedical AI agent. It acts as a central hub, dynamically querying and weaving together data from NCBI, KEGG, Reactome, ClinVar and PubTator. It allows users to ask natural language questions about genes, pathways, diseases, compounds, drugs, mutations and retrieve factual answers directly from official biological databases.

> **⚠️ IMPORTANT OVERVIEW DISCLAIMER ⚠️**
> **Bioweaver is designed as an "Overview Agent." Because biological databases are incredibly vast (pathways can have thousands of genes and metabolites), the script is hardcoded to return TRIMMED and PAGINATED results to prevent the AI from crashing or running out of context memory.**
> **If you need comprehensive, exhaustive data (e.g., all 400 genes in a pathway), you MUST manually tweak the script parameters and use your own heavy-duty API keys.**

## Features & Capabilities

Bioweaver is equipped with tools to query specific biological domains:
*   **Pathways & Metabolism (KEGG & Reactome):** Fetches comprehensive lists of genes, compounds, and drugs involved in specific biological pathways. It dynamically merges data from both databases and deduplicates the results.
*   **Diseases & Modules (KEGG & Wikipedia):** Searches the KEGG Disease database to find genes linked to specific diseases, and automatically falls back to Wikipedia to provide human-readable overviews of rare diseases.
*   **Clinical Variants (ClinVar):** Uses a local DuckDB engine to rapidly search the massive ClinVar database for known pathogenic mutations associated with specific genes.
*   **Advanced Literature NLP (PubTator3):** Integrates with the official NCBI PubTator3 API to scan the latest PubMed literature. It uses advanced NLP to automatically extract mentioned chemicals, genes, diseases, cell lines, species, variants, and biological relationships.
*   **Gene Deep-Dives (NCBI HGNC):** Retrieves official gene summaries, descriptions, and organism data directly from NCBI.
*   **Reactions & Enzymes (KEGG Offline TSVs):** Local parsing of KEGG Reactions, Enzymes, Modules, and Orthology (KO) using offline `.tsv` databases to isolate specific substrates and products without API rate limits.
*   **Local Pathway Analytics:** A custom Python module (`kegg_analytics.py`) that performs set operations and pathway analysis locally (e.g., finding isozymes, shared compounds, orphan reactions). This reduces token usage by preprocessing the datasets before passing context to the LLM.
*   **Query Routing:** The ReAct agent differentiates between database queries and general biological concepts. Standard conceptual questions (e.g., "What are the 10 steps of glycolysis?") are answered using the LLM's internal knowledge, saving API calls.

## Upcomings!!

*   **Reactome-based pathway analysis** is currently in development and will be released in upcoming commits!

## Setup Instructions

**Important Note on Data Files:**
Almost all required database files (HGNC, KEGG Pathways, KEGG Diseases, KEGG Organisms, and Reactome) are fully included in this repository in their respective folders. 

The **only** file that is not included is the ClinVar database (`variant_summary.txt.gz`). This file is over 150MB, and GitHub strictly blocks the uploading of any files larger than 100MB.

Before running the application for the first time, you must download the ClinVar database by running this simple script in your terminal:
```bash
python setup_dbs.py
```
*This will safely download the file directly from the NCBI FTP servers and place it in the `Clinvar_files/` directory.*

**API Key Setup:**
This application is powered by a Google Gemini LLM. You must provide your API key for the AI to function.
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

### Literature & Advanced NLP (PubTator) Queries
17. *"Look up literature on TP53 and Glioblastoma. Extract any FDA-approved drugs or chemicals mentioned."*
18. *"Search the literature for 'BRAF Melanoma'. Extract the specific Cell Lines they used for their experiments and the Species studied."*
19. *"Find recent papers on 'Carbamazepine Epilepsy'. Look for any explicitly extracted NLP Relationships between the chemicals and the diseases."*
20. *"Look up 'CRISPR Cas9 Duchenne Muscular Dystrophy'. Tell me exactly which Species were mentioned in the papers and what specific variants they targeted."*

### Reaction & Enzyme Queries
21. *"Find all reactions where pyruvate is a product."*
22. *"What are the substrates and products of reaction R00224?"*
23. *"List enzymes that involve ATP as a substrate."*
24. *"What is a KEGG module and how is it related to glycolysis?"*

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
*   **NCBI / NLM:** Gene summaries, PubTator3 literature mining, and ClinVar variant data are sourced from the National Center for Biotechnology Information.
*   **Reactome:** Additional pathway data is integrated from the open-source Reactome knowledgebase. Note: Bioweaver uses a custom local Machine Learning model (TF-IDF vectorizer) to perform intelligent, fuzzy semantic matching of Reactome pathway names, preventing exact-character mismatch errors.
