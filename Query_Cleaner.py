import sys
sys.dont_write_bytecode = True
import os
import json
from langchain_core.prompts import ChatPromptTemplate
from langchain_google_genai import ChatGoogleGenerativeAI
from dotenv import load_dotenv, find_dotenv

load_dotenv()



def query_cleaner(input_from_user):
    # input_from_user: The raw user question string to extract biomedical entities from.

    prompt = ChatPromptTemplate.from_template("""
    You are a biomedical entity extraction and normalization agent for NCBI and KEGG.

    Extract entities from the user's query, correct any spelling mistakes, and normalize informal phrases into standard medical search keywords.

    Categories:
    - Gene: official gene symbols, names, or synonyms (e.g., TP53, BRCA1, tumor protein p53)
    - Pathway: biological or metabolic pathway names
    - Disease: standard medical disease names (e.g., convert "ovary related diseases" to "ovarian cancer" or "ovary"; convert "lung related cancer" to "lung cancer")
    - Symptom: clinical symptoms, phenotypes, physical signs, or clinical findings (e.g., seizures, microcephaly, developmental delay, hypotonia, dyspraxia)
    - Organism: species or organism names (e.g., human, mouse)

    Rules:
    - Automatically fix spelling mistakes and typos.
    - Convert descriptive informal phrases (like "ovary related diseases") into formal medical search keywords (like "ovarian cancer", "ovary").
    - If the user asks for a broad category (e.g., "carbohydrate related pathways"), extract specific, well-known pathway names belonging to that category (e.g., "glycolysis", "gluconeogenesis", "pentose phosphate pathway", "TCA cycle") instead of just the broad category name.
    - Analyze the user's question to intelligently determine ONLY the necessary databases to fulfill their request. Options are: NCBI (for gene details), KEGG (for pathways, diseases, compounds, drugs), ClinVar (for mutations/variants), PubMed (for research papers). If the user asks for a specific database (e.g., "from ncbi"), list only that one. If they do not specify, infer the minimum required databases based on the entities they are asking about (e.g. asking for "compounds of a pathway" requires ONLY KEGG. Asking for "mutations in TP53" requires ClinVar and NCBI). Do not list all databases unless absolutely necessary.
    - Do not answer the question.
    - Do not explain or infer anything.
    - Remove duplicates.
    - Return an empty array when no entity is found.
    - Return ONLY valid JSON.
    - Use exactly these keys.

    Output format:
    {{
    "Target_Databases": [],
    "Gene": [],
    "Pathway": [],
    "Disease": [],
    "Symptom": [],
    "Organism": []
    }}

    Question:
    {question}
    """)

    final_prompt = prompt.invoke({
        "question": input_from_user
    })

    api = os.getenv("LLM_API")

    llm = ChatGoogleGenerativeAI(
        model="gemini-flash-lite-latest",
        temperature=0, api_key = api
    )

    response = llm.invoke(final_prompt)

    raw_text = response.content[0]['text']
    print(f"Extracted Entities Raw: {raw_text}")

    clean_text = raw_text.strip()
    if clean_text.startswith("```"):
        lines = clean_text.split("\n")
        if lines[0].startswith("```"):
            lines = lines[1:]
        if len(lines) > 0 and lines[-1].strip().startswith("```"):
            lines = lines[:-1]
        clean_text = "\n".join(lines).strip()

    extracted_data = json.loads(clean_text)
    return extracted_data

# Detailed Script Workings:
# 1. query_cleaner takes a user input question and formats a ChatPromptTemplate.
# 2. It sends the prompt to ChatGoogleGenerativeAI (Gemini) to extract biomedical entity categories (Gene, Pathway, Disease, Organism).
# 3. It parses the JSON output string from Gemini into a Python dictionary and returns it.
    