import os
import json
import logging
from dotenv import load_dotenv
from langchain_core.prompts import PromptTemplate
from langchain_core.output_parsers import JsonOutputParser
from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_groq import ChatGroq
from langchain_huggingface import ChatHuggingFace, HuggingFaceEndpoint
from pydantic import BaseModel, Field

load_dotenv()

class RouterConfig(BaseModel):
    targets: list = Field(description="A list of target gene symbols (e.g., ['TP53', 'BRCA1']) or RS IDs (e.g., ['121918399']). If only one target, put it in a list.")
    target_type: str = Field(description="'gene' if the targets are gene symbols, or 'rs_id' if they are RS IDs.")
    mutation_type: str = Field(description="The specific type of mutation requested. If target_type is 'gene', default to 'single nucleotide variant'. If target_type is 'rs_id', default to an empty string ''.")
    fallback_allowed: bool = Field(description="True if we can fallback to other mutation types if the requested one is not found.")

def route_clinvar_query(user_query, default_targets=None):
    """
    Sends the user's raw query to the Gemini Router.
    Returns a parsed JSON configuration with routing instructions for DuckDB.
    
    Parameters:
    user_query: The natural language string asked by the user.
    default_targets: A list of backup gene symbols to use if the AI cannot find any in the query.
    """
    if default_targets is None:
        default_targets = []
        
    fallback_config = {
        "targets": default_targets,
        "target_type": "gene",
        "mutation_type": "single nucleotide variant",
        "fallback_allowed": True
    }
    
    api_key = os.getenv("LLM_API")
    groq_api = os.getenv("GROQ_API_KEY")
    hf_token = os.getenv("HF_TOKEN")
    use_hf = os.getenv("USE_HF", "False").lower() in ("true", "1", "yes")
    use_groq = os.getenv("USE_GROQ", "False").lower() in ("true", "1", "yes")

    if not api_key and not use_groq and not use_hf:
        print("Warning: LLM_API not found in .env. Falling back to default SNV search.")
        return fallback_config

    try:
        from ai_model_wrapper import get_llm
        llm = get_llm()
        
        parser = JsonOutputParser(pydantic_object=RouterConfig)
        
        prompt = PromptTemplate(
            template="You are a strict JSON API router for a biomedical database called ClinVar.\n"
                     "Analyze the user's query and extract the target gene or RS ID, and the exact mutation type they want.\n"
                     "{format_instructions}\n\nUser Query: {query}\n",
            input_variables=["query"],
            partial_variables={"format_instructions": parser.get_format_instructions()},
        )
        
        chain = prompt | llm | parser
        config = chain.invoke({"query": user_query})
        
        # Sanitize targets for rs_id
        if config.get("target_type") == "rs_id":
            sanitized = []
            for t in config.get("targets", []):
                if isinstance(t, str) and t.lower().startswith("rs"):
                    sanitized.append(t[2:])
                else:
                    sanitized.append(t)
            config["targets"] = sanitized
            
        return config
        
    except Exception as e:
        print(f"Router failed: {e}. Falling back to default SNV search.")
        return fallback_config

if __name__ == "__main__":
    print("Test 1 (Gene):", route_clinvar_query("Are there any copy number loss mutations for TP53?", ["TP53"]))
    print("Test 2 (RS ID):", route_clinvar_query("Search ClinVar for rs121918399", []))
