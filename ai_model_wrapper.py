import os
from dotenv import load_dotenv
from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_groq import ChatGroq
from langchain_huggingface import ChatHuggingFace, HuggingFaceEndpoint

# Load environment variables
load_dotenv()

def get_llm():
    """
    Reads .env configuration and initializes the appropriate LLM.
    Returns the instantiated LLM object.
    """
    api_key = os.getenv("LLM_API")
    groq_api = os.getenv("GROQ_API_KEY")
    hf_token = os.getenv("HF_TOKEN")
    
    use_hf = os.getenv("USE_HF", "False").lower() in ("true", "1", "yes")
    use_groq = os.getenv("USE_GROQ", "False").lower() in ("true", "1", "yes")

    if not api_key and not use_groq and not use_hf:
        print("Warning: No LLM configuration found in .env. Attempting fallback to Gemini.")

    if use_hf:
        if not hf_token:
            raise ValueError("USE_HF is True, but HF_TOKEN is missing from .env")
        llm_endpoint = HuggingFaceEndpoint(
            repo_id="openai/gpt-oss-20b",
            huggingfacehub_api_token=hf_token,
            task="text-generation",
            temperature=0.2,
            max_new_tokens=4096
        )
        llm = ChatHuggingFace(llm=llm_endpoint)
        print("Backend: Loaded HUGGING FACE (openai/gpt-oss-20b)")
    elif use_groq:
        if not groq_api:
            raise ValueError("USE_GROQ is True, but GROQ_API_KEY is missing from .env")
        llm = ChatGroq(
            model="openai/gpt-oss-120b",
            temperature=0.2,
            api_key=groq_api,
            max_retries=0, max_tokens=4096
        )
        print("Backend: Loaded GROQ (openai/gpt-oss-120b)")
    else:
        if not api_key:
            raise ValueError("LLM_API (Gemini key) is missing from .env")
        llm = ChatGoogleGenerativeAI(
            model="gemini-flash-lite-latest",
            temperature=0.2,
            api_key=api_key,
            max_retries=0
        )
        print("Backend: Loaded GEMINI (gemini-3.5-flash)")
        
    return llm
