import sys
sys.dont_write_bytecode = True
import streamlit as st
from agent_graph import build_ncbi_kegg_graph
import sys
import random

# Must be the first Streamlit command
st.set_page_config(page_title="Bioweaver", page_icon="", layout="wide")

st.title("Bioweaver")
st.markdown("Ask anything about genes, pathways, diseases, and mutations. The AI will dynamically query NCBI, KEGG, ClinVar, and PubMed to answer.")

import os
# Check if databases are installed
base_dir = os.path.dirname(os.path.abspath(__file__))
clinvar_path = os.path.join(base_dir, "Clinvar_files", "variant_summary.txt.gz")

missing_dbs = []
if not os.path.exists(clinvar_path):
    missing_dbs.append("ClinVar")

if missing_dbs:
    st.error(f"**Missing Databases Detected:** {', '.join(missing_dbs)}")
    st.warning("Please run `python setup_dbs.py` in your terminal to download the required datasets before using the agent!")
    st.stop() # Stops execution so the app doesn't crash later

# Initialize chat history
if "messages" not in st.session_state:
    st.session_state.messages = []

# Initialize the ReAct agent
if "agent" not in st.session_state:
    # We load the agent once and store it in session state for speed
    st.session_state.agent = build_ncbi_kegg_graph()

# Display chat messages from history on app rerun
for message in st.session_state.messages:
    with st.chat_message(message["role"]):
        st.markdown(message["content"])

recommendations = [
    "What are the mutations for BRCA1?",
    "Genes involved in glycolysis",
    "What drugs treat Alzheimer's disease?",
    "Show me the TCA cycle pathway",
    "Are there any mutations in the TP53 gene?"
]

ph = f"E.g., '{random.choice(recommendations)}'"

# React to user input. We MUST provide a static 'key' here so Streamlit doesn't think 
# we are creating a brand new input box every time the random placeholder changes!
if prompt := st.chat_input(ph, key="main_chat_input"):
    
    # Display user message in chat message container
    st.chat_message("user").markdown(prompt)
    
    # Add user message to chat history
    st.session_state.messages.append({"role": "user", "content": prompt})

    # Display assistant response in chat message container
    with st.chat_message("assistant"):
        with st.spinner("Analyzing ..."):
            try:
                # We will use an empty placeholder to write text continuously
                placeholder = st.empty()
                full_text = ""
                
                # Stream the agent's execution using 'updates' for perfectly stable chunking
                for update in st.session_state.agent.stream(
                    {"messages": [("user", prompt)]},
                    stream_mode="updates"
                ):
                    if "agent" in update:
                        messages = update["agent"].get("messages", [])
                        if messages:
                            # Get the completed message for this chunk
                            last_msg = messages[-1]
                            
                            if hasattr(last_msg, 'content') and last_msg.content:
                                chunk_text = ""
                                if isinstance(last_msg.content, str):
                                    chunk_text = last_msg.content
                                elif isinstance(last_msg.content, list):
                                    chunk_text = "".join(block.get("text", "") for block in last_msg.content if isinstance(block, dict))
                                
                                if chunk_text:
                                    full_text += chunk_text + "\n\n"
                                    # Update the UI instantly in real-time
                                    placeholder.markdown(full_text)
                
                # Final render without the blinking cursor
                placeholder.markdown(full_text)
                
                # Save to history
                st.session_state.messages.append({"role": "assistant", "content": full_text})
                
            except Exception as e:
                st.error(f"An error occurred: {e}")
