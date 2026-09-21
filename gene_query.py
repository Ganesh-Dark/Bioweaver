import sys
sys.dont_write_bytecode = True
from agent_graph import build_ncbi_kegg_graph
import argparse



def run_gene_query_agent(user_query):
    # user_query: The natural language question string asked by the user.

    print(f"Starting LangGraph Gene Query Agent for query: '{user_query}'")

    agent = build_ncbi_kegg_graph()

    print("\n--- [AGENT INTERNAL MONOLOGUE] ---")

    print("\n" + "=" * 50)
    print("FINAL LANGGRAPH SYNTHESIZED RESPONSE")
    print("=" * 50)

    # Stream the agent's execution step by step using updates mode
    for update in agent.stream(
        {"messages": [("user", user_query)]},
        stream_mode="updates"
    ):
        if "agent" in update:
            messages = update["agent"].get("messages", [])
            if messages:
                last_msg = messages[-1]
                
                # Print tool call decisions from the AI to the monologue
                if hasattr(last_msg, 'tool_calls') and last_msg.tool_calls:
                    for call in last_msg.tool_calls:
                        print(f"AI Decision: Calling tool '{call['name']}' with args: {call['args']}")

                # Print the final answer text
                if hasattr(last_msg, 'content') and last_msg.content:
                    chunk_text = ""
                    if isinstance(last_msg.content, str):
                        chunk_text = last_msg.content
                    elif isinstance(last_msg.content, list):
                        chunk_text = "".join(block.get("text", "") for block in last_msg.content if isinstance(block, dict))
                    
                    if chunk_text:
                        print(chunk_text + "\n", end="", flush=True)

        if "tools" in update:
            messages = update["tools"].get("messages", [])
            if messages:
                last_msg = messages[-1]
                if hasattr(last_msg, 'name') and last_msg.name:
                    print(f"Tool '{last_msg.name}' completed and returned data to the AI.")

    print("\n" + "=" * 50 + "\n")


parser = argparse.ArgumentParser(description="Run Your Disease related Query")
parser.add_argument(
    "-q",
    "--query",
    nargs="+",
    required=True,
    help="The natural language query to ask the agent."
)
args = parser.parse_args()

if args.query:
    user_query = " ".join(args.query)
    run_gene_query_agent(user_query)
else:
    parser.print_help()

# Detailed Script Workings:
# 1. run_gene_query_agent builds the ReAct tool-calling agent from agent_graph.py.
# 2. It wraps the user query in a messages list and invokes the agent.
# 3. The agent runs the Reasoning->Acting->Observation loop internally until it has an answer.
# 4. The final message (the agent's answer) is extracted and printed to the terminal.
