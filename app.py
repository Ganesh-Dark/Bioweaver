import sys
sys.dont_write_bytecode = True

import os
from fastapi import FastAPI, Request
from fastapi.responses import StreamingResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
from agent_graph import build_ncbi_kegg_graph
from pydantic import BaseModel

app = FastAPI(title="Bioweaver API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount the static directory for the HTML frontend
base_dir = os.path.dirname(os.path.abspath(__file__))
static_dir = os.path.join(base_dir, "static")
if not os.path.exists(static_dir):
    os.makedirs(static_dir)

app.mount("/static", StaticFiles(directory=static_dir), name="static")

# Build the global agent instance
agent = build_ncbi_kegg_graph()

class ChatRequest(BaseModel):
    message: str

@app.get("/")
async def serve_frontend():
    return FileResponse(os.path.join(static_dir, "index.html"))

@app.post("/chat")
async def chat_endpoint(req: ChatRequest):
    """
    Streams the AI's response chunks back to the client using Server-Sent Events (SSE).
    """
    prompt = req.message

    async def stream_generator():
        try:
            for update in agent.stream(
                {"messages": [("user", prompt)]},
                stream_mode="updates"
            ):
                if "agent" in update:
                    messages = update["agent"].get("messages", [])
                    if messages:
                        last_msg = messages[-1]
                        if hasattr(last_msg, 'content') and last_msg.content:
                            chunk_text = ""
                            if isinstance(last_msg.content, str):
                                chunk_text = last_msg.content
                            elif isinstance(last_msg.content, list):
                                chunk_text = "".join(block.get("text", "") for block in last_msg.content if isinstance(block, dict))
                            
                            if chunk_text:
                                import json
                                # Use JSON encoding so newlines in markdown aren't broken by the SSE protocol
                                payload = json.dumps({"chunk": chunk_text})
                                yield f"data: {payload}\n\n"
        except Exception as e:
            import json
            yield f"data: {json.dumps({'error': str(e)})}\n\n"
            
    return StreamingResponse(stream_generator(), media_type="text/event-stream")

if __name__ == "__main__":
    import uvicorn
    # How to run: python app.py or uv run uvicorn app:app --reload
    print("Starting Bioweaver Backend Server...")
    uvicorn.run("app:app", host="0.0.0.0", port=8000, reload=True)
