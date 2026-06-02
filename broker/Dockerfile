# Hugging Face Docker Space for the upload broker (#43).
FROM python:3.12-slim

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Bundle the broker plus a copy of the lease verifier (single source of truth
# lives at screen/lease.py; keep this copy in sync — tests/test_lease.py pins
# the cross-language contract).
COPY app.py lease.py ./

# HF Spaces route to port 7860.
EXPOSE 7860
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "7860"]
