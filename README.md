# Fully Automated Primer Designer with Async BLAST

Input a gene or pathogen name → Output primers, in silico analysis, and async BLAST.

## Features
- Gene name → NCBI nucleotide sequence
- Primer + probe generation
- GC%, Tm, hairpin/dimer check
- Async BLAST with Celery + Redis

## Setup

```bash
sudo apt update
sudo apt install redis-server -y
sudo systemctl enable redis --now
pip install -r requirements.txt
```

## Run

```bash
# Terminal 1
celery -A blast_worker.celery worker --loglevel=info

# Terminal 2
python app.py
```
