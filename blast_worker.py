from celery import Celery
from Bio.Blast import NCBIWWW, NCBIXML
import json

celery = Celery("blast_worker", broker="redis://localhost:6379/0")

@celery.task
def blast_task(seq, job_id):
    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", seq)
        blast_record = NCBIXML.read(result_handle)
        top_hit = blast_record.alignments[0].hit_def if blast_record.alignments else "No hits found"
    except Exception as e:
        top_hit = f"Error: {str(e)}"

    with open("blast_cache.json", "r+") as f:
        data = json.load(f)
        data[job_id] = top_hit
        f.seek(0)
        json.dump(data, f)
        f.truncate()
