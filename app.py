from flask import Flask, request, render_template, jsonify
from Bio import Entrez, SeqIO
from uuid import uuid4
import primer3
import json
import os
from blast_worker import blast_task

app = Flask(__name__)
BLAST_RESULT_FILE = "blast_cache.json"
Entrez.email = "your_email@example.com"  # Required by NCBI

# Ensure BLAST result cache
if not os.path.exists(BLAST_RESULT_FILE):
    with open(BLAST_RESULT_FILE, "w") as f:
        json.dump({}, f)

def save_blast_result(job_id, result):
    with open(BLAST_RESULT_FILE, "r+") as f:
        data = json.load(f)
        data[job_id] = result
        f.seek(0)
        json.dump(data, f)
        f.truncate()

def get_blast_result(job_id):
    with open(BLAST_RESULT_FILE) as f:
        data = json.load(f)
    return data.get(job_id)

def fetch_sequence_from_ncbi(term):
    try:
        handle = Entrez.esearch(db="nucleotide", term=term, retmode="xml")
        record = Entrez.read(handle)
        ids = record["IdList"]
        if not ids:
            return None
        fetch_handle = Entrez.efetch(db="nucleotide", id=ids[0], rettype="fasta", retmode="text")
        seq_record = SeqIO.read(fetch_handle, "fasta")
        return str(seq_record.seq)
    except Exception as e:
        return None

@app.route("/", methods=["GET", "POST"])
def index():
    result = {}
    if request.method == "POST":
        name = request.form["name"].strip()
        sequence = fetch_sequence_from_ncbi(name)
        if not sequence:
            return render_template("index.html", error="Could not fetch sequence for: " + name)

        primers = primer3.bindings.designPrimers(
            {"SEQUENCE_TEMPLATE": sequence},
            {
                'PRIMER_OPT_SIZE': 20,
                'PRIMER_PRODUCT_SIZE_RANGE': [[80, 150]],
                'PRIMER_PICK_INTERNAL_OLIGO': 1,
                'PRIMER_NUM_RETURN': 1
            }
        )

        fwd = primers["PRIMER_LEFT_0_SEQUENCE"]
        rev = primers["PRIMER_RIGHT_0_SEQUENCE"]
        probe = primers["PRIMER_INTERNAL_0_SEQUENCE"]

        result["primers"] = {
            "Forward": {"seq": fwd, "job_id": str(uuid4())},
            "Reverse": {"seq": rev, "job_id": str(uuid4())},
            "Probe": {"seq": probe, "job_id": str(uuid4())}
        }

        # Launch BLAST jobs
        for primer in result["primers"].values():
            blast_task.apply_async(args=[primer["seq"], primer["job_id"]])

        return render_template("index.html", result=result, name=name)

    return render_template("index.html", result=None)

@app.route("/blast_result/<job_id>")
def blast_result(job_id):
    result = get_blast_result(job_id)
    if result:
        return jsonify({"ready": True, "result": result})
    return jsonify({"ready": False})

if __name__ == "__main__":
    import os
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port)
