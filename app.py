from flask import Flask, request, render_template, jsonify
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import primer3

app = Flask(__name__)
Entrez.email = "aditya@retrobiotech.in"

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
            return render_template("index_fast.html", error="Could not fetch sequence", name=name)

        primers = primer3.bindings.designPrimers(
            {"SEQUENCE_TEMPLATE": sequence},
            {
                'PRIMER_OPT_SIZE': 20,
                'PRIMER_PRODUCT_SIZE_RANGE': [[80, 150]],
                'PRIMER_PICK_INTERNAL_OLIGO': 1,
                'PRIMER_NUM_RETURN': 1
            }
        )

        result = {
            "Forward": {"seq": primers["PRIMER_LEFT_0_SEQUENCE"]},
            "Reverse": {"seq": primers["PRIMER_RIGHT_0_SEQUENCE"]},
            "Probe": {"seq": primers["PRIMER_INTERNAL_0_SEQUENCE"]}
        }

        return render_template("index_fast.html", result=result, name=name)

    return render_template("index_fast.html")

@app.route("/blast")
def blast():
    from flask import request
    seq = request.args.get("seq")
    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", seq)
        blast_record = NCBIXML.read(result_handle)
        top_hit = blast_record.alignments[0].hit_def if blast_record.alignments else "No hits"
    except Exception as e:
        top_hit = f"Error: {str(e)}"
    return jsonify({"result": top_hit})

if __name__ == "__main__":
    import os
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port)
