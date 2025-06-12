from flask import Flask, request, render_template
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import primer3

app = Flask(__name__)
Entrez.email = "aditya@retrobiotech.in"  # Use your actual email

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

def run_blast(seq):
    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", seq)
        blast_record = NCBIXML.read(result_handle)
        return blast_record.alignments[0].hit_def if blast_record.alignments else "No hits found"
    except Exception as e:
        return f"Error: {str(e)}"

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

        result["Forward"] = {"seq": fwd, "blast": run_blast(fwd)}
        result["Reverse"] = {"seq": rev, "blast": run_blast(rev)}
        result["Probe"] = {"seq": probe, "blast": run_blast(probe)}

        return render_template("index.html", result=result, name=name)

    return render_template("index.html", result=None)

if __name__ == "__main__":
    import os
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port)
