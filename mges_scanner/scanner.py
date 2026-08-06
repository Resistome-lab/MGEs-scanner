import importlib.resources
import os
import pyhmmer
import pyrodigal

def get_hmm_db_path():
    """Dynamically locate the bundled HMM database."""
    try:
        path = importlib.resources.files("mges_scanner.db").joinpath("mge_signatures.hmm")
        return str(path)
    except AttributeError:
        with importlib.resources.path("mges_scanner.db", "mge_signatures.hmm") as p:
            return str(p)

def scan_contigs(fasta_path, evalue_cutoff=1e-5):
    """Scan assembly contigs for ORFs and match against MGE HMM profiles."""
    db_path = get_hmm_db_path()
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"HMM database not found at {db_path}")

    gene_finder = pyrodigal.GeneFinder()
    records = []

    # Read FASTA records
    with open(fasta_path, "r") as f:
        current_header = None
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header:
                    records.append((current_header, "".join(current_seq)))
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            records.append((current_header, "".join(current_seq)))

    results = []

    # Load pre-pressed HMM database
    with pyhmmer.plan7.HMMFile(db_path) as hmm_file:
        hmms = list(hmm_file)

    for contig_id, seq in records:
        seq_bytes = seq.encode("ascii")
        
        # Try finding genes using meta mode (ideal for contigs/metagenomes)
        # Fallback to training on sequence if single contig is sufficiently long
        try:
            preds = gene_finder.find_genes(seq_bytes, meta=True)
        except Exception:
            try:
                training_info = gene_finder.train(seq_bytes)
                preds = gene_finder.find_genes(seq_bytes, training_info)
            except Exception:
                preds = []

        proteins = [
            pyhmmer.easel.DigitalSequence(
                pyhmmer.easel.Alphabet.amino(),
                name=f"{contig_id}_{i+1}".encode("utf-8"),
                sequence=pred.translate().encode("utf-8")
            )
            for i, pred in enumerate(preds)
        ]

        detected_mges = set()
        if proteins:
            for top_hits in pyhmmer.hmmsearch(hmms, proteins):
                for hit in top_hits:
                    if hit.evalue <= evalue_cutoff:
                        detected_mges.add(top_hits.query_name.decode("utf-8"))

        mge_count = len(detected_mges)
        if mge_count >= 2:
            mrs, tier = 8.5, "HIGH"
        elif mge_count == 1:
            mrs, tier = 5.0, "MEDIUM"
        else:
            mrs, tier = 0.5, "LOW"

        mge_str = ";".join(sorted(detected_mges)) if detected_mges else "None"
        results.append({
            "Contig_ID": contig_id,
            "Detected_MGE_Signatures": mge_str,
            "Mobility_Risk_Score": mrs,
            "Risk_Tier": tier
        })

    return results
