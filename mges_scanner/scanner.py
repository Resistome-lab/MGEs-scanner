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

    # Step 1: Predict ORFs using Pyrodigal
    orf_finder = pyrodigal.OrfFinder()
    records = []
    
    # Simple FASTA parser
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
    
    # Step 2: HMM Scanning with PyHMMER
    with pyhmmer.plan7.HMMFile(db_path) as hmm_file:
        hmms = list(hmm_file)

    for contig_id, seq in records:
        preds = orf_finder.find_genes(seq)
        proteins = [
            pyhmmer.easel.DigitalSequence(
                pyhmmer.easel.Alphabet.amino(),
                name=f"{contig_id}_{i+1}".encode(),
                sequence=pred.translate().encode()
            )
            for i, pred in enumerate(preds)
        ]
        
        detected_mges = set()
        if proteins:
            for top_hits in pyhmmer.hmmsearch(hmms, proteins):
                for hit in top_hits:
                    if hit.evalue <= evalue_cutoff:
                        detected_mges.add(top_hits.query_name.decode())

        # Step 3: Compute Mobility Risk Score (MRS)
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
