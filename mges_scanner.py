#!/usr/bin/env python3

"""
MGEs-scanner: Mobile Genetic Element Profiler
Scans assemblies/contigs for transposons, integrons, and plasmid replicon markers 
to assign a Mobility Risk Score (MRS) to detected genomic loci.
"""

import sys
import os
import argparse
import subprocess
import pandas as pd
from pathlib import Path

# --- Mobility Risk Scoring Matrix ---
FEATURE_WEIGHTS = {
    "conjugative_plasmid": 5.0,   # Conjugative plasmid / Inc group
    "integron_intI": 4.5,         # Integron integrase (intI1, intI2, etc.)
    "transposase_IS": 3.5,        # Insertion Sequences / Transposases
    "plasmid_non_conj": 2.5,      # Non-conjugative plasmid replicon
    "chromosomal_marker": 0.5     # Core chromosomal signatures
}

def parse_args():
    parser = argparse.ArgumentParser(
        description="MGEs-scanner: Profile Mobile Genetic Elements and calculate Mobility Risk Scores."
    )
    parser.add_argument("-i", "--input", required=True, help="Input contigs or assemblies (FASTA format)")
    parser.add_argument("-m", "--mge-db", required=True, help="Path to MGE HMM database (.hmm file)")
    parser.add_argument("-p", "--plsdb", help="Path to indexed PLSDB DIAMOND database (.dmnd)", default=None)
    parser.add_argument("-o", "--output-dir", default="mge_scanner_output", help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of CPU threads")
    parser.add_argument("-w", "--window", type=int, default=10000, help="Window size (bp) to calculate AMR-MGE proximity")
    return parser.parse_args()

def run_command(cmd, log_file=None):
    """Executes system commands with safety checks."""
    try:
        if log_file:
            with open(log_file, "w") as f:
                subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True)
        else:
            subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        print(f"[ERROR] Command failed: {' '.join(cmd)}", file=sys.stderr)
        sys.exit(1)

def predict_orfs(input_fasta, output_dir, threads):
    """Step 1: Predict ORFs using Prodigal."""
    print("[1/4] Predicting Open Reading Frames (ORFs) with Prodigal...")
    proteins = os.path.join(output_dir, "predicted_proteins.faa")
    gff = os.path.join(output_dir, "orfs.gff")
    
    cmd = ["prodigal", "-i", input_fasta, "-a", proteins, "-f", "gff", "-o", gff, "-q"]
    run_command(cmd)
    return proteins, gff

def scan_mge_hmms(protein_fasta, hmm_db, output_dir, threads):
    """Step 2: Profile MGE markers using HMMER (hmmscan)."""
    print("[2/4] Scanning protein sequences against MGE HMM profiles...")
    hmm_out = os.path.join(output_dir, "hmmer_results.tbl")
    
    cmd = [
        "hmmscan",
        "--cpu", str(threads),
        "--tblout", hmm_out,
        "-E", "1e-5",
        hmm_db,
        protein_fasta
    ]
    run_command(cmd, log_file=os.path.join(output_dir, "hmmer.log"))
    return hmm_out

def parse_hmmer_results(hmm_tbl):
    """Parses HMMER tabular output into a structured DataFrame."""
    results = []
    if not os.path.exists(hmm_tbl):
        return pd.DataFrame()

    with open(hmm_tbl, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) >= 15:
                target_name = parts[0]
                query_name = parts[2]
                e_value = float(parts[4])
                score = float(parts[5])
                
                results.append({
                    "ORF_ID": query_name,
                    "MGE_Hit": target_name,
                    "E_Value": e_value,
                    "BitScore": score
                })
    return pd.DataFrame(results)

def scan_plsdb(input_fasta, plsdb_path, output_dir, threads):
    """Step 3: Align against PLSDB using DIAMOND (Optional)."""
    if not plsdb_path or not os.path.exists(plsdb_path):
        print("[3/4] PLSDB database not provided or found. Skipping plasmid blast alignment.")
        return None
        
    print("[3/4] Aligning contigs against PLSDB for plasmid classification...")
    plsdb_out = os.path.join(output_dir, "plsdb_hits.tsv")
    
    cmd = [
        "diamond", "blastx",
        "-d", plsdb_path,
        "-q", input_fasta,
        "-o", plsdb_out,
        "-k", "1",
        "-f", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
        "-p", str(threads),
        "--quiet"
    ]
    run_command(cmd)
    
    if os.path.exists(plsdb_out) and os.path.getsize(plsdb_out) > 0:
        return pd.read_csv(plsdb_out, sep="\t", names=["Contig_ID", "PLSDB_Hit", "pident", "length", "evalue", "bitscore"])
    return None

def calculate_mobility_risk(df_mge, df_plsdb, output_dir):
    """Step 4: Combine hits and calculate Risk Score per Contig/Locus."""
    print("[4/4] Calculating Mobility Risk Scores (MRS)...")
    
    scores = {}
    
    # Process HMM hits
    if not df_mge.empty:
        for _, row in df_mge.iterrows():
            contig_id = "_".join(row["ORF_ID"].split("_")[:-1]) # Extract contig prefix
            hit_name = row["MGE_Hit"].lower()
            
            weight = FEATURE_WEIGHTS["transposase_IS"]
            if "intI" in hit_name or "integron" in hit_name:
                weight = FEATURE_WEIGHTS["integron_intI"]
            elif "inc" in hit_name or "rep" in hit_name:
                weight = FEATURE_WEIGHTS["plasmid_non_conj"]
            elif "tra" in hit_name or "trb" in hit_name or "virb" in hit_name:
                weight = FEATURE_WEIGHTS["conjugative_plasmid"]

            if contig_id not in scores:
                scores[contig_id] = {"Raw_Score": 0.0, "Elements": []}
            
            scores[contig_id]["Raw_Score"] += weight
            scores[contig_id]["Elements"].append(row["MGE_Hit"])

    # Process PLSDB hits
    if df_plsdb is not None and not df_plsdb.empty:
        for _, row in df_plsdb.iterrows():
            contig_id = row["Contig_ID"]
            if contig_id not in scores:
                scores[contig_id] = {"Raw_Score": 0.0, "Elements": []}
            scores[contig_id]["Raw_Score"] += FEATURE_WEIGHTS["conjugative_plasmid"]
            scores[contig_id]["Elements"].append(f"PLSDB:{row['PLSDB_Hit']}")

    # Build Summary DataFrame
    summary = []
    for contig, data in scores.items():
        raw_score = data["Raw_Score"]
        
        norm_score = min(10.0, round(raw_score, 2))
        if norm_score >= 7.0:
            risk_level = "HIGH (High Mobilizability Potential)"
        elif norm_score >= 3.5:
            risk_level = "MEDIUM (Moderate Mobilizability Potential)"
        else:
            risk_level = "LOW (Likely Chromosomal / Low Dissemination Risk)"

        summary.append({
            "Contig_ID": contig,
            "Mobility_Risk_Score": norm_score,
            "Risk_Category": risk_level,
            "Detected_Elements_Count": len(data["Elements"]),
            "MGE_Signatures": ";".join(set(data["Elements"]))
        })

    summary_df = pd.DataFrame(summary)
    
    if summary_df.empty:
        print("[INFO] No MGE signatures detected in input sequences.")
        summary_df = pd.DataFrame(columns=["Contig_ID", "Mobility_Risk_Score", "Risk_Category", "Detected_Elements_Count", "MGE_Signatures"])
        
    out_file = os.path.join(output_dir, "mge_mobility_summary.tsv")
    summary_df.to_csv(out_file, sep="\t", index=False)
    return out_file

def main():
    args = parse_args()
    
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("==================================================")
    print("      MGEs-scanner: Mobile Genetic Element Profiler")
    print(f" Input: {args.input}")
    print(f" Threads: {args.threads}")
    print(f" Output Directory: {out_dir.resolve()}")
    print("==================================================")

    # 1. ORF Prediction
    prot_fasta, gff_file = predict_orfs(args.input, out_dir, args.threads)
    
    # 2. HMMER Profiling
    hmm_tbl = scan_mge_hmms(prot_fasta, args.mge_db, out_dir, args.threads)
    df_mge = parse_hmmer_results(hmm_tbl)

    # 3. PLSDB Search
    df_plsdb = scan_plsdb(args.input, args.plsdb, out_dir, args.threads)

    # 4. Mobility Risk Scoring
    summary_tsv = calculate_mobility_risk(df_mge, df_plsdb, out_dir)

    print("==================================================")
    print(" Pipeline Completed Successfully!")
    print(f" Summary File: {summary_tsv}")
    print("==================================================")

if __name__ == "__main__":
    main()
