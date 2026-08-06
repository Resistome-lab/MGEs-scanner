import argparse
import sys
from mges_scanner.scanner import scan_contigs

def main():
    parser = argparse.ArgumentParser(
        description="MGEs-scanner: Mobile Genetic Element Profiler and Risk Scorer"
    )
    parser.add_argument(
        "--input", "-i", required=True, help="Path to input contigs FASTA file"
    )
    parser.add_argument(
        "--output", "-o", default="mge_summary.tsv", help="Output TSV path (default: mge_summary.tsv)"
    )
    parser.add_argument(
        "--evalue", "-e", type=float, default=1e-5, help="E-value cutoff threshold (default: 1e-5)"
    )

    args = parser.parse_args()

    print(f"[+] Scanning {args.input} for MGE signatures...")
    try:
        results = scan_contigs(args.input, evalue_cutoff=args.evalue)
    except Exception as e:
        print(f"[-] Error during scanning: {e}", file=sys.stderr)
        sys.exit(1)

    with open(args.output, "w") as f:
        f.write("Contig_ID\tDetected_MGE_Signatures\tMobility_Risk_Score\tRisk_Tier\n")
        for r in results:
            f.write(f"{r['Contig_ID']}\t{r['Detected_MGE_Signatures']}\t{r['Mobility_Risk_Score']}\t{r['Risk_Tier']}\n")

    print(f"[+] Scanning complete! Results saved to {args.output}")

if __name__ == "__main__":
    main()
