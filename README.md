# MGEs-scanner

Mobile Genetic Element Profiler & Mobility Risk Scoring Engine for Genomes and Metagenomic Contigs.

`MGEs-scanner` scans assemblies and contigs for transposons, integrons, insertion sequences, and plasmid replicon markers to calculate a quantitative **Mobility Risk Score (MRS)** for detected genomic loci.

---

## PART 1: INSTALLATION & SETUP

### Install MGEs-scanner
Install `mges-scanner` directly from GitHub using `pip3` in your Linux/WSL environment:

```bash
pip3 install git+https://github.com/Resistome-lab/MGEs-scanner.git
```
Note: mges-scanner comes pre-packaged with a bundled HMM signature database (mge_signatures.hmm), so no additional database preparation, Conda packages, or profile downloads are required.

## PART 2: HOW TO USE IT

### Step 1: Run mges-scanner on assembly contigs

Pass your input assemblies (.fasta, .fa, or .fna) to mges-scanner:
```bash
mges-scanner -i contigs.fasta -o mge_results.tsv
```
You can also customize the E-value threshold:
```bash
mges-scanner --input contigs.fasta --output mge_results.tsv --evalue 1e-10
```

### Step 2: Run MGEs-scanner on your assemblies/contigs

To profile contigs (.fasta or .fa), pass your input assemblies and database paths to mges-scanner:
```bash
mges-scanner --input contigs.fasta --mge-db MGE_profiles.hmm --plsdb plsdb.dmnd --output-dir mge_results --threads 8
```
Note: You can adjust parameters and file paths according to your system environment.

### Step 3: Look for the output files

mges-scanner generates a single tab-delimited summary report (.tsv) containing mobility classifications for each analyzed contig:

1) Contig_ID: Identifier of the analyzed contig/assembly fragment.
2) Detected_MGE_Signatures: Semicolon-delimited list of identified Pfam MGE signatures (or None).
3) Mobility_Risk_Score: Quantitative Mobility Risk Score (0.5 to 8.5+).
4) Risk_Tier: Categorized mobility risk level (HIGH, MEDIUM, or LOW).
