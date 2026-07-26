# MGEs-scanner

Mobile Genetic Element Profiler & Mobility Risk Scoring Engine for Genomes and Metagenomic Contigs.

`MGEs-scanner` scans assemblies and contigs for transposons, integrons, insertion sequences, and plasmid replicon markers to calculate a quantitative **Mobility Risk Score (MRS)** for detected genomic loci.

---

## PART 1: INSTALLATION & SETUP

### Step 1: Install MGEs-scanner
Install `mges-scanner` directly from GitHub using `pip3` in your Linux/WSL environment:

```bash
pip3 install git+[https://github.com/Resistome-lab/MGEs-scanner.git](https://github.com/Resistome-lab/MGEs-scanner.git)

### Step 2: Install Prerequisites

MGEs-scanner requires Prodigal, HMMER, and DIAMOND. Install them via Conda:
```bash
conda install -c bioconda prodigal hmmer diamond
```

## PART 2: HOW TO USE IT
### Step 1: Prepare your MGE HMM Database

Prepare a concatenated HMM database file containing MGE profiles (e.g., ISfinder, INTEGRALL, Pfam MGEs):
```bash
hmmpress MGE_profiles.hmm
```
### Step 2: Run MGEs-scanner on your assemblies/contigs

To profile contigs (.fasta or .fa), pass your input assemblies and database paths to mges-scanner:
```bash
mges-scanner --input contigs.fasta --mge-db MGE_profiles.hmm --plsdb plsdb.dmnd --output-dir mge_results --threads 8
```
Note: You can adjust parameters and file paths according to your system environment.

### Step 3: Look for the output files

MGEs-scanner automatically generates results inside your specified output directory (--output-dir):

    MGEs_mobility_summary.tsv: Quantified mobility summary per contig/locus:

        Contig_ID: Identifier of the analyzed contig.

        Mobility_Risk_Score: Calculated score from 0.0 to 10.0.

        Risk_Category: Categorized risk level (HIGH, MEDIUM, or LOW).

        Detected_Elements_Count: Total count of distinct MGE features detected.

        MGE_Signatures: Semicolon-delimited list of identified elements.

    predicted_proteins.faa: Translated ORFs extracted by Prodigal.

    hmmer_results.tbl: Tabular output generated directly by HMMER.
