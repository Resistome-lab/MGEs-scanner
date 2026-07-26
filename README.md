\# MGEs-scanner (Mobile Genetic Element Profiler)



`MGEs-scanner` is an automated bioinformatics CLI tool designed to scan genomic assemblies or contigs for mobile genetic elements (MGEs)—such as transposons, insertion sequences (IS), integrons, and plasmid replicon markers—and assign a \*\*Mobility Risk Score (MRS)\*\* to detected loci.



\---



\## Why Context Matters



An Antimicrobial Resistance (AMR) gene residing on a bacterial chromosome poses a significantly lower risk of rapid horizontal spread than the exact same gene located on a conjugative plasmid flanked by transposases. 



`MGEs-scanner` quantifies this dissemination potential by assessing the genetic neighborhood around target sequences, helping researchers prioritize high-risk mobilizable AMR threats.



\---



\## Features \& Workflow



1\. \*\*ORF Prediction:\*\* Translates input nucleotide contigs/assemblies into open reading frames using \*\*Prodigal\*\*.

2\. \*\*HMM Profile Scanning:\*\* Runs \*\*HMMER (`hmmscan`)\*\* against curated MGE databases (ISfinder, INTEGRALL, Pfam MGE profiles, or MobileOG-db).

3\. \*\*Plasmid Marker Verification:\*\* Performs fast protein search alignment against \*\*PLSDB\*\* using \*\*DIAMOND\*\* (Optional).

4\. \*\*Mobility Risk Scoring Engine:\*\* Aggregates detected MGE features per locus and computes a weighted Mobility Risk Score (MRS) scaled from 0.0 to 10.0.



\---



\## Prerequisites \& Installation



\### 1. System Dependencies

Install required bioinformatics tools via Conda/Mamba:



```bash

conda install -c bioconda prodigal hmmer diamond

```



\### 2. Python Dependencies

```bash

pip install pandas biopython

```

