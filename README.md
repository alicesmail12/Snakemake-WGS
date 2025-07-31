# Snakemake-WGS
Here I have created a snakemake workflow that performs WGS on several samples at once from FASTQ files.

### Steps
**Step 1: FastQC**

FastQC analyses FASTQ files and creates a HTML quality check report. This includes quality scores for each read and base position, as well as assessments of GC content, nucleotide distribution, read length, overrepresented sequences and adapter content. In a snakemake workflow on a HPC this is how FastQC can be run:
```
module load FastQC
fastqc {input} -o {FASTQC_DIR}
```
Where input is `{File}_R{Read}.fastq` and the output is a directory.
