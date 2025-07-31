# Snakemake-WGS
Here I have created a snakemake workflow that performs WGS on several samples at once from FASTQ files.

### Steps
**Step 1: FastQC**
FastQC analyses FASTQ files and creates a HTML quality check report. This includes quality scores for each read and base position, as well as assessing GC content, nucleotide distribution, read length, overrepresented sequences and adapter content.
```
fastqc {input} -o {FASTQC_DIR}
```

