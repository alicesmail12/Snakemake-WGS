# Snakemake-WGS
Here I have created a snakemake workflow that performs WGS on several samples at once from fastq files.

### Steps
**1: FastQC**
'''
fastqc {input} -o {FASTQC_DIR}
'''

