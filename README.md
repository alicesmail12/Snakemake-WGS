# Snakemake-WGS
Here I have created a snakemake workflow that performs WGS on several samples at once from FASTQ files. It can take a really long time, depending on HPC capacity!

### Steps
**Step 1: FastQC**

**FastQC** analyses FASTQ files and creates a **HTML quality check report**. This includes quality scores for each read and base position, as well as assessments of GC content, nucleotide distribution, read length, overrepresented sequences and adapter content. 

In a snakemake workflow on a HPC this is how FastQC can be run:
```
module load FastQC
fastqc {File}_R{Read}.fastq -o {FASTQC_DIR}
```

**Step 2: BWA Alignment**

**BWA mem** takes a reference genome and **aligns reads** from a FASTQ file using the Burrows-Wheeler Aligner algorithm. **SAMtools sort** then sorts the resulting **BAM file** by coordinate, and **SAMtools index** indexes the output BAM file to generate a **BAI file**. 

**SAMtools flagstat** and **SAMtools idxstats** can also be called to get a summary of the BWA alignment (for example, how many **total reads were aligned**, and how many were aligned to each chromosome).

```
module load ncurses SAMtools BWA picard
bwa mem -M -t 12 {input.fasta} {input.R1} {input.R2} | samtools sort - -O bam | tee {output.BAM} | samtools index - {output.BAI}
samtools flagstat {output.BAM} > {output.stat} 
samtools idxstats {output.BAM} > {output.idx} 
```

