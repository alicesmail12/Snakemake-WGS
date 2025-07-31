# Snakemake-WGS
Here I have created a snakemake workflow that performs WGS on several samples at once from FASTQ files. It can take a really long time, depending on HPC capacity!

### Steps
**Step 1: FastQC**

**FastQC** analyses FASTQ files and creates a **HTML quality check report**. This includes quality scores for each read and base position, as well as assessments of GC content, nucleotide distribution, read length, overrepresented sequences and adapter content. 

In a snakemake workflow on a HPC this is how FastQC can be run:
```python
# Modules
module load FastQC

# Run
fastqc {File}_R{Read}.fastq -o {FASTQC_DIR}
```

**Step 2: BWA Alignment**

**BWA mem** takes a reference genome and **aligns reads** from a FASTQ file using the Burrows-Wheeler Aligner algorithm. **SAMtools sort** then sorts the resulting **BAM file** by coordinate, and **SAMtools index** indexes the output BAM file to generate a **BAI file**. 

**SAMtools flagstat** and **SAMtools idxstats** can also be called to get a summary of the BWA alignment (for example, how many **total reads were aligned**, and how many were aligned to each chromosome).

```python
# Modules
module load ncurses SAMtools BWA picard

# Run
bwa mem -M -t 12 {input.fasta} {input.R1} {input.R2} | samtools sort - -O bam | tee {output.BAM} | samtools index - {output.BAI}
samtools flagstat {output.BAM} > {output.stat} 
samtools idxstats {output.BAM} > {output.idx} 
```

**Step 3: Removing Duplicates**

**PICARD MarkDuplicates** identifies and marks **duplicate reads**, and retains the read with the highest base quality scores. This is a way to try and correct for any sequencing errors. **SAMtools index** then indexes the output BAM file to generate a new **BAI file**.

```python
# Modules
module load SAMtools BWA picard

# Deduplicate
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I={input.BAM} \
    O={output.BAMdedup} \
    REMOVE_DUPLICATES=true \
    M={output.metrics}

# Index
samtools index {output.BAMdedup}
```

**Step 4: Fixing Read Groups**

I have run this in my workflow so that I can assign sample names to every read in a BAM file. This means I can make a group VCF file later on. **PICARD AddOrReplaceReadGroups** assigns RG values to all reads in a file.

```python
# Modules
module load Java picard

# Read groups
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    I={input.BAMdedup} \
    O={output.BAMRG} \
    SORT_ORDER=coordinate \
    RGPL=illumina \
    RGSM={params.File} \
    CREATE_INDEX=True \
    RGPU={params.File} 
```

**Step 5: Haplotype Calling**

**GATK HaplotypeCaller** calls SNPs and indels from each sample and creates a **GVCF file** which can then be used by GenotypeGVCFs to genotype multiple samples at a time, in a single file. The GVCF output contains the genotype at every genomic location, not just where there is variation from the reference.

Here I have used the list flag (L) to generate a GVCF file for each chromosome, for every sample. This means this command can multi-thread and it runs a lot quicker.

```python
# Modules
module load Java/17.0.6

# Haplotype calling
/Software/gatk-4.6.2.0/gatk --java-options "-Xmx80g" HaplotypeCaller \
    -R {params.Index} \
    -I {input.BAMdedup} \
    -O {output} \
    -L {params.Chr} \
    -ERC GVCF \
    --do-not-run-physical-phasing true \
    --output-mode EMIT_VARIANTS_ONLY \
    --native-pair-hmm-threads 12
```

**Step 6: Combine GVCFs**

Not only have I got multiple samples, I also have created GVCFs for every chromosome in each sample: so I need to combine them using **CombineGVCFs**. First I make a list of all the files I want to combine.

```python
# Generate list
with open('{Run}_gvcfs.list', 'w') as outfile:
      outfile.write('\n'.join(str(i) for i in '{File}_{Chr}.single.g.vcf.gz'))
```
```bash
# Modules
module load picard Java

# Combine
/Software/gatk-4.6.2.0/gatk --java-options "-Xmx80g -XX:ParallelGCThreads=6 -Djava.io.tmpdir=./" CombineGVCFs \
    -R {params.Index} \
    --variant {input.GVCFList} \
    -O {output.GVCFCom}
```

**Step 7: Convert to VCF**
Now I have a massive GVCF file, I can use **GenotypeGVCFS** to convert it into a smaller VCF file that just contains locations where there is a variant in at least one sample. 

```python
# Modules
module load picard Java

# Genotype
/Software/gatk-4.6.2.0/gatk --java-options "-Xmx80g" GenotypeGVCFs \
    -R {params.Index} \
    -V {input.GVCF} \
    -O {output.VCF}
```




