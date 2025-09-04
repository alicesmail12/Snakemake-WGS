# Snakemake-WGS
Here I have created a snakemake workflow that performs WGS on several samples at once from FASTQ files. It can take a really long time, depending on HPC capacity! Here I have just listed what each step does. Variables in curly braces are specified in my snakemake rules (see `Snakefile.py`).

![alt text](https://github.com/alicesmail12/Snakemake-WGS/blob/main/WGS-Snakemake-DAG.png?raw=true)

The resulting file looks like this:

![alt text](https://github.com/alicesmail12/Snakemake-WGS/blob/main/WGS-Out-Table.png?raw=true)

### Steps
**Step 1: FastQC**

First I used `FastQC`, which analyses FASTQ files and creates a **HTML quality check report**, including quality scores for each read and base position, as well as assessments of GC content, nucleotide distribution, read length, overrepresented sequences and adapter content. 

In a snakemake workflow on a HPC this is how FastQC can be run:
```bash
# Modules
module load FastQC

# Run
fastqc {File}_R{Read}.fastq -o {FASTQC_DIR}
```

**Step 2: BWA Alignment**

Then I used `BWA mem`, which takes a reference genome and **aligns reads** from a FASTQ file using the Burrows-Wheeler Aligner algorithm. `SAMtools sort` then sorts the resulting **BAM file** by coordinate, and `SAMtools index` indexes the output BAM file to generate a **BAI file**. 

`SAMtools flagstat` and `SAMtools idxstats` can also be called to get a summary of the BWA alignment (for example, how many **total reads were aligned**, and how many were aligned to each chromosome).

```bash
# Modules
module load ncurses SAMtools BWA picard

# Run
bwa mem -M -t 12 {input.fasta} {input.R1} {input.R2} | samtools sort - -O bam | tee {output.BAM} | samtools index - {output.BAI}
samtools flagstat {output.BAM} > {output.stat} 
samtools idxstats {output.BAM} > {output.idx} 
```

**Step 3: Removing Duplicates**

`MarkDuplicates` identifies and marks **duplicate reads**, and retains the read with the highest base quality scores. This is a way to try and correct for any sequencing errors. `SAMtools index` then indexes the output BAM file to generate a new **BAI file**.

```bash
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

I have run `AddOrReplaceReadGroups` in my workflow so that I can assign sample names to every read in a BAM file. This means I can make a group VCF file later on. 

```bash
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

`GATK HaplotypeCaller` calls SNPs and indels from each sample and creates a **GVCF file** which can then be used by GenotypeGVCFs to genotype multiple samples at a time, in a single file. The GVCF output contains the genotype at **every genomic location**, not just where there is variation from the reference.

Here I have used the list flag (`L`) to generate a GVCF file for each chromosome, for every sample. This means this command can multi-thread and it runs a lot quicker.

```bash
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

Not only have I got multiple samples, I also have created GVCFs for every chromosome in each sample: so I need to **combine** them using `CombineGVCFs` to make a **really large GVCF file**. First I make a list of all the files I want to combine, and then input this into `CombineGVCFs` using the `--variant` flag.

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

Now I have a massive GVCF file, I can use `GenotypeGVCFS` to convert it into a smaller **VCF file** that just contains locations where there is a variant in at least one sample. 

```bash
# Modules
module load picard Java

# Genotype
/Software/gatk-4.6.2.0/gatk --java-options "-Xmx80g" GenotypeGVCFs \
    -R {params.Index} \
    -V {input.GVCF} \
    -O {output.VCF}
```

**Step 8: Split Variant Types**

From the vcf file, I can then use `SelectVariants` to **split the vcf** into indels and SNPs for later filtering.

```bash
# Modules
module load picard Java

# Get SNPs    
{params.GATKDir}/gatk SelectVariants \
    -V {input.VCF} \
    -select-type SNP \
    -O {output.VCFSNP}

# Get indels
{params.GATKDir}/gatk SelectVariants \
    -V {input.VCF} \
    -select-type INDEL \
    -O {output.VCFIndel}
```

**Step 9: Filter Variant Types**

With the split variant files I can apply different **quality control filters** using `VariantFiltration`.
- `QD` is a metric representing normalised variant quality. In the below code variants with a QD below 2 are marked with 'QD2'.
- `QUAL` is the variant confidence (QD is derived from QUAL).
- `SOR` is the Strand Odds Ratio, an estimation of strand bias that takes into account the ratio of reads covering two alleles.
- `FS` is called Fisher strand, representing the Phred scaled probability of strand bias.
- `MQ` gives a value for mapping quality.
- `MQRankSum` compares mapping quality scores that support the reference allele or the alternate allele.
- `ReadPosRankSum` compares whether the positions of the reference and alternate alleles are different across reads. A negative value means that the alternate allele is found at the ends of each read more often than the reference allele (and vice versa).

```bash
# Modules
module load picard Java

# Filtering
{params.GATKDir}/gatk VariantFiltration \
    -V {input.VCFSNP} \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O {output.VCFSNP}
    
{params.GATKDir}/gatk VariantFiltration \
    -V {input.VCFIndel} \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"  \
    -O {output.VCFIndel}

# Merge
java -jar $EBROOTPICARD/picard.jar MergeVcfs \
    I={output.VCFSNP} \
    I={output.VCFIndel} \
    O={output.VCF}
```

**Step 10: Normalise VCF**

`BCFTools norm` left-aligns the variants and normalises indels, as well as checking if reference alleles match. Here I add `-m -` to split sites that are multi-allelic into rows that are biallelic, and `-Oz` to specify that the output is a compressed vcf. I can also run `BCFTools stats` to get some information about the vcf file, such as non-reference allele frequency, depth and read quality.

```bash
# Modules
module load BCFtools texlive

# Normalise
bcftools norm -m - -Oz {input.VCF} > {output.VCF}

# Stats
bcftools stats -s - {input.VCF} > {output.VCFStats}
```

**Step 11: VEP**

Here I use `VEP` locally to gather **extra information** about each variant, such as the nearest gene, PolyPhen/Sift/CADD scores and gnomAD allele frequencies. To run this step I downloaded VEP (https://grch38.ensembl.org/info/docs/tools/vep/script/vep_download.html) and some files required for the plugins. VEP can take a really long time depending on how many annotations you want to make. 

```bash
# Modules
module load Perl tabix Bio-DB-HTS DBD-mysql OpenSSL

# Set directory to VEP
cd ./Software/ensembl-vep-release-114/
    
# Plugins
export PERL5LIB=$PERL5LIB:/Software/ensembl-vep-release-114/Plugins
    
# VEP
{params.VEPSoftware}/vep --cache \
    --dir_cache "./Software/ensembl-vep-release-114/.vep" \
    --input_file {input.VCF} \
    --output_file {output.VCF} \
    --fasta {params.Fasta} \
    --offline \
    --species homo_sapiens \
    --force_overwrite \
    --format vcf \
    --vcf \
    --no_check_variants_order \
    --check_existing \
    --freq_pop gnomAD \
    --assembly GRCh38 \
    --hgvs \
    --variant_class \
    --keep_csq \
    --af_gnomad \
    --polyphen b \
    --sift b \
    --symbol \
    --total_length \
    --plugin LoFtool \
    --plugin CADD,/Software/ensembl-vep-release-114/Plugins/cadd/whole_genome_SNVs.tsv.gz,/Software/ensembl-vep-release-114/Plugins/cadd/gnomad.genomes.r3.0.indel.tsv.gz \
    --fields "Uploaded_variation,Location,Allele,Gene,Feature,SYMBOL,Existing_variation,VARIANT_CLASS,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,HGVSc,HGVSp,BIOTYPE,IMPACT,CLIN_SIG,PolyPhen,SIFT,MAX_AF,gnomAD_AF,AF,CADD_PHRED,CADD_RAW" \
    --max_af \
    --pick \
    --pick_order rank,canonical,tsl \
    --buffer_size 50000 \
    --fork 8
```

**Step 12: Format VEP Output**

I have a **Python script** that takes in each VCF line and formats it, generating a **csv file** with 'HET', 'REF', 'HOM' or 'unknown' for each variant for each sample. The variables I pass to the Python file include the VCF file, the sample names (what I want the output to be called), an OMIM file (to annotate the VCF with OMIM entries), and the working directory where the output will be created.

```bash
# Modules
module load Python zlib pandas
    
# Python
python {params.FormatPy} --VCFFile {input.VCF} \
    --SampleName {params.Run} \
    --OMIMFile {params.OMIM} \
    --WDir {params.VEP_Dir}
```

Within the Python file I can apply several **filters** (such as for high CADD scores or low AFs).

```python
# CADD
vcfCombinedFilt = vcfCombinedFilt.loc[((vcfCombinedFilt['CADD_PHRED'])>11)]

# Max AF
vcfCombinedFilt = vcfCombinedFilt.loc[((vcfCombinedFilt['MAX_AF'])<0.01)]
```











