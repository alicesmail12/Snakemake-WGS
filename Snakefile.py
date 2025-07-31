# Filename wildcards
FILES = ['Sample1', 'Sample2'] 
READS = [1, 2]
RUN = ['Sample1-Sample2']

# Chromosomes
numbers = list(range(1, 23))
CHRS = ["chr" + str(num) for num in numbers]
CHRS.extend(['chrM', 'chrX', 'chrY'])

# Directories
CAT_FASTQ_DIR='/Fastq-Files/'
FASTQC_DIR='/FastQC-Output/'
BWA_DIR='/BWA-Output/'
SAMTOOLS_DIR='/Samtools-Output/'
VCF_DIR='/VCF-Output/'
VEP_DIR='/VEP-Output/'

# Rule all
rule all:
  input:
   expand(os.path.join(FASTQC_DIR, '{File}_R{Read}_fastqc.html'), File=FILES, Read=READS),
   expand(os.path.join(VEP_DIR, '{Run}-Formatted-Filtered.csv'), Run=RUN)

# Step 1: FastQC
rule Get_FastQC:
  input:
    os.path.join(CAT_FASTQ_DIR, '{File}_R{Read}.fastq')
  output:
    os.path.join(FASTQC_DIR, '{File}_R{Read}_fastqc.{Filetype}')
  shell:
    """
    #!/bin/bash
    echo "Running FastQC"
    module load FastQC/0.11.9-Java-1.8.0_144 | fastqc {input} -o {FASTQC_DIR}
    echo "FastQC completed"
    """

# Step 2: BWA Alignment
rule BWA_Align:
  input:
    fasta = 'hg38.fa',
    R1 = os.path.join(CAT_FASTQ_DIR, '{File}_R1.fastq'),
    R2 = os.path.join(CAT_FASTQ_DIR, '{File}_R2.fastq')
  output:
    BAM = os.path.join(BWA_DIR, '{File}_bwa_output.bam'),
    BAI = os.path.join(BWA_DIR, '{File}_bwa_output.bam.bai')
  threads: 12
  shell:
    """
    #!/bin/bash
    "Starting BWA alignment"
    module load ncurses/6.4-GCCcore-12.3.0 SAMtools BWA/0.7.17-foss-2017b picard/2.20.6-Java-1.8.0_144 
    bwa mem -M -t {threads} {input.fasta} {input.R1} {input.R2} | samtools sort - -O bam | tee {output.BAM} | samtools index - {output.BAI}
    echo "BWA alignment completed"
    """  
  
# Step 3: BWA Statistics
rule BWA_Stats:
  input:
    BAM = os.path.join(BWA_DIR, '{File}_bwa_output.bam')
  output:
    idx = os.path.join(SAMTOOLS_DIR, '{File}_idxstat_bwa_output.txt'),
    stat = os.path.join(SAMTOOLS_DIR, '{File}_stat_bwa_output.txt')
  shell:
    """
    #!/bin/bash
    echo "Generating BWA statistics"
    module load ncurses/6.4-GCCcore-12.3.0 SAMtools/1.17-GCC-12.2.0 | samtools flagstat {input.BAM} > {output.stat} 
    samtools idxstats {input.BAM} > {output.idx} 
    echo "BWA statistics generated"
    """

# Step 4: Remove duplicate reads
rule DeDup:
  input:
    BAM = os.path.join(BWA_DIR, '{File}_bwa_output.bam')
  output:
    BAMdedup = os.path.join(BWA_DIR, '{File}_bwa_dedup.bam'),
    metrics = os.path.join(BWA_DIR, '{File}_dup_metrics.txt')
  shell:
    """
    #!/bin/bash
    module load SAMtools/1.17-GCC-12.2.0 BWA/0.7.17-foss-2017b picard/2.20.6-Java-1.8.0_144
    echo "Removing duplicates using Picard"
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I={input.BAM} \
    O={output.BAMdedup} \
    REMOVE_DUPLICATES=true \
    M={output.metrics}
    samtools index {output.BAMdedup}
    echo "Duplicates removed"
    """
    
# Step 5: Fix read groups
rule RG:
  input:
    BAMdedup = os.path.join(BWA_DIR, '{File}_bwa_dedup.bam')
  params:
    File=lambda wc: wc.get("File")
  output:
    BAMRG = os.path.join(BWA_DIR, '{File}_bwa_dedup_with_RG.bam')
  shell:
    """
    #!/bin/bash
    echo 'Fixing read groups'
    module load Java/17.0.6 picard/2.20.6-Java-1.8.0_144
    java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    I={input.BAMdedup} \
    O={output.BAMRG} \
    SORT_ORDER=coordinate \
    RGID=foo \
    RGLB=bar \
    RGPL=illumina \
    RGSM={params.File} \
    CREATE_INDEX=True \
    RGPU={params.File} 
    echo 'Read groups fixed'
    """

# Step 6: Haplotype caller
rule HapCaller:
  input:
    BAMdedup = os.path.join(BWA_DIR, '{File}_bwa_dedup_with_RG.bam')
  params:
    GATKDir = "/Software/gatk-4.6.2.0",
    Index = "hg38.fa",
    Chr = lambda wc: wc.get("Chr")
  resources:
    threads=12
  threads: 1
  output:
    os.path.join(VCF_DIR, '{File}_{Chr}.single.g.vcf.gz')
  shell:
    """
    #!/bin/bash
    module load Java/17.0.6 
    echo "Running GATK"
    {params.GATKDir}/gatk --java-options "-Xmx80g" HaplotypeCaller \
    -R {params.Index} \
    -I {input.BAMdedup} \
    -O {output} \
    -L {params.Chr} \
    -ERC GVCF \
    --do-not-run-physical-phasing true \
    --output-mode EMIT_VARIANTS_ONLY \
    --native-pair-hmm-threads {threads}
    echo "GATK Finished"
    """  

# Step 7: Make GVCF list
rule GVCFList:
  input:
    expand(os.path.join(VCF_DIR, '{File}_{Chr}.single.g.vcf.gz'), Chr=CHRS, File=FILES)
  output:
    expand(os.path.join(VCF_DIR, '{Run}_gvcfs.list'), Run=RUN)
  run:
    with open(output[0], 'w') as outfile:
      outfile.write('\n'.join(str(i) for i in input))

# Step 8: Combine CHR VCFS
rule CombineGVCFS:
  input:
    GVCFList = os.path.join(VCF_DIR, '{Run}_gvcfs.list')
  params:
    GATKDir = "/Software/gatk-4.6.2.0",
    Index = "hg38.fa"
  resources:
    threads=12
  output:
    GVCFCom = os.path.join(VCF_DIR, '{Run}.combined.g.vcf.gz')
  shell:
    """
    #!/bin/bash
    module load picard/2.20.6-Java-1.8.0_144 Java/17.0.6
    echo 'Combining gVCFs'
    {params.GATKDir}/gatk --java-options "-Xmx80g -XX:ParallelGCThreads=6 -Djava.io.tmpdir=./" CombineGVCFs \
    -R {params.Index} \
    --variant {input.GVCFList} \
    -O {output.GVCFCom}
    """
  
# Step 9: Convert to VCF
rule VCFConvert:
  input:
    GVCF = os.path.join(VCF_DIR, '{Run}.combined.g.vcf.gz')
  params:
    GATKDir = /Software/gatk-4.6.2.0",
    Index = "hg38.fa"
  output:
    VCF = os.path.join(VCF_DIR, '{Run}.converted.vcf.gz')
  shell:
    """
    #!/bin/bash
    module load picard/2.20.6-Java-1.8.0_144 Java/17.0.6
    echo 'Calling genotypes'
    {params.GATKDir}/gatk --java-options "-Xmx80g" GenotypeGVCFs \
    -R {params.Index} \
    -V {input.GVCF} \
    -O {output.VCF}
    """
    
# Step 10: Split variant types
rule SplitVCF:
  input:
    VCF = os.path.join(VCF_DIR, '{Run}.converted.vcf.gz')
  params:
    GATKDir = "/Software/gatk-4.6.2.0"
  output:
    VCFSNP = os.path.join(VCF_DIR, '{Run}_SNPs.vcf.gz'),
    VCFIndel = os.path.join(VCF_DIR, '{Run}_indels.vcf.gz')
  shell:
    """
    #!/bin/bash
    module load picard/2.20.6-Java-1.8.0_144 Java/17.0.6
    
    {params.GATKDir}/gatk SelectVariants \
    -V {input.VCF} \
    -select-type SNP \
    -O {output.VCFSNP}
    
    {params.GATKDir}/gatk SelectVariants \
    -V {input.VCF} \
    -select-type INDEL \
    -O {output.VCFIndel}
    """

# Step 11: Filter variants
rule FilterVCF:
  input:
    VCFSNP = os.path.join(VCF_DIR, '{Run}_SNPs.vcf.gz'),
    VCFIndel = os.path.join(VCF_DIR, '{Run}_indels.vcf.gz')
  params:
    GATKDir = "/Software/gatk-4.6.2.0"
  output:
    VCFSNP = os.path.join(VCF_DIR, '{Run}_SNPs_Filtered.vcf.gz'),
    VCFIndel = os.path.join(VCF_DIR, '{Run}_indels_Filtered.vcf.gz'),
    VCF = os.path.join(VCF_DIR, '{Run}_Filtered.vcf.gz')
  shell:
    """
    #!/bin/bash
    module load picard/2.20.6-Java-1.8.0_144 Java/17.0.6
    
    echo 'Filtering variants'
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
    
    module load picard/2.20.6-Java-1.8.0_144 Java/17.0.6
    java -jar $EBROOTPICARD/picard.jar MergeVcfs \
    I={output.VCFSNP} \
    I={output.VCFIndel} \
    O={output.VCF}
    """

# Step 12: Normalise VCF
rule NormaliseVCF:
  input:
    VCF = os.path.join(VCF_DIR, '{Run}_Filtered.vcf.gz')
  output:
    VCF = os.path.join(VCF_DIR, '{Run}_Filtered_BCFTools.vcf.gz')
  shell:
    """
    #!/bin/bash
    module load BCFtools/1.10.2-foss-2019b
    bcftools norm -m - -Oz {input.VCF} > {output.VCF}
    """

# Step 13: VCF stats
rule VCFStats:
  input:
    VCF = os.path.join(VCF_DIR, '{Run}_Filtered_BCFTools.vcf.gz')
  output:
    VCFStats = os.path.join(VCF_DIR, '{Run}_vcf_stats.vchk')
  shell:
    """
    #!/bin/bash
    module load BCFtools/1.10.2-foss-2019b texlive/20200406-GCCcore-10.2.0
    bcftools stats -s - {input.VCF} > {output.VCFStats}
    """

# Step 14: VEP
rule VEP:
  input:
    VCF = os.path.join(VCF_DIR, '{Run}_Filtered_BCFTools.vcf.gz')
  params:
    VEPSoftware = "/Software/ensembl-vep-release-114",
    Fasta = "hg38.fa"
  resources:
    threads=12
  output:
    VCF = os.path.join(VEP_DIR, '{Run}_VEP.vcf'),
    stat = os.path.join(VEP_DIR, '{Run}_VEP.vcf_summary.html')
  shell:
    """
    #!/bin/bash
    echo 'Loading modules...'
    module load Perl/5.34.0-GCCcore-11.2.0 tabix/0.2.6-GCCcore-10.2.0 Bio-DB-HTS/3.01-GCC-11.2.0 DBD-mysql/4.050-GCC-11.2.0 OpenSSL/1.1
    cd /Software/ensembl-vep-release-114/
    
    # Plugins
    export PERL5LIB=$PERL5LIB:/Software/ensembl-vep-release-114/Plugins
    
    # VEP
    echo 'VEPing...'
    {params.VEPSoftware}/vep --cache \
    --dir_cache "/Software/ensembl-vep-release-114/.vep" \
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
    """

# Step 15: Format VEP Output
rule VEPFormat:
  input:
    VCF = os.path.join(VEP_DIR, '{Run}_VEP.vcf')
  params:
    FormatPy = "Format-WGS-Output.py",
    OMIM = "Full-OMIM-Table.txt",
    VEP_Dir = "/VEP-Output/",
    Run = lambda wc: wc.get("Run")
  output:
    Excel = os.path.join(VEP_DIR, '{Run}-Formatted-Filtered.csv')
  shell:
    """
    #!/bin/bash
    echo 'Formatting VEP File...'
    module load Python/3.9.6-GCCcore-11.2.0 zlib/1.2.11-GCCcore-11.2.0 pandas/1.1.2-foss-2020a-Python-3.8.2
    
    # Python
    python {params.FormatPy} --VCFFile {input.VCF} \
    --SampleName {params.Run} \
    --OMIMFile {params.OMIM} \
    --WDir {params.VEP_Dir}
    """
