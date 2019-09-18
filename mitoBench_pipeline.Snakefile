################################################################################
# Snakefile for the processing of ancient DNA human mtDNA capture data
#
# Alex Huebner, 21/02/19
################################################################################

from snakemake.utils import min_version
from snakemake.utils import R

from glob import glob
import os
import pandas as pd

min_version("5.0")
shell.executable("/bin/bash")

workdir: config['tmpdir']

# Create directory for cluster file if it doesn't exist

if not os.path.isdir(f"{config['tmpdir']}/cluster_logs"):
    os.makedirs(f"{config['tmpdir']}/cluster_logs")

# Infer expected sample ids and paths to the corresponding BAM files from config 

SAMPLES = [line.rstrip() for line in open(config['samplelist'], "rt")]
BAMS = {sample: "{}/{}.{}".format(config['bamdir'], sample, config['bamsuffix'])
        for sample in SAMPLES}

# Extract project directory
PROJDIR = config['projdir']

# Auxilliary functions

def mixemt_downsampling(flagstatfn):
    ''' Determines the number of reads in a BAM file based on samtools flagstat
        and calculates the fraction of reads necessary to obtain 40,000 reads,
        the suggested input into mixEMT. The fraction will be used as input into
        samtools view -s for subsampling with SEED 0.
    '''
    if os.path.isfile(flagstatfn):
        with open(flagstatfn, "rt") as flagstatfile:
            nreads = int(next(flagstatfile).split(" ")[0])
        return "{:.4f}".format(30000 / nreads)
    else:
        return 1.0

# Snakemake rules
wildcard_constraints:
    sample = config['sampleIDconstraint']

localrules: link_index, seqdepth, contamMix_align_against_consensus, summary, copy_tmp_to_proj

rule all:
    input: expand("{projdir}/summary_table.csv", projdir=[PROJDIR]),
           expand("{projdir}/fasta/{sample}.fa", projdir=[PROJDIR], sample=SAMPLES)

# Prepare sorted, duplicate removed BAM files aligned only against the MT genome

rule link_index:
    # Link the raw data to file and index the BAM file in order to retrieve only
    # the MT reads
    output:
        bam = "rawdata/{sample}.bam",
        bai = "rawdata/{sample}.bam.bai"
    message: "Link the BAM file to the temp directory and index it: {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.2.2"
    params:
        bam = lambda wildcards: BAMS[wildcards.sample]
    threads: 1
    shell:
        """
        ln -s {params.bam} ${{PWD}}/rawdata/{wildcards.sample}.bam
        samtools index rawdata/{wildcards.sample}.bam
        """


rule extract_MT_reads:
    # Extract all reads that aligned against the MT genome without considering
    # the mapping quality
    input:
        bam = "rawdata/{sample}.bam",
        bai = "rawdata/{sample}.bam.bai"
    output:
        temp("bam/{sample}_MT.bam")
    message: "Extract the MT reads from sample {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.2.2"
    params: 
        region = "MT"
    threads: 1
    shell:
        """
        samtools view -bh \
                -o {output} {input.bam} {params.region}
        """


rule bwa_aln:
    input:
        "bam/{sample}_MT.bam"
    output:
        temp("bam/{sample}_MT.{i}.sai")
    message: "Align reads of type {wildcards.i} of sample {wildcards.sample} to only the MT genome with 1000bp overhang with MPI EVA BWA ancient settings (no seed)"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.2.2"
    params: reffa = f"{workflow.basedir}/resources/NC_012920_1000.fa"
    threads: 8
    shell:
        """
        bwa aln \
            -t {threads} \
            -n 0.01 -o 2 -l 16500 \
            -b -{wildcards.i} \
            -f {output} \
            {params.reffa} \
            {input}
        """


rule bwa_sampe:
    input:
        bam = "bam/{sample}_MT.bam",
        sais = lambda wildcards: [f"bam/{wildcards.sample}_MT.{i}.sai" for i in range(0, 3)]
    output:
        temp("bam/{sample}_MT_12.bam")
    message: "Generate BAM for non-merged reads of sample {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.2.2"
    params: reffa = f"{workflow.basedir}/resources/NC_012920_1000.fa",
            readgroup = lambda wildcards: r'@RG\tID:{sample}\tSM:{sample}\tPL:illumina'.format(sample = wildcards.sample)
    threads: 2
    shell:
        """
        bwa sampe \
            -r '{params.readgroup}' \
            -f /dev/stdout \
            {params.reffa} \
            {input.sais[1]} \
            {input.sais[2]} \
            {input.bam} \
            {input.bam} | \
        samtools view -hb - > {output}
        """


rule bwa_samse:
    input:
        bam = "bam/{sample}_MT.bam",
        sais = lambda wildcards: [f"bam/{wildcards.sample}_MT.{i}.sai" for i in range(0, 3)]
    output:
        temp("bam/{sample}_MT_0.bam")
    message: "Generate BAM for merged reads of sample {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.2.2"
    params: reffa = f"{workflow.basedir}/resources/NC_012920_1000.fa",
            readgroup = lambda wildcards: r'@RG\tID:{sample}\tSM:{sample}\tPL:illumina'.format(sample = wildcards.sample)
    threads: 2
    shell:
        """
        bwa samse \
            -r '{params.readgroup}' \
            -f /dev/stdout \
            {params.reffa} \
            {input.sais[0]} \
            {input.bam} | \
        samtools view -hb - > {output}
        """

rule bam_merge_wrap_sort:
    # This command merges the paired reads and single reads aligned against the
    # MT genome, wraps the alignment to its original length of 16 569 bp,
    # filters for a minimal read length of 30bp and a minimal mapping quality
    # of MQ25 and sorts the output
    input:
        bam_12 = "bam/{sample}_MT_12.bam",
        bam_0 = "bam/{sample}_MT_0.bam"
    output:
        temp("bam/{sample}_MTonly.sorted.bam")
    message: "Merge, wrap and sort the aligned reads of {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.2.2"
    params: merged_bam = "bam/{sample}_MT_120.bam",
            header = "bam/{sample}_MT_120.header",
            reffa = f"{workflow.basedir}/resources/NC_012920.fa",
            realigned_bam = "bam/{sample}_MT_120_realigned.bam"
    threads: 2
    shell:
        """
        if [ -f {params.merged_bam} ]; then
            rm {params.merged_bam}
        fi
        samtools merge -c \
            {params.merged_bam} \
            {input.bam_12} \
            {input.bam_0} 
        java -Xms512M -Xmx1G \
            -jar ${{CONDA_DEFAULT_ENV}}/share/circularmapper-1.93.4-1/realign-1.93.4.jar \
            -e 1000 \
            -r {params.reffa} \
            -f true \
            -x false \
            -i {params.merged_bam}
        samtools view -H {params.merged_bam} > {params.header}
        cat {params.header} <( \
            samtools view -L <(echo -e "MT\t0\t16569") -q 25 {params.realigned_bam} | \
            bioawk -c sam '{{if (length($seq) >= 30){{print}}}}'
        ) | \
        samtools view -Sb - | \
        samtools sort -o {output} -
        rm {params.merged_bam} {params.realigned_bam} {params.header}
        """


rule bam_rmdup:
    input:
        "bam/{sample}_MTonly.sorted.bam"
    output:
        "bam/{sample}_MTonly.sorted.rmdup.bam"
    message: "Remove duplicate reads from {wildcards.sample}"
    log: "logs/dedup/{sample}_dedup.log"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.2.2"
    params: reffa = f"{workflow.basedir}/resources/NC_012920.fa"
    shell:
        """
        java -Xms512M -Xmx1G \
            -jar ${{CONDA_DEFAULT_ENV}}/share/dedup-0.12.3-1/DeDup-0.12.3.jar \
            -i {input} \
            --merged \
            -o bam
        mv bam/{wildcards.sample}_MTonly.sorted.hist logs/dedup/{wildcards.sample}_dedup.hist
        mv bam/{wildcards.sample}_MTonly.sorted.log {log}
        samtools sort bam/{wildcards.sample}_MTonly.sorted_rmdup.bam | \
        samtools calmd -b - {params.reffa} > {output}
        samtools index {output}
        rm bam/{wildcards.sample}_MTonly.sorted_rmdup.bam
        """

# Quality checks

rule damage_profiler:
    input:
        "bam/{sample}_MTonly.sorted.rmdup.bam"
    output:
        "qual/{sample}/identity_histogram.pdf"
    message: "Generate damage profile for {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.2.2"
    params: dir = "qual",
            tmpdir = "qual/{sample}_MTonly.sorted.rmdup",
            outdir = "qual/{sample}",
            reffa = f"{workflow.basedir}/resources/NC_012920.fa"
    shell:
        """
        java -Xms512M -Xmx1G \
            -Djava.awt.headless=true \
            -jar ${{CONDA_DEFAULT_ENV}}/share/damageprofiler-0.4.4-1/DamageProfiler-0.4.4.jar \
            -i {input} \
            -o {params.dir} \
            -r {params.reffa}
        mv {params.tmpdir}/* {params.outdir}/
        rmdir {params.tmpdir}/
        """


rule flagstat:
    input:
        "bam/{sample}_MTonly.sorted.rmdup.bam"
    output:
        "logs/flagstat/{sample}.flagstat"
    message: "Determine the number of aligned reads for {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.2.2"
    shell:
        "samtools flagstat {input} > {output}"


rule mixemt:
    input:
        bam = "bam/{sample}_MTonly.sorted.rmdup.bam",
        flagstat = "logs/flagstat/{sample}.flagstat"
    output:
        "logs/mixemt/{sample}.mixemt.pos.tab"
    message: "Determine mtDNA contamination for {wildcards.sample} using mixEMT"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.2.2"
    priority: -20
    params: 
            subsampling = lambda wildcards: True if float(mixemt_downsampling(f"logs/flagstat/{wildcards.sample}.flagstat")) < 1 else False,
            subsampling_fraction = lambda wildcards: mixemt_downsampling(f"logs/flagstat/{wildcards.sample}.flagstat"),
            subbam = "bam/{sample}_MTonly.sorted.rmdup.downsampled.bam",
            mixemtprefix = "logs/mixemt/{sample}.mixemt"
    shell:
        """
        if [[ -f {params.subbam} ]]; then
            rm {params.subbam}
        fi
        if [[ {params.subsampling} = "True" ]]; then  # subsampling
            samtools view -bh \
                          -s {params.subsampling_fraction} \
                          -o {params.subbam} \
                          {input.bam}
        else
            ln -s ${{PWD}}/{input.bam} {params.subbam}
        fi
        samtools index {params.subbam}
        python ${{CONDA_DEFAULT_ENV}}/bin/mixemt -v -t \
                {params.mixemtprefix} \
                {params.subbam} \
                > {params.mixemtprefix}.log \
                2> {params.mixemtprefix}.stderr || \
        if [[ ($(wc -l < {params.mixemtprefix}.log) = "0") && ($(tail -1 {params.mixemtprefix}.stderr) = "0 contributors passed filtering steps.") ]]; then
            echo -e "hap1\tNA\tNA\t0" > {params.mixemtprefix}.log
            touch {params.mixemtprefix}.pos.tab
        fi
        rm {params.subbam}
        """


rule seqdepth:
    input:
        expand("bam/{sample}_MTonly.sorted.rmdup.bam", sample=SAMPLES)
    output:
        "logs/seqdepth.csv"
    message: "Run samtools depth to determine the coverage across the MT genome"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.2.2"
    params: 
        header = "chr\tpos\t{}".format("\t".join([sm for sm in SAMPLES])),
        bams = " ".join(["bam/{}_MTonly.sorted.rmdup.bam".format(sm) for sm in SAMPLES])
    shell:
        """
        echo "{params.header}" > {output}
        samtools depth -a -r MT:1-16569 {params.bams} >> {output}
        """

rule haplogrep2:
    input:
        "snpAD/{sample}.snpAD.fasta"
    output:
        "logs/haplogrep2/{sample}.hsd"
    message: "Determine mtDNA haplogroup for {wildcards.sample} using HaploGrep2"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.2.2"
    shell:
        """
        java -Xms256m -Xmx1G -jar {workflow.basedir}/resources/haplogrep-2.1.19.jar \
            --in {input} \
            --extend-report \
            --format fasta \
            --out {output}
        """

# Contamination estimate using contamMix

rule contamMix_create_sequencePanel:
    input:
        "snpAD/{sample}.snpAD.fasta"
    output:
        "logs/contamMix/{sample}/sequence_panel.fasta"
    message: "Align consensus sequence of {wildcards.sample} to panel of 311 modern humans"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.2.2"
    params: 
        panel = f"{workflow.basedir}/resources/311hu+rCRS.fas"
    shell:
        """
        cat {input} {params.panel} | \
        mafft - > {output}
        """

rule contamMix_align_against_consensus:
    input:
        flagstat = "logs/flagstat/{sample}.flagstat",
        fas = "snpAD/{sample}.snpAD.fasta",
        bam = "bam/{sample}_MTonly.sorted.rmdup.bam"
    output:
        "logs/contamMix/{sample}/{sample}.consensus_aligned.bam"
    message: "Align sequences of {wildcards.sample} against its consensus sequence"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.2.2"
    params: 
        readgroup = lambda wildcards: r'@RG\tID:{sample}\tSM:{sample}\tPL:illumina'.format(sample = wildcards.sample),
        subsampling = lambda wildcards: True if float(mixemt_downsampling(f"logs/flagstat/{wildcards.sample}.flagstat")) < 1 else False,
        subsampling_fraction = lambda wildcards: mixemt_downsampling(f"logs/flagstat/{wildcards.sample}.flagstat"),
        subbam = "logs/contamMix/{sample}/{sample}.subsampled.bam"
    threads: 4
    shell:
        """
        # Create FastA of consensus sequence with overhang and index with BWA
        bioawk -c fastx '{{
            print ">" $name; \
            print $seq substr($seq,1,1000); 
        }}' {input.fas} > logs/contamMix/{wildcards.sample}/{wildcards.sample}.fasta
        bwa index logs/contamMix/{wildcards.sample}/{wildcards.sample}.fasta
        # Sub-sample BAM file to 30 000 reads
        if [[ {params.subsampling} = "True" ]]; then  # subsampling
            samtools view -bh \
                          -s {params.subsampling_fraction} \
                          -o {params.subbam} \
                          {input.bam}
        else
            ln -s ${{PWD}}/{input.bam} {params.subbam}
        fi
        # BWA aln against sample's consensus sequence
        bwa aln \
            -t {threads} \
            -n 0.01 -o 2 -l 16500 \
            -b -1 \
            -f logs/contamMix/{wildcards.sample}/{wildcards.sample}.1.sai \
            logs/contamMix/{wildcards.sample}/{wildcards.sample}.fasta \
            {params.subbam}
        bwa aln \
            -t {threads} \
            -n 0.01 -o 2 -l 16500 \
            -b -2 \
            -f logs/contamMix/{wildcards.sample}/{wildcards.sample}.2.sai \
            logs/contamMix/{wildcards.sample}/{wildcards.sample}.fasta \
            {params.subbam}
        bwa aln \
            -t {threads} \
            -n 0.01 -o 2 -l 16500 \
            -b -0 \
            -f logs/contamMix/{wildcards.sample}/{wildcards.sample}.0.sai \
            logs/contamMix/{wildcards.sample}/{wildcards.sample}.fasta \
            {params.subbam}
        # BWA sampe and samse to create alignment files
        bwa sampe \
            -r '{params.readgroup}' \
            -f /dev/stdout \
            logs/contamMix/{wildcards.sample}/{wildcards.sample}.fasta \
            logs/contamMix/{wildcards.sample}/{wildcards.sample}.1.sai \
            logs/contamMix/{wildcards.sample}/{wildcards.sample}.2.sai \
            {params.subbam} \
            {params.subbam} | \
        samtools view -hb - > logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT_12.bam
        bwa samse \
            -r '{params.readgroup}' \
            -f /dev/stdout \
            logs/contamMix/{wildcards.sample}/{wildcards.sample}.fasta \
            logs/contamMix/{wildcards.sample}/{wildcards.sample}.0.sai \
            {params.subbam} | \
        samtools view -hb - > logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT_0.bam
        # Merge sam files, re-wrap and sort
        if [ -f logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT_120.bam ]; then
            rm logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT_120.bam
        fi
        samtools merge -c \
            logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT_120.bam \
            logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT_12.bam \
            logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT_0.bam 
        echo "{wildcards.sample}" > logs/contamMix/{wildcards.sample}/{wildcards.sample}.fasta_1000_elongated
        java -Xms512M -Xmx1G \
            -jar ${{CONDA_DEFAULT_ENV}}/share/circularmapper-1.93.4-1/realign-1.93.4.jar \
            -e 1000 \
            -r logs/contamMix/{wildcards.sample}/{wildcards.sample}.fasta \
            -f true \
            -x false \
            -i logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT_120.bam
        samtools view -H logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT_120_realigned.bam \
                > logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT_120.header
        cat logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT_120.header <( \
            samtools view -L <(echo -e "{wildcards.sample}\t0\t16569") -q 25 logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT_120_realigned.bam | \
            bioawk -c sam '{{if (length($seq) >= 30){{print}}}}'
        ) | \
        samtools view -Sb - | \
        samtools sort -o {output} -
        # Clean
        rm logs/contamMix/{wildcards.sample}/{wildcards.sample}.*.sai \
           logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT* \
           logs/contamMix/{wildcards.sample}/{wildcards.sample}.fasta* \
           {params.subbam}
        """

rule contamMix_estimate:
    input:
        aln = "logs/contamMix/{sample}/sequence_panel.fasta",
        bam = "logs/contamMix/{sample}/{sample}.consensus_aligned.bam"
    output:
        "logs/contamMix/{sample}/contamMix_log.txt"
    message: "Use contamMix to estimate contamination of sample {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.2.2"
    priority: -10
    params: 
        tarball = f"{workflow.basedir}/resources/contamMix_1.0-10.tar.gz",
        contamMix_exec = f"{workflow.basedir}/resources/contamMix/exec/estimate.R"
    threads: 3
    shell:
        """
        if [[ ! -f ${{CONDA_DEFAULT_ENV}}/lib/R/etc/ldpaths ]]; then
            cp /etc/R/ldpaths ${{CONDA_DEFAULT_ENV}}/lib/R/etc/ldpaths
        fi
        R -e 'if (!("contamMix" %in% installed.packages())) {{ \
                install.packages(c("coda", "getopt"), , repos="http://cran.us.r-project.org"); \
                install.packages("{params.tarball}", type="source")}}'
        {params.contamMix_exec} \
                --samFn {input.bam} \
                --malnFn {input.aln} \
                --tabOutput TRUE > {output}
        """

# Genotype calling mit snpAD

rule bam2snpAD:
    input: 
        "bam/{sample}_MTonly.sorted.rmdup.bam"
    output:
        temp("snpAD/{sample}.snpad_input")
    message: "Convert BAM file into snpAD input format for sample {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.2.2"
    params:
        bam2snpAD = f"{workflow.basedir}/resources/snpAD-0.3.3/Bam2snpAD",
        reffasta = f"{workflow.basedir}/resources/NC_012920.fa"   
    shell:
        """
        {params.bam2snpAD} \
                -Q 25 \
                -q 30 \
                -r MT \
                -f {params.reffasta} \
                {input} > {output}
        """

rule snpAD_estimation:
    input:
        "snpAD/{sample}.snpad_input"
    output:
        priors = temp("snpAD/{sample}.priors.txt"),
        errors = temp("snpAD/{sample}.errors.txt")
    message: "Estimate the genotype likelihoods using snpAD for sample {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.2.2"
    params:
        snpAD = f"{workflow.basedir}/resources/snpAD-0.3.3/snpAD",
    log: "snpAD/{sample}.snpAD.log"
    shell:
        """
        {params.snpAD} \
                -o {output.priors} \
                -O {output.errors} \
                {input} > {log} 2>&1
        """

rule snpAD_call:
    input:
        snpAD = "snpAD/{sample}.snpad_input",
        priors = "snpAD/{sample}.priors.txt",
        errors = "snpAD/{sample}.errors.txt"
    output:
        "snpAD/{sample}.snpAD.vcf"
    message: "Call the genotypes using snpAD for sample {wildcards.sample} fixing the likelihood of a heterozygous genotype to a very small number"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.2.2"
    params:
        snpADcall = f"{workflow.basedir}/resources/snpAD-0.3.3/snpADCall",
    shell:
        """
        {params.snpADcall} \
                -e <(bioawk -t '$4 < 1e6{{print $1, $2, $3, $4}}' {input.errors}) \
                -p $(cat {input.priors} |perl -F',' -lane 'for ($i=0 ; $i<4 ; $i++) {{ $F[$i]=$F[$i]/($F[0]+$F[1]+$F[2]+$F[3]); }} for ($i=4; $i<@F; $i++) {{ $F[$i] = "1e-320" ; }} print join ",",@F;') \
                {input.snpAD} > {output}
        """

rule snpAD_vcf2fa:
    input:
        "snpAD/{sample}.snpAD.vcf"
    output:
        fa = "snpAD/{sample}.snpAD.fasta",
        tbi = "snpAD/{sample}.snpAD.vcf.gz.tbi"
    message: "Convert snpAD VCF file into a FastA file for sample {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.2.2"
    params:
        vcf2fa = f"{workflow.basedir}/resources/snpAD-0.3.3/vcf2fasta.pl",
    shell:
        """
        {params.vcf2fa} < {input} | \
        sed 's/>MT/>{wildcards.sample}/' > {output.fa}
        bgzip -f {input}
        tabix -f {input}.gz
        """

# Summary

rule summary:
    input: expand("bam/{sample}_MTonly.sorted.rmdup.bam", sample=SAMPLES),
           expand("qual/{sample}/identity_histogram.pdf", sample=SAMPLES),
           expand("logs/mixemt/{sample}.mixemt.pos.tab", sample=SAMPLES),
           expand("logs/contamMix/{sample}/contamMix_log.txt", sample=SAMPLES),
           "logs/seqdepth.csv",
           expand("logs/haplogrep2/{sample}.hsd", sample=SAMPLES),
           expand("snpAD/{sample}.snpAD.fasta", sample=SAMPLES)
    output: "{projdir}/summary_table.csv"
    message: "Summarise the results in a table"
    version: "0.2.2"
    run:
        R("""
          # Install packages if not available
          packages <- c("ape", "data.table", "tidyverse")
          missing_packages <- packages[!(packages %in% installed.packages())]
          if (length(missing_packages) > 0) {{
              install.packages(missing_packages,
                               repos = "http://cran.us.r-project.org")
          }}
          
          # Load libaries
          library(ape)
          library(data.table)
          library(tidyverse)
          
          # Start summary
          ## number of reads
          flagstat_fns <- list.files("logs/flagstat",
                                     pattern = "\\\.flagstat",
                                     full.names = T)
          nreads_MT <- map_df(flagstat_fns, function(fn) {{
                              tibble(sample = str_replace(basename(fn), "\\\.flagstat", ""),
                                     nReads = as.numeric(unlist(str_split(readLines(fn)[1], " "))[1]))
                              }})
          ## coverage
          seqdepth <- fread("logs/seqdepth.csv") %>%
                      gather("sample", "cov", 3:ncol(.)) %>%
                      group_by(sample) %>%
                      summarise(meanCov = mean(cov),
                                medianCov = median(cov),
                                minCov = min(cov),
                                maxCov = max(cov),
                                sdCov = sd(cov),
                                cov_gt_5x = sum(cov >= 5))
          
          ## read length
          mode_of_readlength <- function(rl) {{
              rl %>%
              filter(n == max(n)) %>%
              arrange(desc(sumRL)) %>%
              pull(Length) %>%
              .[1]
          }}
          readlength <- map_df(seqdepth$sample, function(s) {{
                               readlength_dist <- fread(paste("qual", s, "lgdistribution.txt", sep="/")) %>%
                                                  group_by(Length) %>%
                                                  summarise(n = sum(Occurrences)) %>%
                                                  mutate(sumRL = Length * n)
                               tibble(sample = s,
                                      meanReadLength = sum(readlength_dist$sumRL) / sum(readlength_dist$n),
                                      modeReadLength = mode_of_readlength(readlength_dist))
                               }})
          
          ## contamination
          ### mixEMT
          mixemt_fns <- list.files("logs/mixemt",
                                   pattern = "\\\.mixemt\\\.log", full.names = T)
          mixemt <- map_df(mixemt_fns, function(fn) {{
                           fread(fn, header = F,
                                 col.names = c("comp", "hg", "contr", "noReads")) %>%
                           mutate(sample = str_replace(basename(fn), "\\\.mixemt\\\.log", ""),
                                  contr = paste0(contr * 100, "%"),
                                  label = paste0(hg, " (", contr, ",", noReads, ")")) %>%
                           group_by(sample) %>%
                           summarise(mixEMT = str_replace(toString(label), ", ", ";"))
                           }})
          ### contamMix
          contamMix <- map_df(seqdepth$sample, function(s) {{
                              fread(paste("logs/contamMix", s, "contamMix_log.txt", sep="/"),
                                    select = 1:4,
                                    col.names = c("error", "propAuth",
                                                  "propAuth_lowerQuantile",
                                                  "propAuth_upperQuantile")) %>%
                              mutate(sample = s,
                                     label = paste0(round((1 - propAuth) * 100, 2), "% (",
                                                    round((1 - propAuth_upperQuantile) * 100, 2), "-",
                                                    round((1 - propAuth_lowerQuantile) * 100, 2), "%)"),
                                     flagContamination = 1 - propAuth >= 0.1) %>%
                              select(sample, proportionContamination = label, flagContamination)
                              }})
          
          ## Consensus sequence
          fas_snpAD_fns <- list.files("snpAD",
                                      pattern = "\\\.snpAD\\\.fasta", full.names = T)
          fas_snpAD <- lapply(fas_snpAD_fns, function(fn) read.dna(fn, "fasta", as.character = T))
          ### Compare consensus sequences
          consensus_comparison <- map_df(seq(1, length(fas_snpAD)), function(i) {{
                                          tibble(sample = str_replace(basename(fas_snpAD_fns[i]),
                                                                      "\\\.snpAD\\\.fasta", ""),
                                                 noNs = sum(fas_snpAD[[i]] == "n"),
                                                 percentNs = round(noNs * 100 / 16569, 1),
                                                 flagMissingness = percentNs > 1)
                                        }})
          
          ## Haplogroup assignment
          count_no_polys <- function(s) {{
            # Count the number of polymorphisms that are not caused by missing data
            if (is.na(s)) {{
              rawstring <- unlist(str_split(s, " "))
              fltstring <- rawstring[!str_detect(rawstring, "N$")]
              return(length(fltstring))
            }} else return(0)
          }}

          hsd_fns <- list.files("logs/haplogrep2",
                                pattern = "\\\.hsd", full.names = T)
          haplogrep <- map_df(hsd_fns, function(fn) {{
                           fread(fn,
                                 colClasses = list("character"=c("SampleID"))) %>%
                           mutate(noNotFoundPolys = count_no_polys(Not_Found_Polys),
                                  noRemainingPolys = count_no_polys(Remaining_Polys)) %>%
                           select(sample = SampleID,
                                  haplogroup = Haplogroup,
                                  hgQ = Quality,
                                  noNotFoundPolys, noRemainingPolys)
                              }})
                    
          # Merge summary 
          summary_table <- nreads_MT %>%
                           rename(`number of reads` = nReads) %>%
                           left_join(seqdepth, by = "sample") %>%
                           select(-medianCov) %>%
                           mutate(meanCov = round(meanCov, 1),
                                  sdCov = round(sdCov, 0)) %>%
                           rename(`mean coverage` = meanCov,
                                  `std coverage` = sdCov,
                                  `sites >= 5-fold coverage` = cov_gt_5x,
                                  `minimum coverage` = minCov,
                                  `maximum coverage` = maxCov) %>%
                           left_join(readlength, by = "sample") %>%
                           select(-meanReadLength) %>%
                           rename(`mode of read length` = modeReadLength) %>%
                           left_join(consensus_comparison, by = "sample") %>%
                           select(-noNs) %>%
                           rename(`% of Ns in snpAD consensus sequence` = percentNs,
                                  `flagged due to number of missing bases (> 1%)` = flagMissingness) %>%
                           left_join(haplogrep, by = "sample") %>%
                           rename(`haplogroup` = haplogroup,
                                  `haplogroup quality` = hgQ,
                                  `no. of not-found polys` = noNotFoundPolys,
                                  `no. of remaining polys` = noRemainingPolys) %>%
                           left_join(contamMix, by = "sample") %>%
                           rename(`proportion of contamination DNA (contamMix)` = proportionContamination,
                                  `flagged due to contamination (> 10%)` = flagContamination) %>%
                           left_join(mixemt, by = "sample") %>%
                           rename(`haplogroup contributions (mixEMT)` = mixEMT) %>%
                           mutate(quality = ifelse(`flagged due to number of missing bases (> 1%)` | `flagged due to contamination (> 10%)`,
                                                   "low", "high"))
          fwrite(summary_table,
                 sep = "\t",
                 file = "{output}")
        """)

rule copy_tmp_to_proj:
    input: expand("{projdir}/summary_table.csv", projdir=[PROJDIR])
    output: expand("{projdir}/fasta/{sample}.fa", projdir=[PROJDIR], sample=SAMPLES)
    message: "Copy results from tmp to project folder"
    version: "0.2.2"
    shell:
        """
        mkdir -p {PROJDIR}/sample_stats
        for sm in {SAMPLES}; do
            echo ${{sm}}
            # Create per-sample folder
            mkdir -p {PROJDIR}/sample_stats/${{sm}}
            # Copy non-clipped BAM file
            cp -r bam/${{sm}}_MTonly.sorted.rmdup.bam* {PROJDIR}/sample_stats/${{sm}}/
            # DeDup logs
            cp -r logs/dedup/${{sm}}_dedup.* {PROJDIR}/sample_stats/${{sm}}/
            # HaploGrep
            cp logs/haplogrep2/${{sm}}.hsd {PROJDIR}/sample_stats/${{sm}}/
            # MixEMT
            cp logs/mixemt/${{sm}}.* {PROJDIR}/sample_stats/${{sm}}/
            # DamageProfiler
            cp qual/${{sm}}/3pGtoA_freq.txt {PROJDIR}/sample_stats/${{sm}}/${{sm}}_3pGtoA_freq.txt
            cp qual/${{sm}}/5pCtoT_freq.txt {PROJDIR}/sample_stats/${{sm}}/${{sm}}_5pCtoT_freq.txt
            cp qual/${{sm}}/DamagePlot.pdf {PROJDIR}/sample_stats/${{sm}}/${{sm}}_DamagePlot.pdf
            cp qual/${{sm}}/Length_plot.pdf {PROJDIR}/sample_stats/${{sm}}/${{sm}}_Length_plot.pdf
            cp qual/${{sm}}/lgdistribution.txt {PROJDIR}/sample_stats/${{sm}}/${{sm}}_lgdistribution.txt
            cp qual/${{sm}}/misincorporation.txt {PROJDIR}/sample_stats/${{sm}}/${{sm}}_misincorporation.txt
            # snpAD
            cp snpAD/${{sm}}.snpAD.vcf.gz* {PROJDIR}/sample_stats/${{sm}}/
            cp snpAD/${{sm}}.snpAD.fasta {PROJDIR}/fasta/${{sm}}.fa
        done
        cp logs/seqdepth.csv {PROJDIR}/sample_stats/
        """

# Clean temporary output

rule clean_tmp:
    params:
        tmpdir = config['tmpdir']
    shell:
        "rm -r {params.tmpdir}"

