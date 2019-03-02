################################################################################
# Snakefile for the processing of ancient DNA human mtDNA capture data
#
# Alex Huebner, 21/02/19
################################################################################

from snakemake.utils import min_version

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

# Snakemake rules
wildcard_constraints:
    sample = config['sampleIDconstraint']

localrules: link_index

rule all:
    input: expand("bam/{sample}_MTonly.sorted.rmdup.bam", sample=SAMPLES),
           expand("qual/{sample}/identity_histogram.pdf", sample=SAMPLES),
           expand("logs/flagstat/{sample}.flagstat", sample=SAMPLES)


# Prepare sorted, duplicate removed BAM files aligned only against the MT genome

rule link_index:
    # Link the raw data to file and index the BAM file in order to retrieve only
    # the MT reads
    output:
        bam = "rawdata/{sample}.bam",
        bai = "rawdata/{sample}.bam.bai"
    message: "Link the BAM file to the temp directory and index it: {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.1"
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
    version: "0.1"
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
    version: "0.1"
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
    version: "0.1"
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
    version: "0.1"
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
    input:
        bam_12 = "bam/{sample}_MT_12.bam",
        bam_0 = "bam/{sample}_MT_0.bam"
    output:
        temp("bam/{sample}_MTonly.sorted.bam")
    message: "Merge, wrap and sort the aligned reads of {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.1"
    params: merged_bam = "bam/{sample}_MT_120.bam",
            reffa = f"{workflow.basedir}/resources/NC_012920.fa",
            realigned_bam = "bam/{sample}_MT_120_realigned.bam"
    #threads: 4
    shell:
        """
        samtools merge \
            {params.merged_bam} \
            {input.bam_12} \
            {input.bam_0} 
        realignsamfile \
			-e 1000 \
			-r {params.reffa} \
			-f true \
			-x false \
			-i {params.merged_bam}
        samtools sort -o {output} {params.realigned_bam}
		rm {params.merged_bam} {params.realigned_bam}
        """


rule bam_rmdup:
    input:
        "bam/{sample}_MTonly.sorted.bam"
    output:
        "bam/{sample}_MTonly.sorted.rmdup.bam"
    message: "Remove duplicate reads from {wildcards.sample}"
    log: "logs/dedup/{sample}_dedup.log"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.1"
    shell:
        """
        dedup -i {input} \
              --merged \
              -o bam
        mv bam/{wildcards.sample}_MTonly.sorted.hist logs/dedup/{wildcards.sample}_dedup.hist
        mv bam/{wildcards.sample}_MTonly.sorted.log {log}
        mv bam/{wildcards.sample}_MTonly.sorted_rmdup.bam {output}
        """

# Quality checks

rule damage_profiler:
    input:
        "bam/{sample}_MTonly.sorted.rmdup.bam"
    output:
        "qual/{sample}/identity_histogram.pdf"
    message: "Generate damage profile for {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.1"
    params: dir = "qual",
            tmpdir = "qual/{sample}_MTonly.sorted.rmdup",
            outdir = "qual/{sample}",
            reffa = f"{workflow.basedir}/resources/NC_012920.fa"
    shell:
        """
        damageprofiler \
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
    version: "0.1"
    shell:
        "samtools flagstat {input} > {output}"


#rule mixemt:
    #input:
        #"bam/{sample}_MTonly.sorted.rmdup.bam"
    #output:
        #"logs/haplogrep2/{sample}.hsd"
    #message: "Determine mtDNA contamination for {wildcards.sample} using mixEMT"
    #conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    #version: "0.1"
    #params: dir = "qual",
            #tmpdir = "qual/{sample}_MTonly.sorted.rmdup",
            #outdir = "qual/{sample}",
            #reffa = f"{workflow.basedir}/resources/NC_012920.fa"
    #shell:
        #"""
        #READS=$(head -1 {input} | cut -f1 -d' ')
        #if [[ ${{READS}} -gt {config[NREADS]} ]]; then
            #FRAC=$(perl -E "say {config[NREADS]} / ${{READS}}")
            #{config[SAMTOOLS]} view -bh \
                    #-s ${{FRAC}} \
                    #-o {wildcards.id}-sub.bam {wildcards.id}.bam
            #{config[SAMTOOLS]} index {wildcards.id}-sub.bam
            #python2 {config[MIXEMT]} -v \
                    #-t {wildcards.id}.mixemt \
                    #{wildcards.id}-sub.bam > {output.log} 2> {log}
        #else
            #python2 {config[MIXEMT]} -v \
                    #-t {wildcards.id}.mixemt \
                    #{wildcards.id}.realigned.bam > {output.log} 2> {log}
        #fi
        #"""

#rule haplogrep2:
    #input:
        #"bam/{sample}_MTonly.sorted.rmdup.bam"
    #output:
        #"logs/haplogrep2/{sample}.hsd"
    #message: "Determine mtDNA haplogroup for {wildcards.sample} using HaploGrep2"
    #conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    #version: "0.1"
    #params: dir = "qual",
            #tmpdir = "qual/{sample}_MTonly.sorted.rmdup",
            #outdir = "qual/{sample}",
            #reffa = f"{workflow.basedir}/resources/NC_012920.fa"
    #shell:
        #"""
        #damageprofiler \
            #-i {input} \
            #-o {params.dir} \
            #-r {params.reffa}
        #mv {params.tmpdir} {params.outdir}
        #"""


#rule template:
    #input:
    #output:
    #message: ""
    #params: 
    #shell:
        #"""
        #"""

#### Rules
# Starting point BAM file aligned against hg19 plus decoy
# 3. Quality filter
#   * coverage
#   * average read length (maybe in DamageProfiler)
#   * contamination estimates
#        * mixEMT
#        * contamMix
#   * damage profile using DamageProfiler
#   * haplogroup using HaploGrep2
# 4. Determine whether UDG treated?
# 5. Calculate trim cut-off based on damage profile
# 6. Trim the ends of BAM file
# 7. Consensus calling using ANGSD
# 8. Consensus calling using snpAD without trimming
# 9. Summary statistics
