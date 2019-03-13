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
        return "{:.4f}".format(40000 / nreads)
    else:
        return 1.0


def calcuate_trim_threshold(dirname, end):
    ''' Infer the number of bases to trim on both the 5p and 3p end of reads to
        avoid having issues with ancient DNA damage pattern for subsequent
        consensus calling.
        The threshold is determined and by calculating the mean and standard
        deviation of the C > T and G > A substitution frequency for the last 25
        bases as reported by DamageProfiler, respectively. As we expect a lower
        substitution frequency towards the middle of a read, we determine the
        first position, for which the subsitution frequency is within one
        standard deviation of the mean, and remove all bases to the left of it
        on the 5p end and all bases to the right on the 3p end.
    '''
    if os.path.isdir(dirname):
        if end == "5p":  # 5' end
            damageprof = pd.read_csv(dirname + "/5pCtoT_freq.txt", sep="\t",
                                comment="#")
        elif end == "3p":   # 3' end
            damageprof = pd.read_csv(dirname + "/3pGtoA_freq.txt", sep="\t",
                                comment="#")
        damageprof.columns = ['pos', 'freq']

        if damageprof['freq'].max() >= 0.02:  # only cut with C to T frequency >= 2%
            threshold = damageprof['freq'].mean() + damageprof['freq'].std()
            cut_position = damageprof.loc[damageprof['freq'] <= threshold]. \
                    iloc[0].name
        else:
            cut_position = 0
        return cut_position
    else:
        cut_position = 0


# Snakemake rules
wildcard_constraints:
    sample = config['sampleIDconstraint']

localrules: link_index

rule all:
    input: expand("bam/{sample}_MTonly.sorted.rmdup.bam", sample=SAMPLES),
           expand("qual/{sample}/identity_histogram.pdf", sample=SAMPLES),
           expand("logs/mixemt/{sample}.mixemt.pos.tab", sample=SAMPLES),
           "logs/seqdepth.csv",
           expand("fasta/{sample}_angsd.fa", sample=SAMPLES),
           expand("fasta/{sample}_angsd_unclipped.fa", sample=SAMPLES),
           expand("logs/haplogrep2/{sample}.hsd", sample=SAMPLES)


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
    threads: 2
    shell:
        """
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
        java -Xms512M -Xmx1G \
            -jar ${{CONDA_DEFAULT_ENV}}/share/dedup-0.12.3-1/DeDup-0.12.3.jar \
            -i {input} \
            --merged \
            -o bam
        mv bam/{wildcards.sample}_MTonly.sorted.hist logs/dedup/{wildcards.sample}_dedup.hist
        mv bam/{wildcards.sample}_MTonly.sorted.log {log}
        samtools sort -o {output} bam/{wildcards.sample}_MTonly.sorted_rmdup.bam
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
    version: "0.1"
    params: dir = "qual",
            tmpdir = "qual/{sample}_MTonly.sorted.rmdup",
            outdir = "qual/{sample}",
            reffa = f"{workflow.basedir}/resources/NC_012920.fa"
    shell:
        """
        java -Xms512M -Xmx1G \
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
    version: "0.1"
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
    version: "0.1"
    params: 
            subsampling = lambda wildcards: True if float(mixemt_downsampling(f"logs/flagstat/{wildcards.sample}.flagstat")) < 1 else False,
            subsampling_fraction = lambda wildcards: mixemt_downsampling(f"logs/flagstat/{wildcards.sample}.flagstat"),
            subbam = "bam/{sample}_MTonly.sorted.rmdup.downsampled.bam",
            mixemtprefix = "logs/mixemt/{sample}.mixemt"
    shell:
        """
        if [[ {params.subsampling} = "True" ]]; then  # subsampling
            samtools view -bh \
                          -s {params.subsampling_fraction} \
                          -o {params.subbam} \
                          {input.bam}
        else
            ln -s ${{PWD}}/{input.bam} {params.subbam}
        fi
        samtools index {params.subbam}
        mixemt -v -t \
                {params.mixemtprefix} \
                {params.subbam} \
                > {params.mixemtprefix}.log \
                2> {params.mixemtprefix}.stderr
        rm {params.subbam}*
        """


rule seqdepth:
    input:
        expand("bam/{sample}_MTonly.sorted.rmdup.bam", sample=SAMPLES)
    output:
        "logs/seqdepth.csv"
    message: "Run samtools depth to determine the coverage across the MT genome"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.1"
    params: 
        header = "chr\tpos\t{}\n".format("\t".join([sm for sm in SAMPLES])),
        bams = " ".join(["bam/{}_MTonly.sorted.rmdup.bam".format(sm) for sm in SAMPLES])
    shell:
        """
        echo "{params.header}" > {output}
        samtools depth -a -r MT {params.bams} >> {output}
        """

# Consensus calling

rule trim_bam:
    input:
        bam = "bam/{sample}_MTonly.sorted.rmdup.bam",
        damageprof = "qual/{sample}/identity_histogram.pdf"
    output:
        "bam/{sample}_MTonly.sorted.rmdup.clip.bam"
    message: "Soft-clip ends of reads with high amounts of DNA damage: {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.1"
    params:
        p3_nbases = lambda wildcards: calcuate_trim_threshold(f"qual/{wildcards.sample}", "3p"),
        p5_nbases = lambda wildcards: calcuate_trim_threshold(f"qual/{wildcards.sample}", "5p")
    shell:
        """
        bam trimBam {input.bam} {output} -L {params.p5_nbases} -R {params.p3_nbases}
        """

rule angsd_consensus:
    # Call consensus sequence as described by Ehler et al. (2018): amtDB 
    input: "bam/{sample}_MTonly.sorted.rmdup.clip.bam"
    output: "fasta/{sample}_angsd.fa"
    message: "Call consensus sequence using ANGSD following Ehler et al. (2018): {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.1"
    threads: 2
    params: reffasta = f"{workflow.basedir}/resources/NC_012920.fa",
            dir = "fasta"
    shell:
        """
        angsd \
            -i {input} \
            -minMapQ 30 \
            -minQ 20 \
            -doFasta 2 \
            -doCounts 1 \
            -ref {params.reffasta} \
            -out {params.dir}/{wildcards.sample}_angsd.tmp
        bioawk \
            -c fastx '{{print ">{wildcards.sample}"; print $seq}}' \
            {params.dir}/{wildcards.sample}_angsd.tmp.fa.gz > {output}
        rm {params.dir}/{wildcards.sample}_angsd.tmp.*
        """

rule angsd_consensus_unclipped:
    # Call consensus sequence as described by Ehler et al. (2018): amtDB 
    input: "bam/{sample}_MTonly.sorted.rmdup.bam"
    output: "fasta/{sample}_angsd_unclipped.fa"
    message: "Call consensus sequence using ANGSD following Ehler et al. (2018) on unclipped BAM file: {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.1"
    threads: 2
    params: reffasta = f"{workflow.basedir}/resources/NC_012920.fa",
            dir = "fasta"
    shell:
        """
        angsd \
            -i {input} \
            -minMapQ 30 \
            -minQ 20 \
            -doFasta 2 \
            -doCounts 1 \
            -ref {params.reffasta} \
            -out {params.dir}/{wildcards.sample}_angsd_unclipped.tmp
        bioawk \
            -c fastx '{{print ">{wildcards.sample}"; print $seq}}' \
            {params.dir}/{wildcards.sample}_angsd_unclipped.tmp.fa.gz > {output}
        rm {params.dir}/{wildcards.sample}_angsd_unclipped.tmp.*
        """


rule haplogrep2:
    input:
        "fasta/{sample}_angsd.fa"
    output:
        "logs/haplogrep2/{sample}.hsd"
    message: "Determine mtDNA haplogroup for {wildcards.sample} using HaploGrep2"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.1"
    shell:
        """
        java -jar {workflow.basedir}/resources/haplogrep-2.1.19.jar \
            --in {input} \
            --extend-report \
            --format fasta \
            --out {output}
        """

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
