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

localrules: link_index, seqdepth, summary, copy_tmp_to_proj

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
        header = "chr\tpos\t{}".format("\t".join([sm for sm in SAMPLES])),
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

# Contamination estimate using contamMix

rule contamMix_create_sequencePanel:
    input:
        "fasta/{sample}_angsd.fa"
    output:
        "logs/contamMix/{sample}/sequence_panel.fasta"
    message: "Align consensus sequence of {wildcards.sample} to panel of 311 modern humans"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.1"
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
        fas = "fasta/{sample}_angsd.fa",
        bam = "bam/{sample}_MTonly.sorted.rmdup.bam"
    output:
        "logs/contamMix/{sample}/{sample}.consensus_aligned.bam"
    message: "Align sequences of {wildcards.sample} against its consensus sequence"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.1"
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
        samtools sort -o {output} logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT_120_realigned.bam
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
    version: "0.1"
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

rule summary:
    input: expand("bam/{sample}_MTonly.sorted.rmdup.bam", sample=SAMPLES),
           expand("qual/{sample}/identity_histogram.pdf", sample=SAMPLES),
           expand("logs/mixemt/{sample}.mixemt.pos.tab", sample=SAMPLES),
           expand("logs/contamMix/{sample}/contamMix_log.txt", sample=SAMPLES),
           "logs/seqdepth.csv",
           expand("fasta/{sample}_angsd.fa", sample=SAMPLES),
           expand("fasta/{sample}_angsd_unclipped.fa", sample=SAMPLES),
           expand("logs/haplogrep2/{sample}.hsd", sample=SAMPLES)
    output: "{projdir}/summary_table.csv"
    message: "Summarise the results in a table"
    version: "0.1"
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
                                sdCov = sd(cov),
                                cov_gt_5x = sum(cov >= 5))
          
          ## read length
          mode_of_readlength <- function(rl) {{
              rl %>%
              filter(n == max(n)) %>%
              pull(Length)
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
                                     label = paste0(round(propAuth * 100, 2), "% (",
                                                    round(propAuth_lowerQuantile * 100, 2), "-",
                                                    round(propAuth_upperQuantile * 100, 2), "%)")) %>%
                              select(sample, proportionAuthentic = label)
                              }})
          
          ## Number of trimmed bases
          ### 5' end
          fiveP_CtoT <- map_df(seqdepth$sample, function(s) {{
                          damage <- fread(paste("qual", s, "5pCtoT_freq.txt", sep="/")) 
                          threshold <- mean(damage$`5pC>T`) + sd(damage$`5pC>T`)
                          cut_position <- filter(damage, `5pC>T` <= threshold) %>%
                                          slice(1) %>%
                                          pull(pos) - 1
                          tibble(sample = s,
                                 `5'_trBases` = ifelse(max(damage$`5pC>T`) >= 0.02,
                                                       cut_position, 0))
                        }})
          ### 3' end
          threeP_GtoA <- map_df(seqdepth$sample, function(s) {{
                          damage <- fread(paste("qual", s, "3pGtoA_freq.txt", sep="/")) 
                          threshold <- mean(damage$`3pG>A`) + sd(damage$`3pG>A`)
                          cut_position <- filter(damage, `3pG>A` <= threshold) %>%
                                          slice(1) %>%
                                          pull(pos) - 1
                          tibble(sample = s,
                                 `3'_trBases` = ifelse(max(damage$`3pG>A`) >= 0.02,
                                                       cut_position, 0))
                        }})
          
          ## Consensus sequence comparison
          ### Clipped BAM files
          fas_clipped_fns <- list.files("fasta",
                                        pattern = "_angsd\\\.fa", full.names = T)
          fas_clipped <- lapply(fas_clipped_fns, function(fn) read.dna(fn, "fasta", as.character = T))
          
          ### Non-clipped BAM files
          fas_nonclipped_fns <- list.files("fasta",
                                        pattern = "_angsd_unclipped\\\.fa", full.names = T)
          fas_nonclipped <- lapply(fas_nonclipped_fns, function(fn) read.dna(fn, "fasta", as.character = T))
          ### Compare consensus sequences
          consensus_comparison <- map_df(seq(1, length(fas_clipped)), function(i) {{
                                          tibble(sample = str_replace(basename(fas_clipped_fns[i]),
                                                                      "_angsd\\\.fa", ""),
                                                 noNs = sum(fas_clipped[[i]] == "n"),
                                                 percentNs = round(noNs * 100 / 16569, 1),
                                                 noDiffs = sum(fas_clipped[[i]] != fas_nonclipped[[i]]))
                                        }})
          
          ## Haplogroup assignment
          hsd_fns <- list.files("logs/haplogrep2",
                                pattern = "\\\.hsd", full.names = T)
          haplogrep <- map_df(hsd_fns, function(fn) {{
                              fread(fn) %>%
                              mutate(noNotFoundPolys = ifelse(is.na(Not_Found_Polys), "-",
                                                              as.character(str_count(Not_Found_Polys, " ") + 1)),
                                     noRemainingPolys = ifelse(is.na(Remaining_Polys), "-",
                                                              as.character(str_count(Remaining_Polys, " ") + 1))) %>%
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
                                  `sites >= 5-fold coverage` = cov_gt_5x) %>%
                           left_join(readlength, by = "sample") %>%
                           select(-meanReadLength) %>%
                           rename(`mode of read length` = modeReadLength) %>%
                           left_join(fiveP_CtoT, by = "sample") %>%
                           rename(`no. of trimmed bases at 5'-end` = `5'_trBases`) %>%
                           left_join(threeP_GtoA, by = "sample") %>%
                           rename(`no. of trimmed bases at 3'-end` = `3'_trBases`) %>%
                           left_join(consensus_comparison, by = "sample") %>%
                           select(-noNs) %>%
                           rename(`% of Ns in consensus sequence` = percentNs,
                                  `no. of differences clipped vs. non-clipped` = noDiffs) %>%
                           left_join(haplogrep, by = "sample") %>%
                           rename(`haplogroup clipped reads` = haplogroup,
                                  `haplogroup quality` = hgQ,
                                  `no. of not-found polys` = noNotFoundPolys,
                                  `no. of remaining polys` = noRemainingPolys) %>%
                           left_join(contamMix, by = "sample") %>%
                           rename(`proportion of authentic DNA (contamMix)` = proportionAuthentic) %>%
                           left_join(mixemt, by = "sample") %>%
                           rename(`haplogroup contributions (mixEMT)` = mixEMT)   
          fwrite(summary_table,
                 sep = "\t",
                 file = "{output}")
        """)

rule copy_tmp_to_proj:
    output: expand("{projdir}/fasta/{sample}.fa", projdir=[PROJDIR], sample=SAMPLES)
    message: "Copy results from tmp to project folder"
    version: "0.1"
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
            # Unclipped FastA
            cp fasta/${{sm}}_angsd_unclipped.fa {PROJDIR}/sample_stats/${{sm}}/${{sm}}_unclipped_consensus.fa
            # Clipped FastA
            cp fasta/${{sm}}_angsd.fa {PROJDIR}/fasta/${{sm}}.fa
        done
        cp logs/seqdepth.csv {PROJDIR}/sample_stats/
        """

# Clean temporary output

rule clean_tmp:
    params:
        tmpdir = config['tmpdir']
    shell:
        "rm -r {params.tmpdir}"

