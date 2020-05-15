################################################################################
# Snakefile for the processing of sequencing data for human mtDNA
################################################################################

__author__ = ('Alexander Hübner (alexhbnr@gmail.com)')
__version__ = '0.3'
__date__ = '08 May 2020'

from snakemake.utils import min_version

from glob import glob
import os
import re
import subprocess 

import pandas as pd
import pysam

min_version("5.0")
shell.executable("/bin/bash")

workdir: config['tmpdir']

# Create directory for cluster file if it doesn't exist

if not os.path.isdir(f"{config['tmpdir']}/cluster_logs"):
    os.makedirs(f"{config['tmpdir']}/cluster_logs")
if not os.path.isdir(f"{config['tmpdir']}/tmp"):
    os.makedirs(f"{config['tmpdir']}/tmp")

# Read the sample list with sequence type (single-end, paired-end)
# Expected two columns: sample\tseqdatatype
SAMPLESLIST = pd.read_csv(config['samplelist'], sep="\t", index_col=['sample'])
SAMPLES = SAMPLESLIST.index.tolist()

# Extract project directory
PROJDIR = config['projdir']

# Auxilliary functions

def check_state(fn):
    """Check whether sample passed the minimum number of reads."""
    if os.path.isfile(fn):
        return next(open(fn, "rt")).rstrip()
    else:
        return "Pass"

def mixemt_downsampling(flagstatfn, targetreads=config['targetreads']):
    ''' Determines the number of reads in a BAM file based on samtools flagstat
        and calculates the fraction of reads necessary to obtain 40,000 reads,
        the suggested input into mixEMT. The fraction will be used as input into
        samtools view -s for subsampling with SEED 0.
    '''
    if os.path.isfile(flagstatfn) and not os.stat(flagstatfn).st_size == 0:
        with open(flagstatfn, "rt") as flagstatfile:
            nreads = int(next(flagstatfile).split(" ")[0])
        return "{:.4f}".format(targetreads / nreads)
    else:
        return 1.0

# Snakemake rules
wildcard_constraints:
    sample = config['sampleIDconstraint']

localrules: determine_sequencing_type, flag_passedreads, seqdepth, contamMix_align_against_consensus, summary, copy_tmp_to_proj

rule all:
    input:
        expand("{projdir}/{sample}/summary_table.csv", projdir=[PROJDIR], sample=SAMPLES),
        expand("{projdir}/{sample}/{sample}.fa", projdir=[PROJDIR], sample=SAMPLES),

rule prepare_fastq:
    # Convert BAM files to FastQ or link FastQ files 
    output:
        "seqdata/{sample}_1.raw_fastq.gz"
    message: "Prepare FastQ files for sample {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: 0.3
    threads: 1
    params:
        seqdir = config['seqdatadir'],
        datatype = config['seqdatatype'],
        suffix = config['seqdatasuffix'],
        pe2 = "seqdata/{sample}_2.raw_fastq.gz",
        pe0 = "seqdata/{sample}_0.raw_fastq.gz",
        outdir = "seqdata"
    script:
        "scripts/prepare_fastq.py"

rule determine_sequencing_type:
    # Determine if the sequencing data is single-end (1) or paired-end (2) and
    # whether reads have been merged (3)
    input:
        "seqdata/{sample}_1.raw_fastq.gz"
    output:
        "checkpoint/{sample}.seqtype"
    message: "Determine sequencing data type for sample {wildcards.sample}"
    version: 0.3
    threads: 1
    params:
        dir = "seqdata"
    run:
        with open(output[0], "wt") as outfile:
            outfile.write(str(len(glob(f"{params.dir}/{wildcards.sample}*.raw_fastq.gz"))) + "\n")

rule adapter_removal:
    # Run AdapterRemoval to trim adapters and merge overlapping reads
    input:
        fq = "seqdata/{sample}_1.raw_fastq.gz",
        chkpt = "checkpoint/{sample}.seqtype"
    output: 
        pe0 = temp("seqdata/{sample}_0.fastq.gz")
    message: "Remove adapter sequences collapsing overlapping sequences: {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: 0.3
    threads: 8
    log: "logs/adapterremoval/{sample}.adapterremoval.log"
    params: 
        seqtype = lambda wildcards: SAMPLESLIST.loc[wildcards.sample, 'seqdatatype'] if config['seqdatatype'] == "fastq" else check_seqdatatype(f"checkpoint/{sample}.seqtype"),
        input_pe1 = "seqdata/{sample}_1.raw_fastq.gz",
        input_pe2 = "seqdata/{sample}_2.raw_fastq.gz",
        basename = "tmp/{sample}",
        output_pe1 = "seqdata/{sample}_1.fastq.gz",
        output_pe2 = "seqdata/{sample}_2.fastq.gz",
    shell:
        """
        if [[ "{params.seqtype}" = "2" ]]; then
            AdapterRemoval \
                    --file1 {params.input_pe1} \
                    --file2 {params.input_pe2} \
                    --output1 {params.output_pe1} \
                    --output2 {params.output_pe2} \
                    --outputcollapsed {output.pe0} \
                    --settings {log} \
                    --basename {params.basename} \
                    --trimns --trimqualities \
                    --minlength 30 \
                    --minquality 20 \
                    --minadapteroverlap 1 \
                    --collapse \
                    --gzip \
                    --threads {threads}
            rm {params.input_pe2}
        else
            AdapterRemoval \
                    --file1 {input.fq} \
                    --output1 {output.pe0} \
                    --settings {log} \
                    --basename {params.basename} \
                    --trimns --trimqualities \
                    --minlength 30 \
                    --minquality 20 \
                    --gzip \
                    --threads {threads}
        fi
        rm -r {params.basename}*
        """

rule bwa_aln:
    # Align reads against human MT genome with 1,000 bp overhang 
    input:
        "seqdata/{sample}_0.fastq.gz"
    output:
        temp("seqdata/{sample}_MT.{i}.sai")
    message: "Align reads of type {wildcards.i} of sample {wildcards.sample} to only the MT genome with 1000bp overhang with MPI EVA BWA ancient settings (no seed)"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    params:
        seqtype = lambda wildcards: SAMPLESLIST.loc[wildcards.sample, 'seqdatatype'] if config['seqdatatype'] == "fastq" else check_seqdatatype(f"checkpoint/{sample}.seqtype"),
        reffa = f"{workflow.basedir}/resources/NC_012920_1000.fa",
        fastq = "seqdata/{sample}_{i}.fastq.gz"
    threads: 8
    shell:
        """
        if [[ "{params.seqtype}" = "single" ]] && [[ "{wildcards.i}" -ne 0 ]]; then
            touch {output}
        else
            bwa aln \
                -t {threads} \
                -n 0.01 -o 2 -l 16500 \
                -f {output} \
                {params.reffa} \
                {params.fastq}
        fi
        """

rule bwa_sampe:
    # Extract the alignment coordinates of paired reads and convert alignment
    # into BAM format filtering for reads mapped as a proper pair and have a
    # sequence length >= 30 bp
    input:
        sais = lambda wildcards: [f"seqdata/{wildcards.sample}_MT.{i}.sai" for i in range(0, 3)],
        pe0 = "seqdata/{sample}_0.fastq.gz"
    output:
        temp("seqdata/{sample}_MT_12.bam")
    message: "Generate BAM for non-merged reads of sample {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    params: 
        seqtype = lambda wildcards: SAMPLESLIST.loc[wildcards.sample, 'seqdatatype'] if config['seqdatatype'] == "fastq" else check_seqdatatype(f"checkpoint/{sample}.seqtype"),
        reffa = f"{workflow.basedir}/resources/NC_012920_1000.fa",
        readgroup = lambda wildcards: r'@RG\tID:{sample}\tSM:{sample}\tPL:illumina'.format(sample = wildcards.sample),
        input_pe1 = "seqdata/{sample}_1.fastq.gz",
        input_pe2 = "seqdata/{sample}_2.fastq.gz"
    threads: 2
    shell:
        """
        if [[ "{params.seqtype}" = "paired" ]]; then
            bwa sampe \
                -r '{params.readgroup}' \
                -f /dev/stdout \
                {params.reffa} \
                {input.sais[1]} \
                {input.sais[2]} \
                {params.input_pe1} \
                {params.input_pe2} | \
            sambamba view -S \
                -f bam -F 'proper_pair and not unmapped and not mate_is_unmapped and sequence_length >= 30' \
                /dev/stdin > {output} 
            rm {params.input_pe1} {params.input_pe2}
        else
           touch {output}
        fi
        """

rule bwa_samse:
    # Extract the alignment coordinates of single reads and convert alignment
    # into BAM format filtering for reads are mapped and have a sequence length
    # >= 30 bp
    input:
        sais = lambda wildcards: [f"seqdata/{wildcards.sample}_MT.{i}.sai" for i in range(0, 3)],
        pe0 = "seqdata/{sample}_0.fastq.gz"
    output:
        temp("seqdata/{sample}_MT_0.bam")
    message: "Generate BAM for merged reads of sample {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    params:
        reffa = f"{workflow.basedir}/resources/NC_012920_1000.fa",
        readgroup = lambda wildcards: r'@RG\tID:{sample}\tSM:{sample}\tPL:illumina'.format(sample = wildcards.sample)
    threads: 2
    shell:
        """
        bwa samse \
            -r '{params.readgroup}' \
            -f /dev/stdout \
            {params.reffa} \
            {input.sais[0]} \
            {input.pe0} | \
        sambamba view -S \
            -f bam -F 'not paired and not unmapped and sequence_length >= 30' \
            /dev/stdin > {output}
        """

rule bam_merge_wrap_sort:
    # This command merges the paired reads and single reads aligned against the
    # MT genome, wraps the alignment to its original length of 16 569 bp, and
    # sorts the output; the code block using bioawk is necessary to fix the SAM
    # ValidationError raised by HTSJDK used by damageprofiler and dedup
    input:
        bam_12 = "seqdata/{sample}_MT_12.bam",
        bam_0 = "seqdata/{sample}_MT_0.bam"
    output:
        temp("seqdata/{sample}_MTonly.sorted.bam")
    message: "Merge, wrap and sort the aligned reads of {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    params:
        seqtype = lambda wildcards: SAMPLESLIST.loc[wildcards.sample, 'seqdatatype'] if config['seqdatatype'] == "fastq" else check_seqdatatype(f"checkpoint/{sample}.seqtype"),
        reffa = f"{workflow.basedir}/resources/NC_012920_1000.fa",
    threads: 5
    shell:
        """
        if [[ "{params.seqtype}" = "single" ]]; then
            samtools view -bh \
                {input.bam_0} | \
            bam-rewrap MT:16569 | \
            samtools view -h - | \
            bioawk -t '{{
                if ($2 == "0" || $2 == "16") {{gsub("16569","0",$8)}}; print
            }}' - | \
            samtools view -Sb - | \
            samtools sort - | \
            samtools calmd -Q - {params.reffa} > {output}
        else
            samtools merge -c \
                - \
                {input.bam_12} \
                {input.bam_0} | \
            bam-rewrap MT:16569 | \
            samtools view -h - | \
            bioawk -t '{{
                if ($2 == "0" || $2 == "16") {{gsub("16569","0",$8)}}; print
            }}' - | \
            samtools view -Sb - | \
            samtools sort - | \
            samtools calmd -Q - {params.reffa} > {output}
        fi
        """

rule bam_rmdup:
    input:
        "seqdata/{sample}_MTonly.sorted.bam"
    output:
        bam = "seqdata/{sample}_MTonly.sorted.rmdup.bam",
        flagstat = temp("seqdata/{sample}_MTonly.sorted.rmdup.flagstat")
    message: "Remove duplicate reads from {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    log: "logs/dedup/{sample}_dedup.log"
    shell:
        """
        dedup \
            -i {input} \
            --merged \
            -o seqdata
        mv seqdata/{wildcards.sample}_MTonly.sorted.hist logs/dedup/{wildcards.sample}_dedup.hist
        mv seqdata/{wildcards.sample}_MTonly.sorted.log {log}
        samtools sort -o {output.bam} seqdata/{wildcards.sample}_MTonly.sorted_rmdup.bam
        samtools index {output.bam}
        rm seqdata/{wildcards.sample}_MTonly.sorted_rmdup.bam
        samtools flagstat {output.bam} > {output.flagstat}
        """

rule flag_passedreads:
    input:
        "seqdata/{sample}_MTonly.sorted.rmdup.flagstat"
    output:
        "results/{sample}_nReads.flag"
    message: "Determine whether to continue with the processing: {wildcards.sample}"
    run:
        with open(output[0], "wt") as outfile:
            with open(input[0], "rt") as flagstatfile:
                nreads = int(next(flagstatfile).split(" ")[0])
                if nreads >= 500:
                    outfile.write("Pass\n")
                else:
                    outfile.write("Fail\n")

################################################################################

#### Genotype calling mit snpAD ################################################

rule bam2snpAD:
    input: 
        bam = "seqdata/{sample}_MTonly.sorted.rmdup.bam",
        flag = "results/{sample}_nReads.flag"
    output:
        temp("snpAD/{sample}.snpad_input")
    message: "Convert BAM file into snpAD input format for sample {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    params:
        bam2snpAD = f"{workflow.basedir}/resources/snpAD-0.3.3/Bam2snpAD",
        reffasta = f"{workflow.basedir}/resources/NC_012920_1000.fa",
        state = lambda wildcards: check_state(f"results/{wildcards.sample}_nReads.flag")
    shell:
        """
        if [[ {params.state} = "Pass" ]]; then
            {params.bam2snpAD} \
                    -Q 25 \
                    -q 30 \
                    -r MT \
                    -f {params.reffasta} \
                    {input.bam} > {output}
        else
            touch {output}
        fi
        """

rule snpAD_estimation:
    input:
        "snpAD/{sample}.snpad_input"
    output:
        priors = temp("snpAD/{sample}.priors.txt"),
        errors = temp("snpAD/{sample}.errors.txt")
    message: "Estimate the genotype likelihoods using snpAD for sample {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    log: "snpAD/{sample}.snpAD.log"
    params:
        snpAD = f"{workflow.basedir}/resources/snpAD-0.3.3/snpAD",
        state = lambda wildcards: check_state(f"results/{wildcards.sample}_nReads.flag")
    threads: 4
    shell:
        """
        if [[ {params.state} = "Pass" ]]; then
            {params.snpAD} \
                    --cpus={threads} \
                    -o {output.priors} \
                    -O {output.errors} \
                    {input} > {log} 2>&1
        else
            touch {output.priors}
            touch {output.errors}
        fi
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
    version: "0.3"
    params:
        snpADcall = f"{workflow.basedir}/resources/snpAD-0.3.3/snpADCall",
        state = lambda wildcards: check_state(f"results/{wildcards.sample}_nReads.flag")
    shell:
        """
        if [[ {params.state} = "Pass" ]]; then
            {params.snpADcall} \
                    -e <(bioawk -t '$4 < 1e6{{print $1, $2, $3, $4}}' {input.errors}) \
                    -p $(cat {input.priors} |perl -F',' -lane 'for ($i=0 ; $i<4 ; $i++) {{ $F[$i]=$F[$i]/($F[0]+$F[1]+$F[2]+$F[3]); }} for ($i=4; $i<@F; $i++) {{ $F[$i] = "1e-320" ; }} print join ",",@F;') \
                    {input.snpAD} > {output}
        else
            touch {output}
        fi
        """

rule snpAD_vcf2fa:
    input:
        "snpAD/{sample}.snpAD.vcf"
    output:
        fa = "snpAD/{sample}.snpAD.fasta",
        tbi = "snpAD/{sample}.snpAD.vcf.gz.tbi"
    message: "Convert snpAD VCF file into a FastA file for sample {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    params:
        vcf2fa = f"{workflow.basedir}/resources/snpAD-0.3.3/vcf2fasta.pl",
        state = lambda wildcards: check_state(f"results/{wildcards.sample}_nReads.flag")
    shell:
        """
        if [[ {params.state} = "Pass" ]]; then
            {params.vcf2fa} < {input} | \
            sed 's/>MT/>{wildcards.sample}/' > {output.fa}
            bgzip -f {input}
            tabix -f {input}.gz
        else
            touch {output.fa}
            touch {output.tbi}
        fi
        """

################################################################################

#### Down-sample BAM file ######################################################
# Down-sample BAM file for contamMix, mixEMT, and damageprofiler in order to
# not waste a lot of time and resources on deeply sequenced samples

rule namesort_norewrap:
    # Sort by name to be able to revert bam-rewrap's introduction of duplicated
    # reads
    input:
        bam = "seqdata/{sample}_MTonly.sorted.bam",
        flag = "results/{sample}_nReads.flag"
    output:
        temp("seqdata/{sample}_MTonly.nsorted.bam")
    message: "Sort BAM file by name: {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    params:
        state = lambda wildcards: check_state(f"results/{wildcards.sample}_nReads.flag")
    shell:
        """
        if [[ {params.state} = "Pass" ]]; then
            samtools sort -l 0 -n -o {output} {input.bam}
        else
            touch {output}
        fi
        """

rule revert_bamrewrap:
    # Revert bam-wrap
    input:
        "seqdata/{sample}_MTonly.nsorted.bam"
    output:
        temp("seqdata/{sample}_MTonly.nsorted.norerwap.bam")
    message: "Revert bam-wrap by removing multiple occurrences of reads for {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    params:
        state = lambda wildcards: check_state(f"results/{wildcards.sample}_nReads.flag")
    script:
        "scripts/revert_bamrewrap.py"

rule sort_norewrap:
    # Sort the BAM file after reverting bam-wrap by position 
    input:
        "seqdata/{sample}_MTonly.nsorted.norerwap.bam"
    output:
        temp("seqdata/{sample}_MTonly.sorted.norerwap.bam")
    message: "Sort the BAM file after reverting bam-wrap by position: {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    params:
        state = lambda wildcards: check_state(f"results/{wildcards.sample}_nReads.flag")
    shell:
        """
        if [[ {params.state} = "Pass" ]]; then
            samtools sort -o {output} {input}
        else
            touch {output}
        fi
        """

rule flagstat:
    # Determine number of aligned reads for sub-sampling
    input:
        "seqdata/{sample}_MTonly.sorted.norerwap.bam"
    output:
        "logs/flagstat/{sample}.flagstat"
    message: "Determine the number of aligned reads for {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    params:
        state = lambda wildcards: check_state(f"results/{wildcards.sample}_nReads.flag")
    shell:
        """
        if [[ {params.state} = "Pass" ]]; then
            samtools flagstat {input} > {output}
        else
            touch {output}
        fi
        """

################################################################################

#### General stats #############################################################

rule subsample_damageprofiler:
    # Sub-sample BAM file to at most 50,000 aligned reads in order to avoid
    # waiting wasting time and resources
    input:
        bam = "seqdata/{sample}_MTonly.sorted.rmdup.bam",
        flagstat = "logs/flagstat/{sample}.flagstat",
        flag = "results/{sample}_nReads.flag"
    output:
        temp("seqdata/{sample}_MTonly.sorted.rmdup.subsampled.bam")
    message: "Subsample MT aligned reads to 50,000 reads for sample: {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    params:
        subsampling_fraction = lambda wildcards: mixemt_downsampling(f"logs/flagstat/{wildcards.sample}.flagstat", targetreads=50000),
        state = lambda wildcards: check_state(f"results/{wildcards.sample}_nReads.flag")
    shell:
        """
        if [[ {params.state} = "Pass" ]]; then
            samtools view -bh \
                        -s {params.subsampling_fraction} \
                        -o {output} \
                        {input.bam}
        else
            touch {output}
        fi
        """

rule damage_profiler:
    input:
        "seqdata/{sample}_MTonly.sorted.rmdup.subsampled.bam"
    output:
        "qual/{sample}/identity_histogram.pdf"
    message: "Generate damage profile for {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    params: 
        dir = "qual",
        tmpdir = "qual/{sample}_MTonly.sorted.rmdup.subsampled",
        outdir = "qual/{sample}",
        reffa = f"{workflow.basedir}/resources/NC_012920.fa",
        state = lambda wildcards: check_state(f"results/{wildcards.sample}_nReads.flag")
    shell:
        """
        if [[ {params.state} = "Pass" ]]; then
            damageprofiler \
                -i {input} \
                -o {params.dir} \
                -r {params.reffa}
            mv {params.tmpdir}/* {params.outdir}/
            rmdir {params.tmpdir}/
        else
            touch {output}
        fi
        """

rule seqdepth:
    input:
        "seqdata/{sample}_MTonly.sorted.rmdup.bam"
    output:
        "logs/seqdepth/{sample}.seqdepth.csv"
    message: "Run samtools depth to determine the coverage across the MT genome: {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    params: 
        header = lambda wildcards: f"chr\tpos\t{wildcards.sample}"
    shell:
        """
        echo "{params.header}" > {output}
        samtools depth -a -r MT {input} >> {output}
        """

rule haplogrep2:
    input:
        fa = "snpAD/{sample}.snpAD.fasta",
        flag = "results/{sample}_nReads.flag"
    output:
        "logs/haplogrep2/{sample}.hsd"
    message: "Determine mtDNA haplogroup for {wildcards.sample} using HaploGrep2"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    params:
        haplogrep = f"{workflow.basedir}/resources/haplogrep",
        state = lambda wildcards: check_state(f"results/{wildcards.sample}_nReads.flag")
    shell:
        """
        if [[ {params.state} = "Pass" ]]; then
             {params.haplogrep} classify \
                --in {input.fa} \
                --extend-report \
                --format fasta \
                --hits 5 \
                --out {output}
        else
            touch {output}
        fi
        """

################################################################################

#### Contamination estimate using mixEMT #######################################

rule mixemt:
    input:
        bam = "seqdata/{sample}_MTonly.sorted.norerwap.bam",
        flagstat = "logs/flagstat/{sample}.flagstat",
        flag = "results/{sample}_nReads.flag"
    output:
        "logs/mixemt/{sample}.mixemt.log"
    message: "Determine mtDNA contamination for {wildcards.sample} using mixEMT"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    priority: -100
    resources: mixemt = 1
    params: 
        state = lambda wildcards: check_state(f"results/{wildcards.sample}_nReads.flag"),
        subsampling = lambda wildcards: True if float(mixemt_downsampling(f"logs/flagstat/{wildcards.sample}.flagstat")) < 1 else False,
        subsampling_fraction = lambda wildcards: mixemt_downsampling(f"logs/flagstat/{wildcards.sample}.flagstat"),
        subbam = "seqdata/{sample}_MTonly.sorted.rmdup.downsampled.bam",
        mixemtprefix = "logs/mixemt/{sample}.mixemt"
    shell:
        """
        if [[ {params.state} = "Pass" ]]; then
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
            mixemt -q 25 -v -t \
                {params.mixemtprefix} \
                {params.subbam} \
                > {output} \
                2> {params.mixemtprefix}.stderr || \
            if [[ ($(wc -l < {params.mixemtprefix}.log) = "0") && ($(tail -1 {params.mixemtprefix}.stderr) = "0 contributors passed filtering steps.") ]]; then
                echo -e "hap1\tNA\tNA\t0" > {params.mixemtprefix}.log
            else
            rm {params.subbam}* {params.mixemtprefix}.*.tab
            fi
        else
            echo -e "hap1\tNA\tNA\t0" > {params.mixemtprefix}.log
        fi
        """
################################################################################

#### Contamination estimate using contamMix ####################################

rule contamMix_create_sequencePanel:
    input:
        fa = "snpAD/{sample}.snpAD.fasta",
        flag = "results/{sample}_nReads.flag"
    output:
        "logs/contamMix/{sample}/sequence_panel.fasta"
    message: "Align consensus sequence of {wildcards.sample} to panel of 311 modern humans"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    resources: local_jobs = 1
    params: 
        panel = f"{workflow.basedir}/resources/311hu+rCRS.fas",
        state = lambda wildcards: check_state(f"results/{wildcards.sample}_nReads.flag")
    shell:
        """
        if [[ {params.state} = "Pass" ]]; then
            cat {input.fa} {params.panel} | \
            mafft - > {output}
        else
            touch {output}
        fi
        """

rule contamMix_align_against_consensus:
    input:
        flagstat = "logs/flagstat/{sample}.flagstat",
        bam = "seqdata/{sample}_MTonly.sorted.norerwap.bam",
        fas = "snpAD/{sample}.snpAD.fasta",
        flag = "results/{sample}_nReads.flag"
    output:
        "logs/contamMix/{sample}/{sample}.consensus_aligned.bam"
    message: "Align sequences of {wildcards.sample} against its consensus sequence"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    params: 
        state = lambda wildcards: check_state(f"results/{wildcards.sample}_nReads.flag"),
        readgroup = lambda wildcards: r'@RG\tID:{sample}\tSM:{sample}\tPL:illumina'.format(sample = wildcards.sample),
        subsampling = lambda wildcards: True if float(mixemt_downsampling(f"logs/flagstat/{wildcards.sample}.flagstat")) < 1 else False,
        subsampling_fraction = lambda wildcards: mixemt_downsampling(f"logs/flagstat/{wildcards.sample}.flagstat")
    threads: 4
    shell:
        """
        if [[ {params.state} = "Pass" ]]; then
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
                     {input.bam} | \
                samtools sort -n - | \
                samtools fastq \
                    -1 logs/contamMix/{wildcards.sample}/{wildcards.sample}_1.fastq.gz \
                    -2 logs/contamMix/{wildcards.sample}/{wildcards.sample}_2.fastq.gz \
                    -0 logs/contamMix/{wildcards.sample}/{wildcards.sample}_0.fastq.gz \
                    -
            else
                samtools sort -n {input.bam} | \
                samtools fastq \
                    -1 logs/contamMix/{wildcards.sample}/{wildcards.sample}_1.fastq.gz \
                    -2 logs/contamMix/{wildcards.sample}/{wildcards.sample}_2.fastq.gz \
                    -0 logs/contamMix/{wildcards.sample}/{wildcards.sample}_0.fastq.gz \
                    -
            fi
            # BWA aln against sample's consensus sequence
            bwa aln \
                -t {threads} \
                -n 0.01 -o 2 -l 16500 \
                -f logs/contamMix/{wildcards.sample}/{wildcards.sample}.1.sai \
                logs/contamMix/{wildcards.sample}/{wildcards.sample}.fasta \
                logs/contamMix/{wildcards.sample}/{wildcards.sample}_1.fastq.gz
            bwa aln \
                -t {threads} \
                -n 0.01 -o 2 -l 16500 \
                -f logs/contamMix/{wildcards.sample}/{wildcards.sample}.2.sai \
                logs/contamMix/{wildcards.sample}/{wildcards.sample}.fasta \
                logs/contamMix/{wildcards.sample}/{wildcards.sample}_2.fastq.gz
            bwa aln \
                -t {threads} \
                -n 0.01 -o 2 -l 16500 \
                -f logs/contamMix/{wildcards.sample}/{wildcards.sample}.0.sai \
                logs/contamMix/{wildcards.sample}/{wildcards.sample}.fasta \
                logs/contamMix/{wildcards.sample}/{wildcards.sample}_0.fastq.gz
            # BWA sampe and samse to create alignment files
            bwa sampe \
                -r '{params.readgroup}' \
                -f logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT_12.bam \
                logs/contamMix/{wildcards.sample}/{wildcards.sample}.fasta \
                logs/contamMix/{wildcards.sample}/{wildcards.sample}.1.sai \
                logs/contamMix/{wildcards.sample}/{wildcards.sample}.2.sai \
                logs/contamMix/{wildcards.sample}/{wildcards.sample}_1.fastq.gz \
                logs/contamMix/{wildcards.sample}/{wildcards.sample}_2.fastq.gz
            bwa samse \
                -r '{params.readgroup}' \
                -f logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT_0.bam \
                logs/contamMix/{wildcards.sample}/{wildcards.sample}.fasta \
                logs/contamMix/{wildcards.sample}/{wildcards.sample}.0.sai \
                logs/contamMix/{wildcards.sample}/{wildcards.sample}_0.fastq.gz
            # Merge sam files, re-wrap and sort
            if [ -f logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT_120.bam ]; then
                rm logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT_120.bam
            fi
            samtools merge -c \
                - \
                logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT_12.bam \
                logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT_0.bam | \
            bam-rewrap '{wildcards.sample}:16569' | \
            samtools view -bhq 25 - | \
            samtools sort - | \
            samtools calmd -Q - logs/contamMix/{wildcards.sample}/{wildcards.sample}.fasta > {output}
            # Clean
            rm logs/contamMix/{wildcards.sample}/{wildcards.sample}.*.sai \
               logs/contamMix/{wildcards.sample}/{wildcards.sample}_MT* \
               logs/contamMix/{wildcards.sample}/{wildcards.sample}.fasta* \
               logs/contamMix/{wildcards.sample}/{wildcards.sample}_{{0,1,2}}.fastq.gz
        else
            touch {output}
        fi
        """

rule contamMix_estimate:
    input:
        aln = "logs/contamMix/{sample}/sequence_panel.fasta",
        bam = "logs/contamMix/{sample}/{sample}.consensus_aligned.bam"
    output:
        "logs/contamMix/{sample}/contamMix_log.txt"
    message: "Use contamMix to estimate contamination of sample {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    priority: -50
    resources: contamMix = 1
    params:
        estimateR = f"{workflow.basedir}/resources/contamMix/exec/estimate.R",
        state = lambda wildcards: check_state(f"results/{wildcards.sample}_nReads.flag")
    threads: 36
    shell:
        """
        if [[ {params.state} = "Pass" ]]; then
            Rscript {params.estimateR} \
                    --samFn {input.bam} \
                    --malnFn {input.aln} \
                    --nIter 15000 \
                    --tabOutput TRUE > {output}
        else
            touch {output}
        fi
        """

################################################################################

#### Summary ###################################################################

rule summary:
    input: 
        flagstat = "logs/flagstat/{sample}.flagstat",
        seqdepth = "logs/seqdepth/{sample}.seqdepth.csv",
        damageprofiler = "qual/{sample}/identity_histogram.pdf",
        mixemt = "logs/mixemt/{sample}.mixemt.log",
        contammix = "logs/contamMix/{sample}/contamMix_log.txt",
        snpAD = "snpAD/{sample}.snpAD.fasta",
        haplogrep = "logs/haplogrep2/{sample}.hsd"
    output: "{projdir}/{sample}/summary_table.csv"
    message: "Summarise the results in a table: {wildcards.sample}"
    conda: f"{workflow.basedir}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    params:
        readlengthdist = "qual/{sample}/lgdistribution.txt"
    script:
        "scripts/summarise_sample.R"

rule copy_tmp_to_proj:
    input:
        "{projdir}/{sample}/summary_table.csv"
    output:
        "{projdir}/{sample}/{sample}.fa"
    message: "Copy results from tmp to project folder: {wildcards.sample}"
    version: "0.3"
    shell:
        """
        mkdir -p {PROJDIR}/{wildcards.sample}
        # Copy non-clipped BAM file
        cp -r seqdata/{wildcards.sample}_MTonly.sorted.rmdup.bam* {PROJDIR}/{wildcards.sample}/
        # DeDup logs
        cp -r logs/dedup/{wildcards.sample}_dedup.* {PROJDIR}/{wildcards.sample}/
        # HaploGrep
        cp logs/haplogrep2/{wildcards.sample}.hsd {PROJDIR}/{wildcards.sample}/
        # MixEMT
        cp logs/mixemt/{wildcards.sample}.* {PROJDIR}/{wildcards.sample}/
        # DamageProfiler
        cp qual/{wildcards.sample}/3pGtoA_freq.txt {PROJDIR}/{wildcards.sample}/{wildcards.sample}_3pGtoA_freq.txt
        cp qual/{wildcards.sample}/5pCtoT_freq.txt {PROJDIR}/{wildcards.sample}/{wildcards.sample}_5pCtoT_freq.txt
        cp qual/{wildcards.sample}/DamagePlot.pdf {PROJDIR}/{wildcards.sample}/{wildcards.sample}_DamagePlot.pdf
        cp qual/{wildcards.sample}/Length_plot.pdf {PROJDIR}/{wildcards.sample}/{wildcards.sample}_Length_plot.pdf
        cp qual/{wildcards.sample}/lgdistribution.txt {PROJDIR}/{wildcards.sample}/{wildcards.sample}_lgdistribution.txt
        cp qual/{wildcards.sample}/misincorporation.txt {PROJDIR}/{wildcards.sample}/{wildcards.sample}_misincorporation.txt
        # snpAD
        cp snpAD/{wildcards.sample}.snpAD.vcf.gz* {PROJDIR}/{wildcards.sample}/
        cp snpAD/{wildcards.sample}.snpAD.fasta {PROJDIR}/{wildcards.sample}/{wildcards.sample}.fa
        """

# Clean temporary output

rule clean_tmp:
    params:
        tmpdir = config['tmpdir']
    shell:
        "rm -r {params.tmpdir}"

