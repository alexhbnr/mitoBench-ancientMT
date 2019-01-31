################################################################################
# Ancient MT pipeline - Schmutzi
# 
# Schmutzi (Renaud et al., 2016) is a tool to estimate the amount of present-day
# contamination in ancient DNA samples that was introduced at different stages
# during the processing of the sample, e.g. excavation, DNA extraction, or
# library preparation. Schmutzi is based on a Bayesian maximum a posteriori
# algorithm and is able to both estimate the amount of contamination and call a
# MT consensus sequence from the reads that are inferred as endogenous.
#
# Using Schmutzi for consensus calling gives us two advantages:
#   1. We do not have to worry about filtering reads for deamination patterns or
#      mask these patterns for subsequent analysis because Schmutzi includes all
#      reads into its analysis and evaluates the deamination pattern using an
#      error model.
#   2. The contamination estimate is included into the analysis and we do not
#      have to debate whether to use all or deaminated-only reads because
#      Schmutzi aims to call a consensus only on the reads that it inferred as
#      being endogenous.
#
# Alex Huebner, 05/10/18
################################################################################

from glob import glob
import os
import pandas as pd


workdir: config['tmpdir']

#### SAMPLES ###################################################################
SAMPLES = [line.rstrip() for line in open(config['sampleprefixes'], "rt")]
BAMS = {sample: "{}/{}.{}".format(config['bamdir'], sample, config['bamsuffix'])
        for sample in SAMPLES}
################################################################################

#### TOOLS #####################################################################
SAMTOOLS = "/mnt/genotyping/sk_pipelines/source/Ychr_bin/samtools-1.9"
BAMFIXPAIR = "/mnt/genotyping/sk_pipelines/source/Ychr_bin/bam-fixpair"
BWA = "/home/public/usr/bin/bwa"
BAMREWRAP = "/mnt/genotyping/sk_pipelines/source/Ychr_bin/bam-rewrap"
BAMSAMPLE = "/mnt/genotyping/sk_pipelines/source/bin/bam-sample/bam-sample.py"
VCF2FASTA = "/mnt/genotyping/sk_pipelines/source/bin/vcf2fasta.py"
FA2HSD = "/mnt/genotyping/sk_pipelines/source/bin/mtDNA.haplogroup.py"
################################################################################


#### SAMPLE PREPARATION ########################################################

localrules: link_index, evaluate_substprofile

rule sample_preparation:
    input:  "schmutzi/samples_passedQC.txt",
            "schmutzi/samples_failedQC.txt"

rule link_index:
    # Link the raw data to file and index the BAM file in order to retrieve only
    # the MT reads
    output:  bam = "rawdata/{sample}.bam",
             bai = "rawdata/{sample}.bam.bai"
    message: "Link the BAM file to the temp directory and index it: {wildcards.sample}"
    params: bam = lambda wildcards: BAMS[wildcards.sample]
    shell:
        """
        ln -s {params.bam} ${{PWD}}/rawdata/{wildcards.sample}.bam
        {SAMTOOLS} index rawdata/{wildcards.sample}.bam
        """

rule extract_MT_reads:
    input:  bam = "rawdata/{sample}.bam",
            bai = "rawdata/{sample}.bam.bai"
    output: "bam/{sample}.MT.bam"
    message: "Extract the MT reads from sample {wildcards.sample}"
    params:  region = "MT"
    shell:
        """
        {SAMTOOLS} view -bh \
                -o {output} {input.bam} {params.region}
        """

rule realign_MT_long:
    input:  "bam/{sample}.MT.bam"
    output: bam = "bam/{sample}.MT_long.bam",
            bai = "bam/{sample}.MT_long.bam.bai"
    message: "Align {wildcards.sample} to the human MT reference genome with a 1000 overlap with ancient settings"
    threads: 8
    shell:
        """
        {BAMFIXPAIR} {input} | \
        {BWA} bam2bam \
                -g /mnt/solexa/Genomes/human_MT_long/bwa-0.4.9 \
                -n 0.01 -o 2 -l 16500 \
                -t {threads} - | \
        {BAMREWRAP} 'gi|251831106|ref|NC_012920.1|:16569' | \
        {SAMTOOLS} sort -o {output.bam} -
        {SAMTOOLS} index {output.bam}
        """

rule schmutzi_prof:
    # BAM2prof generates a frequency file of observed substitutions compared
    # to the MT genome; for samples with little data "-nan" values are written
    # that lead to subsequent failures in running Schmutzi. Therefore, I create
    # these files manually to use them as check points for subsequent analysis.
    input:  bam = "bam/{sample}.MT_long.bam",
            bai = "bam/{sample}.MT_long.bam.bai"
    output: end5 = "schmutzi/{sample}.endo.5p.prof",
            end3 = "schmutzi/{sample}.endo.3p.prof"
    message: "Run bam2prof on sample {wildcards.sample}"
    params:  lengthDeam = config['lengthDeam'],
             library =  config['library'],
             reffasta = "/mnt/solexa/Genomes/human_MT/whole_genome.fa"
    shell:
        """
        /mnt/genotyping/sk_pipelines/source/bin/schmutzi/src/bam2prof \
                -length {params.lengthDeam} \
                -endo \
                -{params.library} \
                -5p {output.end5} \
                -3p {output.end3} {input.bam}
        """

rule evaluate_substprofile:
    # Evaluate the substitution profiles (*.prof) of a sample for the presence
    # of NaN values that are an indicator for the lack of data that is
    # required for subsequent analysis like contamination estimation or
    # consensus calling. If we find samples that failed, we will flag these
    # samples and run majority consensus calling ignoring C > T and G > A
    # substitutions on either forward or reverse strand instead.
    input:  expand("schmutzi/{sample}.endo.5p.prof", sample=SAMPLES)
    output: passed="schmutzi/samples_passedQC.txt",
            failed="schmutzi/samples_failedQC.txt"
    message: "Evaluate the substitution profiles for lack of data"
    run:
        PASSED = []
        FAILED = []
        for sm in SAMPLES:
            p3_profile = pd.read_csv("schmutzi/{}.endo.3p.prof".format(sm),
                                     sep="\t")
            p5_profile = pd.read_csv("schmutzi/{}.endo.5p.prof".format(sm),
                                     sep="\t")
            if (p3_profile.isnull().values.sum() > 0 or 
                p5_profile.isnull().values.sum() > 0):
                FAILED.append(sm)
            else:
                PASSED.append(sm)
        with open(output.passed, "wt") as outfile:
            for sample in PASSED:
                outfile.write(sample + "\n")
        with open(output.failed, "wt") as outfile:
            for sample in FAILED:
                outfile.write(sample + "\n")

################################################################################

#### CONSENSUS CALLING #########################################################

SCHMUTZI_SAMPLES = []
BAMSAMPLE_SAMPLES = []
if (os.path.isfile("schmutzi/samples_passedQC.txt") and
        os.path.isfile("schmutzi/samples_failedQC.txt")):
    SCHMUTZI_SAMPLES = [line.rstrip() for line in open("schmutzi/samples_passedQC.txt")]
    BAMSAMPLE_SAMPLES = [line.rstrip() for line in open("schmutzi/samples_failedQC.txt")]

wildcard_constraints:
    schmutzi = config['constraints'],
    bamsample = config['constraints']

rule consensus_calling:
    input:  expand("schmutzi/{schmutzi}-npred_final_endo.fa", schmutzi=SCHMUTZI_SAMPLES),
            expand("schmutzi/{schmutzi}-npred_final_endo.hsd", schmutzi=SCHMUTZI_SAMPLES),
            expand("schmutzi/{bamsample}_bamsample.vcf.gz", bamsample=BAMSAMPLE_SAMPLES+SCHMUTZI_SAMPLES),
            expand("schmutzi/{bamsample}_angsd.fa", bamsample=BAMSAMPLE_SAMPLES+SCHMUTZI_SAMPLES)

rule contDeam:
    output: "schmutzi/{schmutzi}.cont.est"
    message: "Run ContDeam on schmutzi {wildcards.schmutzi}"
    params:  bam = lambda wildcards: "bam/{schmutzi}.MT_long.bam".format(schmutzi=wildcards.schmutzi),
             lengthDeam = config['lengthDeam'],
             library =  config['library'],
             reffasta = "/mnt/solexa/Genomes/human_MT/whole_genome.fa"
    shell:
        """
        /mnt/genotyping/sk_pipelines/source/bin/schmutzi/src/contDeam.pl \
                --lengthDeam {params.lengthDeam} \
                --library {params.library} \
                --out schmutzi/{wildcards.schmutzi} \
                {params.reffasta} {params.bam}
        """

rule schmutzi_nopred:
    input: "schmutzi/{schmutzi}.cont.est"
    output: "schmutzi/{schmutzi}-npred_final_endo.fa"
    message: "Run schmutzi without predicting the contamination: {wildcards.schmutzi}"
    params: reffasta = "/mnt/solexa/Genomes/human_MT/whole_genome.fa",
            bam = lambda wildcards: "bam/{schmutzi}.MT_long.bam".format(schmutzi=wildcards.schmutzi)
    shell:
        """
        /mnt/genotyping/sk_pipelines/source/bin/schmutzi/src/schmutzi.pl \
                --uselength \
                --notusepredC \
                --name {wildcards.schmutzi} \
                --namec {wildcards.schmutzi}_cont \
                --logindel 1 \
                --ref {params.reffasta} \
                --out schmutzi/{wildcards.schmutzi}-npred \
                schmutzi/{wildcards.schmutzi} \
                /mnt/genotyping/sk_pipelines/source/bin/schmutzi/share/schmutzi/alleleFreqMT/eurasian/freqs/ \
                {params.bam}
        """

rule haplogrep:
    input: "schmutzi/{schmutzi}-npred_final_endo.fa"
    output: "schmutzi/{schmutzi}-npred_final_endo.hsd"
    message: "Run Haplogroup 2: {wildcards.schmutzi}"
    shell:
        """
        java -Xmx3g -jar /mnt/genotyping/sk_pipelines/source/bin/haplogrep/haplogrep-2.1.16.jar \
                --format fasta \
                --in {input} \
                --out {output} \
                --phylotree 17
        """

rule bam_sample:
    output: "schmutzi/{bamsample}_bamsample.vcf.gz"
    message: "Majority consensus call ignoring C>T substitutions: {wildcards.bamsample}"
    params: bam = lambda wildcards: "bam/{bamsample}.MT_long.bam".format(bamsample=wildcards.bamsample),
            reffasta = "/mnt/solexa/Genomes/human_MT/whole_genome.fa",
            mincov = config['mincov']
    shell:
        """
        {BAMSAMPLE} --bam {params.bam} \
                    --ref {params.reffasta} \
                    --format VCF \
                    --sample-name {wildcards.bamsample} \
                    --method majority \
                    --strand-check non-USER_all \
                    --minbq 30 \
                    --mincov {params.mincov} | bgzip > {output}
        """

rule angsd_consensus:
    # Call consensus sequence as described by Ehler et al. (2018): amtDB 
    output: "schmutzi/{bamsample}_angsd.fa"
    message: "Call consensus sequence using ANGSD following Ehler et al. (2018): {wildcards.bamsample}"
    params: bam = lambda wildcards: "bam/{bamsample}.MT_long.bam".format(bamsample=wildcards.bamsample),
            reffasta = "/mnt/solexa/Genomes/human_MT/whole_genome.fa",
            dir = "schmutzi"
    shell:
        """
        /mnt/genotyping/sk_pipelines/source/bin/angsd-0.923/angsd \
                -i {params.bam} \
                -minMapQ 30 \
                -minQ 20 \
                -doFasta 2 \
                -doCounts 1 \
                -ref {params.reffasta} \
                -out {params.dir}/{wildcards.bamsample}_angsd.tmp
        /mnt/genotyping/sk_pipelines/source/bin/bioawk \
                -c fastx '{{print ">{wildcards.bamsample}"; print $seq}}' \
                {params.dir}/{wildcards.bamsample}_angsd.tmp.fa.gz > {output}
        """

################################################################################


#### Back-up & Clean ###########################################################

localrules: backup, clean

rule backup:
    message: "Back-up files to project directory"
    params: projdir = config['projdir']
    shell:
        """
        mkdir -p {params.projdir}
        echo "Move BAM files with MT reads only aligned to MT reference genome\n"
        mkdir -p {params.projdir}/bam
        rsync -av bam/*MT_long.bam* {params.projdir}/bam/
        echo "Move bam_sample.py VCFs"
        mkdir -p {params.projdir}/vcf
        rsync -av schmutzi/*_bamsample.vcf.gz {params.projdir}/vcf/
        echo "Move Schmutzi's consensus FastAs"
        mkdir -p {params.projdir}/fasta
        rsync -av schmutzi/*-npred_final_endo.{{fa,log}} {params.projdir}/fasta/
        rsync -av schmutzi/*_angsd.fa {params.projdir}/fasta/
        echo "Move Schmutzi's contamination estimates"
        mkdir -p {params.projdir}/contamination
        rsync -av schmutzi/*-npred_final.cont.{{est,pdf}} {params.projdir}/contamination/
        echo "Move Haplogrep's haplogroup assignment"
        mkdir -p {params.projdir}/hsd
        rsync -av schmutzi/*-npred_final_endo.hsd {params.projdir}/hsd/
        """

rule clean:
    message: "Delete the temporary folder"
    params: tmpdir = config['tmpdir']
    shell:
        """
        rm -r {params.tmpdir}
        """

################################################################################
