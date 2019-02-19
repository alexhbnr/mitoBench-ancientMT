##############################################################################
# MT pipeline of the Stoneking group
#
# authors: Enrico Macholdt and Alex Huebner
# date: 30/01/18
##############################################################################

from glob import glob
import os.path
from datetime import datetime
import subprocess
import pysam

workdir: config['tmpdir']


def print_report():
    ''' Generate report about files and software versions used
    '''
    logstr = []
    logstr.append("\n# MtDNA processing report\n")
    logstr.append("sequencing run: {}".format(config['bam']))
    logstr.append("Processing started: {}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    git_version = subprocess.Popen("cd /mnt/genotyping/sk_pipelines/source/mtPipeline; git rev-parse HEAD", shell=True, stdout=subprocess.PIPE).stdout.read().decode().rstrip()
    logstr.append("Git version of Snakefile: {}\n".format(git_version))
    logstr.append("### Version of programs used")
    samtools_version = subprocess.Popen("{} --version | head -1 | cut -f2 -d' '".format(config['SAMTOOLS']), shell=True, stdout=subprocess.PIPE).stdout.read().decode().rstrip()
    logstr.append("* samtools: {}".format(samtools_version))
    bamfixpair_version = subprocess.Popen("{} --version  2>&1 >/dev/null | cut -f3 -d' '".format(config['BAMFIXPAIR']), shell=True, stdout=subprocess.PIPE).stdout.read().decode().rstrip()
    logstr.append("* bam-fixpair: {}".format(bamfixpair_version))
    bamrmdup_version = subprocess.Popen("{} --version  2>&1 >/dev/null | cut -f3 -d' '".format(config['BAMRMDUP']), shell=True, stdout=subprocess.PIPE).stdout.read().decode().rstrip()
    logstr.append("* bam-rmdup: {}".format(bamrmdup_version))
    bwa_version = subprocess.Popen("{} --version  2>&1 >/dev/null | cut -f3 -d' '".format(config['BAMRMDUP']), shell=True, stdout=subprocess.PIPE).stdout.read().decode().rstrip()
    logstr.append("* bwa: {}".format(bwa_version))
    picard_version = subprocess.Popen("java -jar {} AddOrReplaceReadGroups --version 2>&1".format(config['PICARD']), shell=True, stdout=subprocess.PIPE).stdout.read().decode().rstrip()
    logstr.append("* Picard: {}".format(picard_version))
    gatk_version = subprocess.Popen("java -jar {} --version".format(config['GATK']), shell=True, stdout=subprocess.PIPE).stdout.read().decode().rstrip()
    logstr.append("* GATK: {}".format(gatk_version))
    bcftools_version = subprocess.Popen("{} --version | head -1 | cut -f2 -d' '".format(config['BCFTOOLS']), shell=True, stdout=subprocess.PIPE).stdout.read().decode().rstrip()
    logstr.append("* BCFTOOLS: {}".format(bcftools_version))
    mtConsensus_version = subprocess.Popen("cd /home/alexander_huebner/Documents/PythonScripts/; git rev-parse HEAD".format(config['MTCONSENSUS']), shell=True, stdout=subprocess.PIPE).stdout.read().decode().rstrip()
    logstr.append("* vcf.MTconsensus.py: {}".format(mtConsensus_version))
    haplogrep_version = subprocess.Popen("java -jar {} --version | head -1 | cut -f9 -d' ' | sed 's/[()]//g'".format(config['HAPLOGREP']), shell=True, stdout=subprocess.PIPE).stdout.read().decode().rstrip()
    logstr.append("* HAPLOGREP: {}".format(haplogrep_version))
    logstr.append("\n")
    return "\n".join(logstr)


onstart:
    print(print_report())

onsuccess:
    print("Workflow finished, no error")
    shell("mail -s ""Snakemake finished without an error"" {config[email]}  < {log} ")

onerror:
    print("An error occurred")
    shell("mail -s ""an error occurred""  {config[email]} < {log}")

localrules: extract_MTreads_and_filter, sort_MTreads, remove_duplicates, split_by_rg, sequencing_depth, concatenate_fastas, backup


RGS = [rg['ID'] for rg in pysam.AlignmentFile(config['bam']).header['RG'] if rg['ID'] != "control"]
#RGS = [line.rstrip() for line in open("rg.txt")]

wildcard_constraints:
    id = "[A-Za-z0-9]+"

rule all:
    input:  expand("{id}.hsd", id=RGS),
            expand("{id}.mixemt.log", id=RGS),
            "consensus.fasta",
            "sequencing_depth.txt"


rule extract_MTreads_and_filter:
    # Extract MT reads and known NUMTs and filter for QC fail
    output: "mt_numts.bam"
    message: "Filtering of reads (MT, Numts)"
    threads: 1
    shell:
        """
        {config[BAMFIXPAIR]} --quiet {config[bam]} | \
        {config[SAMTOOLS]} view -bh -L {config[NUMTCATALOGUE]} -F 512 - > {output} 2> /dev/null
        """


rule sort_MTreads:
    # Sort MT and known NUMTs reads
    input:  "mt_numts.bam"
    output: temp("mt_numts.sorted.bam")
    message: "Sort MT reads for duplicate removal"
    threads: 1
    shell:
        """
        {config[SAMTOOLS]} sort -o {output} {input}
        """


rule remove_duplicates:
    # Remove duplicates with BAM-RMDUP
    input:  "mt_numts.sorted.bam"
    output: temp("mt_numts.rmdup.bam")
    message: "Duplicate removal using Bam-Rmdup"
    log: "rmdup.log"
    threads: 1
    shell:
        """
        {config[BAMRMDUP]} -o {output} -p {input} > {log}
        """


rule split_by_rg:
    # Split by read group
    input:  "mt_numts.rmdup.bam" 
    output: expand("{id}.bam", id=RGS)
    message: "Split BAM file by read group"
    threads: 1
    shell:
        """
        {config[SAMTOOLS]} split -f "%!.%." -u unassigned.bam {input}
        rm unassigned.bam
        """


rule convert_to_FQ:
    # Convert all BAM files to FASTQ
    input:  "{id}.bam"
    output: pe0="{id}.0.fq",
            pe1="{id}.1.fq",
            pe2="{id}.2.fq"
    message: "Convert read group {wildcards.id} to FastQ"
    threads: 1
    shell:
        """
        {config[BAMFIXPAIR]} --kill-widows {input} \
                > {wildcards.id}.fixed.bam \
                2> /dev/null
        {config[SAMTOOLS]} fastq -0 {output.pe0} \
                                 -1 {output.pe1} \
                                 -2 {output.pe2} \
                                 {wildcards.id}.fixed.bam
        rm {wildcards.id}.fixed.bam
        """


rule align_BWAmem:
    # Align FastQ files with BWAmem to MT
    input: pe0="{id}.0.fq",
           pe1="{id}.1.fq",
           pe2="{id}.2.fq"
    output: "{id}.sorted.bam"
    message: "Align sample {wildcards.id} with BWA mem to the MT reference genome"
    params: bwaParams = "-B 6 -O 3 -E 2"
    threads: 4
    shell:
        """
        {config[BWA]} mem -t {threads} {params.bwaParams} {config[MTREFOVERLAP]}.fasta \
                {input.pe1} {input.pe2} | \
                {config[SAMTOOLS]} view -Sb - > {wildcards.id}.pe.bam
        {config[BWA]} mem -t {threads} {params.bwaParams} {config[MTREFOVERLAP]}.fasta \
                {input.pe0} | \
                {config[SAMTOOLS]} view -Sb - > {wildcards.id}.se.bam
        {config[SAMTOOLS]} cat -h <({config[SAMTOOLS]} view -H {wildcards.id}.pe.bam) \
                -o /dev/stdout \
                {wildcards.id}.pe.bam \
                {wildcards.id}.se.bam | \
            {config[BAMREWRAP]} MT:16569 | \
            {config[SAMTOOLS]} sort -o {output} -
        rm {wildcards.id}.pe.bam {wildcards.id}.se.bam
        """


rule add_readgroup:
    # Add RG tag to sequencing data
    input: "{id}.sorted.bam"
    output: bam=temp("{id}.RGtag.bam"),
            bai=temp("{id}.RGtag.bam.bai")
    message: "Add RG tag to aligned reads of {wildcards.id}"
    threads: 1
    shell:
        """
        java -Xmx2g -jar {config[PICARD]} \
                AddOrReplaceReadGroups \
                I={input} O={output.bam} \
                SO="coordinate" \
                RGID={wildcards.id} \
                RGLB={wildcards.id} \
                RGPL=illumina \
                RGPU=unknown \
                RGSM={wildcards.id}
        {config[SAMTOOLS]} index {output.bam}
        """


rule indel_realignment:
    # Perform indel realignment using GATK
    input: bam="{id}.RGtag.bam",
           bai="{id}.RGtag.bam.bai"
    output: "{id}.realigned.bam"
    message: "Indel-realignment of sample {wildcards.id}"
    threads: 1
    shell:
        """
        java -Xmx3g -jar {config[GATK]} \
                -T RealignerTargetCreator \
                -R {config[MTREF]}.fasta \
                -o {wildcards.id}.intervals \
                -I {input.bam}
        cat {wildcards.id}.intervals {config[KNOWNINDELS]} > {wildcards.id}-extra.intervals
        java -Xmx3g -jar {config[GATK]} \
                -T IndelRealigner \
                -targetIntervals {wildcards.id}-extra.intervals \
                -o {output} \
                -R {config[MTREF]}.fasta \
                -I {input.bam}
        rm {wildcards.id}.intervals {wildcards.id}-extra.intervals
        """


rule vcf:
    # Generate pileup and call consensus using BCFTOOLS
    input:  "{id}.realigned.bam"
    output: "{id}.vcf.gz"
    message: "Generate pileup and call a consensus sequencing for sample {wildcards.id}"
    params: mq = 50,
            depth = 1000,
            fields = "DP,AD,ADF,ADR"
    threads: 1
    shell:
        """
        {config[SAMTOOLS]} mpileup -C {params.mq} \
                                   -d {params.depth} \
                                   -f {config[MTREF]}.fasta \
                                   -v -t {params.fields} \
                                   -o {wildcards.id}.gl.vcf \
                                   {input}
        {config[BCFTOOLS]} call -o {output} \
                                -A -O z --ploidy 1 -m \
                                {wildcards.id}.gl.vcf
        rm {wildcards.id}.gl.vcf
        """


rule consensus_fasta:
    # Generate FastA file and report about heteroplasmies and indels
    input:  "{id}.vcf.gz"
    output: fasta = "{id}.fasta",
            hsd = temp("{id}.vcf")
    message: "Generate consensus FastA sequence for sample {wildcards.id}"
    threads: 1
    shell:
        """
        {config[MTCONSENSUS]} -i {input} \
                              --hetreport {wildcards.id}.heteroplasmy.txt \
                              --indreport {wildcards.id}.indel.txt \
                              --hsd {output.hsd} \
                              --indqual {config[INDQUAL]} \
                              -o {output.fasta}
        """


rule concatenate_fastas:
    input: expand("{id}.fasta", id=RGS)
    output: "consensus.fasta"
    message: "Concatenate all FastA files"
    shell:
        """
        cat *.fasta > {output}
        """


rule haplogrep:
    input: "{id}.vcf"
    output: "{id}.hsd"
    message: "Haplogroup calling for sample {wildcards.id}"
    threads: 1
    shell:
        """
        java -jar {config[HAPLOGREP]} --format vcf \
                                      --in {input} \
                                      --out {output} \
                                      --phylotree {config[PHYLOTREE]}
        """


rule flagstat:
    input: "{id}.realigned.bam"
    output: temp("{id}.flagstat")
    message: "Run flagstat analysis for sample {wildcards.id}"
    shell:
        "{config[SAMTOOLS]} flagstat {input} > {output}"


rule mixemt:
    input:  "{id}.flagstat"
    output: log="{id}.mixemt.log",
            pat="{id}.mixemt.pos.tab"
    message: "Run mixEMT analysis for sample {wildcards.id}"
    log: "{id}.mixemt.stderr"
    threads: 1
    shell:
        """
        export PYENV_VERSION={config[pyenv]}
        READS=$(head -1 {input} | cut -f1 -d' ')
        if [[ ${{READS}} -gt {config[NREADS]} ]]; then
            FRAC=$(perl -E "say {config[NREADS]} / ${{READS}}")
            {config[SAMTOOLS]} view -bh \
                    -s ${{FRAC}} \
                    -o {wildcards.id}-sub.bam {wildcards.id}.bam
            {config[SAMTOOLS]} index {wildcards.id}-sub.bam
            python2 {config[MIXEMT]} -v \
                    -t {wildcards.id}.mixemt \
                    {wildcards.id}-sub.bam > {output.log} 2> {log}
        else
            python2 {config[MIXEMT]} -v \
                    -t {wildcards.id}.mixemt \
                    {wildcards.id}.realigned.bam > {output.log} 2> {log}
        fi
        """


rule sequencing_depth:
    # Generate coverage analysis for the whole MT genome per sample
    input: expand("{id}.realigned.bam", id=RGS)
    output: "sequencing_depth.txt"
    message: "Generate coverage analysis for the whole MT genome per sample"
    params: header = "\t".join(["chr", "pos"] + [os.path.basename(fn).replace(".realigned.bam", "") for fn in glob("*.realigned.bam")]),
            bams = " ".join([os.path.basename(fn) for fn in glob("*.realigned.bam")])
    threads: 1
    shell:
        """
        echo -e "{params.header}" > {output}
        {config[SAMTOOLS]} depth -a -r MT {params.bams} >> {output}
        """


rule backup:
    # Copy all files over to /mnt/genotyping, generate overall FastA file including all samples and concatenate HaploGrep2 results
    input:  hsd=expand("{id}.hsd", id=RGS),
            mixemt=expand("{id}.mixemt.log", id=RGS),
            cons="consensus.fasta",
            seqdepth="sequencing_depth.txt"
    output: touch("pipeline.done")
    message: "Backup files to {config[projdir]}"
    threads: 1
    shell:
        """
        if [[ -d {config[projdir]} ]]; then
          cp rmdup.log {config[projdir]}/
          mv consensus.fasta {config[projdir]}/
          cp sequencing_depth.txt {config[projdir]}/
          mkdir -p {config[projdir]}/heteroplasmy
          cp *.heteroplasmy.txt {config[projdir]}/heteroplasmy/
          mkdir -p {config[projdir]}/indel
          cp *.indel.txt {config[projdir]}/indel/
          mkdir -p {config[projdir]}/fasta
          cp *.fasta {config[projdir]}/fasta/
          mkdir -p {config[projdir]}/bam
          cp *.realigned.bam {config[projdir]}/bam/
          rename 's/.realigned//g' {config[projdir]}/bam/*.realigned.bam
          mkdir -p {config[projdir]}/vcf
          cp *.vcf.gz {config[projdir]}/vcf/
          mkdir -p {config[projdir]}/haplogrep
          cp *.hsd {config[projdir]}/haplogrep/
          awk 'FNR==1 && NR!=1 {{ while (/^SampleID/) getline; }} 1 {{print}}' {config[projdir]}/haplogrep/*.hsd \
                  > {config[projdir]}/haplogroup_calls.hsd
          mkdir -p {config[projdir]}/mixemt
          cp *.mixemt.{{log,pos.tab,stderr}} {config[projdir]}/mixemt/
          {config[HGCONT]} --path {config[projdir]}/heteroplasmy/ \
                           --phylotree {config[PHYLOTREEPOS]} \
                           -o {config[projdir]}/haplogroup_contamination.txt
        else
          echo "Project directory {config[projdir]} not found. The directory has to exist before running backup."
          exit 125
        fi
        """


rule clean:
    shell:
        """
        rm -r {config[tmpdir]}
        """
