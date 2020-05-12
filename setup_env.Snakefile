################################################################################
# Snakefile for preparing all the resources before being able to run the
# mitoBench ancient DNA MT pipeline.
#
# Alex Huebner, 01/03/19
################################################################################

from snakemake.utils import R

workdir: "/tmp"

# Constants
PATH = workflow.basedir
HUMAN_MT_FAS = f"{PATH}/resources/NC_012920.fa.gz"


localrules: decompress_fasta, bwa_index, samtools_index


rule all:
    input: 
        f"{PATH}/resources/NC_012920.fa.fai",
        f"{PATH}/resources/NC_012920_1000.fa.fai",
        f"{PATH}/resources/NC_012920_1000.fa.ann",
        f"{PATH}/resources/haplogrep",
        f"{PATH}/resources/contamMix/exec/estimate.R",
        f"{PATH}/resources/install_contamMix.done",
        f"{PATH}/resources/install_summary.done"

# Prepare MT genome for analysis

rule decompress_fasta:
    output:
        "{PATH}/resources/NC_012920.fa"
    message: "De-compress the FastA sequence of the human MT genome with 1000 bp extension"
    conda: f"{PATH}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    shell:
        "gunzip -c {HUMAN_MT_FAS} > {output}"

rule extend_fasta:
    input:
        "{PATH}/resources/NC_012920.fa"
    output:
        "{PATH}/resources/NC_012920_1000.fa"
    conda: f"{PATH}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    shell:
        """
        bioawk -c fastx '{{
            print ">" $name; \
            print $seq substr($seq,1,1000); 
        }}' {input} > {output}
        """

rule bwa_index:
    input:
        "{PATH}/resources/NC_012920_1000.fa"
    output:
        "{PATH}/resources/NC_012920_1000.fa.ann"
    message: "BWA index the FastA sequence of the human MT genome with 1000 bp extension"
    conda: f"{PATH}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    shell:
        "bwa index {input}"

rule samtools_index_extended:
    input:
        "{PATH}/resources/NC_012920_1000.fa"
    output:
        "{PATH}/resources/NC_012920_1000.fa.fai"
    message: "Samtools faidx the FastA sequence of the human MT genome with 1000 bp extension"
    conda: f"{PATH}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    shell:
        "samtools faidx {input}"

rule samtools_index:
    input:
        "{PATH}/resources/NC_012920.fa"
    output:
        "{PATH}/resources/NC_012920.fa.fai"
    message: "Samtools faidx the FastA sequence of the human MT genome"
    conda: f"{PATH}/env/mitoBench_bioconda.yaml"
    version: "0.3"
    shell:
        "samtools faidx {input}"

# Download software not available via conda

rule download_haplogrep:
    output:
        "{PATH}/resources/haplogrep"
    message: "Download haplogrep-cmd from GitHub"
    params: url = "https://github.com/seppinho/haplogrep-cmd/releases/download/v2.2.5/haplogrep.zip"
    shell:
        """
        wget -O {PATH}/resources/haplogrep.zip {params.url}
        unzip {PATH}/resources/haplogrep.zip -d {PATH}/resources
        rm {PATH}/resources/haplogrep.zip
        """

rule uncompress_contamMix:
    output:
        "{PATH}/resources/contamMix/exec/estimate.R"
    message: "Uncompress the tar ball of contamMix"
    conda: f"{PATH}/env/mitoBench_bioconda.yaml"
    params:
        tarball = f"{workflow.basedir}/resources/contamMix_1.0-10.tar.gz" 
    shell:
        "tar xvf {params.tarball}"

rule install_contammix:
    output:
        "{PATH}/resources/install_contamMix.done"
    message: "Install R package of contamMix"
    conda: f"{PATH}/env/mitoBench_bioconda.yaml"
    params:
        tarball = f"{workflow.basedir}/resources/contamMix_1.0-10.tar.gz" 
    script:
        "scripts/install_contamMix.R" 

rule install_summary:
    output:
        "{PATH}/resources/install_summary.done"
    message: "Install R package used for the summary"
    conda: f"{PATH}/env/mitoBench_bioconda.yaml"
    script:
        "scripts/install_summary.R" 
