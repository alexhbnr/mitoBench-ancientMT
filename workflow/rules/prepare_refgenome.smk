rule prepare_refgenome:
    input:
        reffa = "{tmpdir}/refgenome/refgenome.fa.fai",
        reffa_ext = "{tmpdir}/refgenome/refgenome_ext.fa.fai",
        index = "{tmpdir}/refgenome/index_building.done"
    output:
        touch("{tmpdir}/refgenome/prepare_refgenome.done")

rule copy_reffa:
    output:
        "{tmpdir}/refgenome/refgenome.fa"
    message: "Copy the reference genome into tmp dir"
    params:
        reffa = config['reffa']
    shell:
        "cp {params.reffa} {output}"

rule index_reffa:
    input:
        "{tmpdir}/refgenome/refgenome.fa"
    output:
        "{tmpdir}/refgenome/refgenome.fa.fai"
    message: "Index the reference genome using samtools"
    wrapper:
        "v1.3.2/bio/samtools/faidx"

rule extend_reffa:
    input:
        "{tmpdir}/refgenome/refgenome.fa"
    output:
        "{tmpdir}/refgenome/refgenome_ext.fa"
    message: "Extend the genome for the alignment due to circular nature of mtDNA"
    params:
        extension = config['circular_extension']
    run:
        with open(output[0], "wt") as outfile:
            for name, seq in pyfastx.Fasta(input[0], build_index=False):
                seq = seq + seq[:1000]
                outfile.write(f">{name}\n{seq}\n")

rule index_extend_reffa:
    input:
        "{tmpdir}/refgenome/refgenome_ext.fa"
    output:
        "{tmpdir}/refgenome/refgenome_ext.fa.fai"
    message: "Index the extended reference genome using samtools"
    wrapper:
        "v1.3.2/bio/samtools/faidx"

if "bwa" in config['aligner']:

    rule bwa_index_reffa:
        input:
            "{tmpdir}/refgenome/refgenome_ext.fa"
        output:
            dummy = touch("{tmpdir}/refgenome/index_building.done"),
            idx = "{tmpdir}/refgenome/refgenome_ext.amb"
        message: "Index the extended reference genome using BWA"
        params:
            prefix = "{tmpdir}/refgenome/refgenome_ext"
        wrapper:
            "v1.3.2/bio/bwa/index"

elif config['aligner'] == "bowtie2":

    rule bowtie2_index_reffa:
        input:
            "{tmpdir}/refgenome/refgenome_ext.fa"
        output:
            touch("{tmpdir}/refgenome/index_building.done")
        message: "Index the extended reference genome using BowTie2"
        wrapper:
            "v1.3.2/bio/bowtie2/build"
