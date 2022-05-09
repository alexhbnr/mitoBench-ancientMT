def return_expected_haplogrep2_res(wildcards):
    status = pd.read_csv(checkpoints.flag_passedreads.get(**wildcards).output[0],
                         sep="\t")
    passed_samples = status.loc[status['status'] == "pass"]['sample'].tolist()
    print([f"{config['tmpdir']}/logs/haplogrep/{sample}.haplogrep2.tsv" for sample in passed_samples])
    return [f"{config['tmpdir']}/logs/haplogrep/{sample}.haplogrep2.tsv" for sample in passed_samples]


rule haplogroup:
    input:
        lambda wildcards: return_expected_haplogrep2_res(wildcards)
    output:
        touch("{tmpdir}/haplogroups.done")

rule download_haplogrep:
    output:
        touch("{tmpdir}/haplogrep.installed") 
    message: "Download and install HaploGrep2"
    params:
        url = "https://github.com/seppinho/haplogrep-cmd/releases/download/v2.4.0/haplogrep.zip",
        path = f"{workflow.basedir}/resources/haplogrep.zip"
    shell:
        """
        wget -O {params.path} {params.url} && unzip {params.path}
        """

rule haplogrep_classify:
    input:
        dummy = "{tmpdir}/haplogrep.installed",
        vcf = "{tmpdir}/consensus_sequence/{sample}.vcf.gz"
    output:
        "{tmpdir}/logs/haplogrep/{sample}.haplogrep2.tsv"
    message: "Infer the most likely mtDNA haplogroups: {wildcards.sample}"
    params:
        path = f"{workflow.basedir}/resources/haplogrep"
    shell:
        """
        {params.path} classify --extend-report --format=vcf --hits=5 --write-fasta \
            --in {input.vcf} \
            --out {output}
        """
