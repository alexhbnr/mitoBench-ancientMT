def determine_se_fastq_files(wildcards):
    sdtype = int(open(checkpoints.determine_seqdata_type.get(**wildcards).output[0],
                      "rt").readline().rstrip())
    fns = []
    if sdtype % 2 == 1:
        fns.append(f"{config['tmpdir']}/seqdata/{wildcards.sample}_0_se.fastp.fastq.gz")
    if sdtype > 1:
        fns.append(f"{config['tmpdir']}/seqdata/{wildcards.sample}_0_pe.fastp.fastq.gz")
    if config['merge_overlapping']:
        fns.append(f"{config['tmpdir']}/seqdata/{wildcards.sample}_0_mg.fastp.fastq.gz")
    return fns


rule adapterremoval:
    input:
        lambda wildcards: [f"{config['tmpdir']}/seqdata/{wildcards.sample}_{i}.fastp.fastq.gz" for i in range(3)]
    output:
        touch("{tmpdir}/seqdata/{sample}.adapterremoval")
    shell:
        "touch {output}"

checkpoint determine_seqdata_type:
    input:
        "{tmpdir}/seqdata/{sample}.prepare_fastq"
    output:
        sdtype = "{tmpdir}/seqdata/{sample}.sdtype"
    message: "Determine the sequencing data type (paired-end/single-end/both): {wildcards.sample}"
    resources:
        mem_gb = 2
    params:
        dir = "{tmpdir}/seqdata"
    run:
        fqfns = [fq
                 for fq in glob(f"{params.dir}/{wildcards.sample}*.raw_fastq.gz")
                 if re.search(rf"{wildcards.sample}_[0-2].raw_fastq.gz",
                              os.path.basename(fq)) and
                    os.stat(fq).st_size > 50]
        with open(output[0], "wt") as outfile:
            outfile.write(str(len(fqfns)) + "\n")

if config['merge_overlapping']:

    rule fastp_pe_merge:
        input:
            dummy = "{tmpdir}/seqdata/{sample}.prepare_fastq",
            sample = lambda wildcards: [f"{wildcards.tmpdir}/seqdata/{wildcards.sample}_{i}.raw_fastq.gz" for i in range(1, 3)]
        output:
            trimmed = ["{tmpdir}/seqdata/{sample}_1.fastp.fastq.gz", "{tmpdir}/seqdata/{sample}_2.fastp.fastq.gz"],
            unpaired = temp("{tmpdir}/seqdata/{sample}_0_pe.fastp.fastq.gz"),
            merged = temp("{tmpdir}/seqdata/{sample}_0_mg.fastp.fastq.gz"),
            html = "{tmpdir}/logs/fastp/{sample}_pe.html",
            json = "{tmpdir}/logs/fastp/{sample}_pe.json"
        message: "Process the paired-end sequencing data using FastP: {wildcards.sample}"
        resources:
            mem_gb = 4
        params:
            extra = "-z 4 --merge --detect_adapter_for_pe -l 30 -q 20 -w 8 -V" if config['fastp_adaptertrimming'] else "-z 4 --merge -A -l 30 -q 20 -w 8 -V"
        threads: 8
        wrapper:
            "v1.3.2/bio/fastp"

else:

    rule fastp_pe_nonmerge:
        input:
            dummy = "{tmpdir}/seqdata/{sample}.prepare_fastq",
            sample = lambda wildcards: [f"{wildcards.tmpdir}/seqdata/{wildcards.sample}_{i}.raw_fastq.gz" for i in range(1, 3)]
        output:
            trimmed = ["{tmpdir}/seqdata/{sample}_1.fastp.fastq.gz", "{tmpdir}/seqdata/{sample}_2.fastp.fastq.gz"],
            unpaired = temp("{tmpdir}/seqdata/{sample}_0_pe.fastp.fastq.gz"),
            html = "{tmpdir}/logs/fastp/{sample}_pe.html",
            json = "{tmpdir}/logs/fastp/{sample}_pe.json"
        message: "Process the paired-end sequencing data using FastP: {wildcards.sample}"
        resources:
            mem_gb = 4
        params:
            extra = "-z 4 --detect_adapter_for_pe -l 30 -q 20 -w 8 -V" if config['fastp_adaptertrimming'] else "-z 4 -A -l 30 -q 20 -w 8 -V"
        threads: 8
        wrapper:
            "v1.3.2/bio/fastp"

rule fastp_se:
    input:
        sample = ["{tmpdir}/seqdata/{sample}_0.raw.fastq.gz"]
    output:
        trimmed = temp("{tmpdir}/seqdata/{sample}_0_se.fastp.fastq.gz"),
        html = "{tmpdir}/logs/fastp/{sample}_se.html",
        json = "{tmpdir}/logs/fastp/{sample}_se.json"
    message: "Process the single-end sequencing data using FastP: {wildcards.sample}"
    resources:
        mem_gb = 4
    params:
        extra = "-z 4 -A"
    threads: 4
    wrapper:
        "v1.3.2/bio/fastp"

rule concat_singleend_fastp:
    input:
        lambda wildcards: determine_se_fastq_files(wildcards)
    output:
        "{tmpdir}/seqdata/{sample}_0.fastp.fastq.gz"
    message: "Concatenate the single-end quality filtered sequencing data: {wildcards.sample}"
    resources:
        mem_gb = 2
    shell:
        "cat {input} > {output}"
