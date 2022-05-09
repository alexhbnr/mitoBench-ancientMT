rule samtools_calmd:
    input:
        aln = "{tmpdir}/seqdata/{sample}_MTonly.aligned.bam",
        ref = "{tmpdir}/refgenome/refgenome.fa"
    output:
        pipe("{tmpdir}/seqdata/{sample}_MTonly.calmd.bam")
    message: "Calculate the MD field: {wildcards.sample}"
    resources:
        mem_gb = 4
    params:
        extra = "-u"
    wrapper:
        "v1.3.2/bio/samtools/calmd"

rule samtools_sort_coord:
    input:
        "{tmpdir}/seqdata/{sample}_MTonly.calmd.bam"
    output:
        "{tmpdir}/seqdata/{sample}_MTonly.sorted.bam"
    message: "Sort aligned reads by coordinate: {wildcards.sample}"
    resources:
        mem_gb = 8
    threads: 4
    wrapper:
        "v1.3.2/bio/samtools/sort"

if config['duplicate_removal']:

    rule picard_markduplicates:
        input:
            "{tmpdir}/seqdata/{sample}_MTonly.sorted.bam"
        output:
            bam = "{tmpdir}/seqdata/{sample}_MTonly.sorted.rmdup.bam",
            metrics = "{tmpdir}/logs/markduplicates/{sample}.markduplicates.log"
        message: "Mark duplicate reads with Picard MarkDuplicates: {wildcards.sample}"
        resources:
            mem_gb = 12
        params:
            extra = "--REMOVE_DUPLICATES true --ASSUME_SORTED true"
        wrapper:
            "v1.3.2/bio/picard/markduplicates"

else:

    rule link_dedup_bam:
        input:
            "{tmpdir}/seqdata/{sample}_MTonly.sorted.bam"
        output:
            "{tmpdir}/seqdata/{sample}_MTonly.sorted.rmdup.bam"
        message: "Link the sorted BAM file because we skip the duplicate removal step: {wildcards.sample}"
        resources:
            mem_gb = 2
        run:
            if not os.path.islink(output[0]):
                if wildcards.tmpdir[0] != "/":
                    os.symlink(f"{os.getcwd()}/{input[0]}", output[0])
                else:
                    os.symlink(input[0], output[0])

rule samtools_index_rmdup:
    input:
        "{tmpdir}/seqdata/{sample}_MTonly.sorted.rmdup.bam"
    output:
        "{tmpdir}/seqdata/{sample}_MTonly.sorted.rmdup.bam.bai"
    message: "Index the final BAM file: {wildcards.sample}"
    wrapper:
        "v1.3.2/bio/samtools/index"

checkpoint samtools_flagstat:
    input:
        "{tmpdir}/seqdata/{sample}_MTonly.sorted.rmdup.bam"
    output:
        "{tmpdir}/logs/flagstat/{sample}.flagstat"
    message: "Run samtools flagstat: {wildcards.sample}"
    wrapper:
        "v1.3.2/bio/samtools/flagstat"

rule samtools_depth:
    input:
        bams = "{tmpdir}/seqdata/{sample}_MTonly.sorted.rmdup.bam"
    output:
        "{tmpdir}/logs/depth/{sample}.depth"
    message: "Run samtools depth: {wildcards.sample}"
    params:
        extra = "-a"
    wrapper:
        "v1.3.2/bio/samtools/depth"

checkpoint summarise_depth:
    input:
        "{tmpdir}/logs/depth/{sample}.depth"
    output:
        "{tmpdir}/logs/depth/{sample}.avg_depth"
    message: "Summarise the depth: {wildcards.sample}"
    run:
        summary = pd.read_csv(input[0], sep="\t", header=None,
                    names=['chrom', 'pos', 'depth']) \
            .groupby(['chrom'])['depth'].agg(['mean', 'std', 'median']) \
            .reset_index()
        summary['median'] = summary['median'].astype(int)
        summary.to_csv(output[0], sep="\t", index=False, float_format="%.2f")

checkpoint flag_passedreads:
    input:
        bai = "{tmpdir}/seqdata/{sample}_MTonly.sorted.rmdup.bam.bai",
        flagstat = "{tmpdir}/logs/flagstat/{sample}.flagstat",
        depth = "{tmpdir}/logs/depth/{sample}.avg_depth"
    output:
        "{tmpdir}/results/{sample}_nReads.flag"
    message: "Determine whether to continue with the processing: {wildcards.sample}"
    run:
        with open(output[0], "wt") as outfile:
            with open(input['flagstat'], "rt") as flagstatfile:
                nreads = int(next(flagstatfile).split(" ")[0])
                if nreads >= 500:
                    outfile.write("Pass\n")
                else:
                    outfile.write("Fail\n")
