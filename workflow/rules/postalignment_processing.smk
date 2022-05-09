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

rule samtools_flagstat:
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

checkpoint flag_passedreads:
    input:
        bai = expand("{tmpdir}/seqdata/{sample}_MTonly.sorted.rmdup.bam.bai", tmpdir=[config['tmpdir']], sample=SAMPLES),
        flagstat = expand("{tmpdir}/logs/flagstat/{sample}.flagstat", tmpdir=[config['tmpdir']], sample=SAMPLES),
        depth = expand("{tmpdir}/logs/depth/{sample}.depth", tmpdir=[config['tmpdir']], sample=SAMPLES)
    output:
        "{tmpdir}/logs/initial_alignment.tsv"
    message: "Determine whether to continue with the processing"
    run:
        # Read flagstat results
        flagstat = {}
        for fn in input.flagstat:
            with open(fn, "rt") as flagstatfile:
                nreads = int(next(flagstatfile).split(" ")[0])
            sample = os.path.basename(fn).replace(".flagstat", "")
            flagstat[sample] = nreads
        flagstat = pd.DataFrame.from_dict(flagstat, orient="index", columns=['nReads']) \
            .reset_index() \
            .rename({'index': 'sample'}, axis=1)

        # Read depth summary 
        depth = pd.concat([pd.read_csv(fn, sep="\t", header=None,
                                       names=['chrom', 'pos', 'depth'])
                           .assign(sample=os.path.basename(fn).replace(".depth", ""))
                           for fn in input.depth])
        depth_summary =  depth \
            .groupby(['sample', 'chrom'])['depth'].agg(['mean', 'std', 'median']) \
            .reset_index() \
            .drop(['chrom'], axis=1)
        depth_summary['median'] = depth_summary['median'].astype(int)

        # Merge
        stats = flagstat.merge(depth_summary, how="left", on="sample")
        stats['status'] = ['pass' if x >= config['n_aln_reads'] else 'fail'
                           for x in stats['nReads']]
        stats[['sample', 'status', 'nReads', 'median', 'mean', 'std']] \
            .sort_values(['sample']) \
            .to_csv(output[0], sep="\t", index=False, float_format="%.6f")
