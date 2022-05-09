def return_expected_consensusseqs(wildcards):
    flag = open(checkpoints.flag_passedreads.get(**wildcards).output[0],
                "rt").readline().rstrip()
    if flag == "Pass":
        return f"{wildcards.tmpdir}/consensus_sequence/{wildcards.sample}.fa"
    else:
        return ""


def chrname_refgenome(fn):
    for name, seq in pyfastx.Fasta(fn, build_index=False):
        break
    return f"{name}:1-{len(seq)}"


def chrlength_refgenome(fn):
    for name, seq in pyfastx.Fasta(fn, build_index=False):
        break
    return len(seq)


def return_subsampling_frac_gt(wildcards):
    depth_summary = pd.read_csv(checkpoints.summarise_depth.get(**wildcards).output[0], sep="\t")
    median = depth_summary.iloc[0]['median']
    if median > config['genotype_median']:
        return f"-s {config['genotype_median']/median:.6f}"
    else:
        return ""


if config['genotype_downsampling']:

    rule downsample_genotyping:
        input:
            "{tmpdir}/seqdata/{sample}_MTonly.sorted.rmdup.bam"
        output:
            "{tmpdir}/seqdata/{sample}_MTonly.gt.bam"
        message: "Down-sample BAM file prior to genotyping: {wildcards.sample}"
        params:
            extra = lambda wildcards: f"-bh {return_subsampling_frac_gt(wildcards)}"
        wrapper:
            "v1.3.2/bio/samtools/view"

else:

    rule link_bam_genotyping:
        input:
            "{tmpdir}/seqdata/{sample}_MTonly.sorted.rmdup.bam"
        output:
            "{tmpdir}/seqdata/{sample}_MTonly.gt.bam"
        message: "Link BAM file prior to genotyping: {wildcards.sample}"
        params:
            fn = lambda wildcards: f"{wildcards.tmpdir}/seqdata/{wildcards.sample}_MTonly.sorted.rmdup.bam" if wildcards.tmpdir.startswith("/") else f"{os.getcwd()}/{wildcards.tmpdir}/seqdata/{wildcards.sample}_MTonly.sorted.rmdup.bam"
        shell:
            "ln -s {params.fn} {output}"


if config['genotyper'] == "snpAD":

    rule bam2snpAD:
        input: 
            bam = "{tmpdir}/seqdata/{sample}_MTonly.gt.bam",
            reffa = "{tmpdir}/refgenome/refgenome.fa"
        output:
            temp("{tmpdir}/snpAD/{sample}.snpad_input")
        message: "Convert BAM file into snpAD input format: {wildcards.sample}"
        resources:
            mem_gb = 8
        params:
            bam2snpAD = "~/github/mitoBench-ancientMT/resources/snpAD-0.3.9/Bam2snpAD/Bam2snpAD",
            bq = config['genotype_minbq'],
            mq = config['genotype_minmq'],
            chrname = chrname_refgenome(config['reffa'])
        shell:
            """
            {params.bam2snpAD} \
                -Q {params.mq} \
                -q {params.bq} \
                -r {params.chrname} \
                -f {input.reffa} \
                {input.bam} > {output}
            """

    rule snpAD_estimation:
        input:
            "{tmpdir}/snpAD/{sample}.snpad_input"
        output:
            priors = temp("{tmpdir}/snpAD/{sample}.priors.txt"),
            errors = temp("{tmpdir}/snpAD/{sample}.errors.txt")
        message: "Estimate the genotype likelihoods using snpAD: {wildcards.sample}"
        log: "{tmpdir}/snpAD/{sample}.snpAD.log"
        params:
            snpAD = "~/github/mitoBench-ancientMT/resources/snpAD-0.3.9/snpAD/snpAD",
        threads: 4
        shell:
            """
            {params.snpAD} \
                --cpus={threads} \
                -o {output.priors} \
                -O {output.errors} \
                {input} > {log} 2>&1
            """

    rule snpAD_modify_errors:
        input:
            "{tmpdir}/snpAD/{sample}.errors.txt"
        output:
            "{tmpdir}/snpAD/{sample}.errors.mod.txt"
        message: "Adapt the errors for haploid genotype calls: {wildcards.sample}"
        shell:
            """
            bioawk -t '$4 < 1e6{{print $1, $2, $3, $4}}' {input} > {output}
            """

    rule snpAD_modify_priors:
        input:
            "{tmpdir}/snpAD/{sample}.priors.txt"
        output:
            "{tmpdir}/snpAD/{sample}.priors.mod.txt"
        message: "Adapt the priors for haploid genotype calls: {wildcards.sample}"
        run:
            with open(output[0], "wt") as outfile:
                with open(input[0], "rt") as infile:
                    priors = infile.readline().rstrip().split(",")
                    sum_homozygous_priors = sum([float(p) for p in priors[:4]])
                    mod_priors = [f"{(float(p) / sum_homozygous_priors):.6f}"
                                  if i < 4 else "1e-320"
                                  for i, p in enumerate(priors)]
                    outfile.write(",".join(mod_priors) + "\n")

    rule snpAD_call:
        input:
            snpAD = "{tmpdir}/snpAD/{sample}.snpad_input",
            priors = "{tmpdir}/snpAD/{sample}.priors.mod.txt",
            errors = "{tmpdir}/snpAD/{sample}.errors.mod.txt"
        output:
            pipe("{tmpdir}/consensus_sequence/{sample}.vcf")
        message: "Call the genotypes using snpAD fixing the likelihood of a heterozygous genotype to a very small number: {wildcards.sample}"
        params:
            snpADcall = "~/github/mitoBench-ancientMT/resources/snpAD-0.3.9/snpAD/snpADCall",
        shell:
            """
            {params.snpADcall} \
                -e {input.errors} \
                -p {input.priors} \
                --name={wildcards.sample} \
                {input.snpAD} > {output}
            """

    rule snpAD_bgzip:
        input:
            "{tmpdir}/consensus_sequence/{sample}.vcf"
        output:
            "{tmpdir}/consensus_sequence/{sample}.vcf.gz"
        message: "Compress the snpAD VCF file: {wildcards.sample}"
        wrapper:
            "v1.3.2/bio/bgzip"

elif config['genotyper'] == "freeBayes":

    rule freebayes:
        input:
            samples = ["{tmpdir}/seqdata/{sample}_MTonly.gt.bam"],
            ref = "{tmpdir}/refgenome/refgenome.fa"
        output:
            "{tmpdir}/consensus_sequence/{sample}.vcf.gz"
        message: "Genotype with freeBayes: {wildcards.sample}"
        params:
            extra = f"--report-monomorphic -p 1 -i -X -u -q {config['genotype_minbq']}"
        wrapper:
            "v1.3.2/bio/freebayes"

#elif config['genotyper'] == "bcftools":

