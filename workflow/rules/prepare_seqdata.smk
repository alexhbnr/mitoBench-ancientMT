if config['seqdatatype'] == "fastq":

    rule link_fastq_files:
        output:
            touch("{tmpdir}/seqdata/{sample}.prepare_fastq")
        message: "Link the FastQ files into the correct folder: {wildcards.sample}"
        resources:
            mem_gb = 2
        params:
            seqdatadir = config['seqdatadir'],
            seqdatasuffix = config['seqdatasuffix'],
            outdir = "{tmpdir}/seqdata"
        run:
            fqfns = [fn for fn in glob(f"{params.seqdatadir}/{wildcards.sample}*.{params.seqdatasuffix}")
                     if re.search(rf"{wildcards.sample}(_[0-2])*\.{params.seqdatasuffix}", os.path.basename(fn))]
            if len(fqfns) == 0:
                print(f"No FastQ files found for sample {wildcards.sample}. "
                      "FastQ files need to have the format "
                      "'SAMPLENAME_[012].fastq.gz'")
                sys.exit(1)
            else:
                for fqfn in fqfns:
                    if any([f"_{i}" in os.path.basename(fqfn) for i in range(3)]):
                        if not os.path.islink(f"{params.outdir}/"
                                              f"{os.path.basename(fqfn).replace('.' + params.seqdatasuffix, '')}"
                                              ".raw_fastq.gz"):
                            os.symlink(f"{fqfn}", f"{params.outdir}/"
                                       f"{os.path.basename(fqfn).replace('.' + params.seqdatasuffix, '')}"
                                       ".raw_fastq.gz")
                    else:
                        if not os.path.islink(f"{params.outdir}/"
                                              f"{os.path.basename(fqfn).replace('.' + params.seqdatasuffix, '')}"
                                              "_1.raw_fastq.gz"):
                            os.symlink(f"{fqfn}", f"{params.outdir}/"
                                       f"{os.path.basename(fqfn).replace('.' + params.seqdatasuffix, '')}"
                                       "_1.raw_fastq.gz")
            Path(output[0]).touch()

elif config['seqdatatype'] == "bam":

    rule samtools_nsort_bamfixpair:
        input:
            lambda wildcards: f"{config['seqdatadir']}/{wildcards.sample}.{config['seqdatasuffix']}"
        output:
            pipe("{tmpdir}/seqdata/{sample}.nsorted.bam")
        message: "Sort the reads of the BAM file by name: {wildcards.sample}"
        resources:
            mem_gb = 8
        params:
            extra = "-nu"
        threads: 4
        wrapper:
            "v1.32/bio/samtools/sort"

    rule samtools_fixmate_bamfixpair:
        input:
            "{tmpdir}/seqdata/{sample}.nsorted.bam"
        output:
            pipe("{tmpdir}/seqdata/{sample}.fixmate.bam")
        message: "Fix the mate flags: {wildcards.sample}"
        resources:
            mem_gb = 2
        threads: 1
        wrapper:
            "v1.32/bio/samtools/fixmate"

    rule samtools_fastq:
        input:
            "{tmpdir}/seqdata/{sample}.fixmate.bam"
        output:
            pe1 = "{tmpdir}/seqdata/{sample}_1.raw_fastq.gz",
            pe2 = "{tmpdir}/seqdata/{sample}_2.raw_fastq.gz",
            pe0 = "{tmpdir}/seqdata/{sample}_0.raw_fastq.gz"
        message: "Convert to FastQ: {wildcards.sample}"
        resources:
            mem_gb = 2
        params:
            params_fastq = "-N -c 0"
        threads: 1
        wrapper:
            "https://github.com/alexhbnr/snakemake-wrappers/raw/bio/samtools/fastq/separate"

    rule prepare_fastq_from_bam:
        input:
            pe1 = "{tmpdir}/seqdata/{sample}_1.raw_fastq.gz",
            pe2 = "{tmpdir}/seqdata/{sample}_2.raw_fastq.gz",
            pe0 = "{tmpdir}/seqdata/{sample}_0.raw_fastq.gz"
        output:
            touch("{tmpdir}/seqdata/{sample}.prepare_fastq")

