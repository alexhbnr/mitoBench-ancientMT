if config['aligner'] == "bwaaln":

    rule bwa_aln_summary:
        input:
            "{tmpdir}/seqdata/{sample}_MTonly.aligned.bam"
        output:
            touch("{tmpdir}/seqdata/{sample}.initial_alignment")
        shell:
            "touch {output}"

    rule bwa_aln_initial:
        input:
            adapterremoval = "{tmpdir}/seqdata/{sample}.adapterremoval",
            prepare_refgenome = "{tmpdir}/refgenome/prepare_refgenome.done",
            fastq = "{tmpdir}/seqdata/{sample}_{pair}.fastp.fastq.gz",
            idx = "{tmpdir}/refgenome/refgenome_ext.amb"
        output:
            temp("{tmpdir}/seqdata/{sample}_{pair}.sai")
        message: "Align with BWA aln: {wildcards.sample} for data type {wildcards.pair}"
        resources:
            mem = 8
        params:
            extra = config['bwa_aln_params']
        threads: 8
        wrapper:
            "v1.3.2/bio/bwa/aln"

    rule bwa_sampe_initial:
        input:
            fastq = ["{tmpdir}/seqdata/{sample}_1.fastp.fastq.gz", "{tmpdir}/seqdata/{sample}_2.fastp.fastq.gz"],
            sai = ["{tmpdir}/seqdata/{sample}_1.sai", "{tmpdir}/seqdata/{sample}_2.sai"],
            idx = "{tmpdir}/refgenome/refgenome_ext.amb"
        output:
            pipe("{tmpdir}/seqdata/{sample}_MTonly_12.bam")
        message: "Generate BAM file for the paired-end data: {wildcards.sample}"
        resources:
            mem = 4
        params:
            sort = "none",
            extra = lambda wildcards: "-r " + r'@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina'.format(sample = wildcards.sample)
        threads: 1
        wrapper:
            "v1.3.2/bio/bwa/sampe"

    rule samtools_fixmate_sampe:
        input:
            "{tmpdir}/seqdata/{sample}_MTonly_12.bam"
        output:
            pipe("{tmpdir}/seqdata/{sample}_MTonly_12.fixmate.bam")
        message: "Fix the mate flags: {wildcards.sample}"
        resources:
            mem_gb = 2
        threads: 1
        wrapper:
            "v1.3.2/bio/samtools/fixmate"

    rule bwa_sampe_initiial_filtering:
        input:
            "{tmpdir}/seqdata/{sample}_MTonly_12.fixmate.bam"
        output:
            temp("{tmpdir}/seqdata/{sample}_MTonly_12.aligned.bam")
        message: "Only keep aligned read pairs: {wildcards.sample}"
        resources:
            mem_gb = 2
        params:
            extra = "-uh -e 'flag.proper_pair && !flag.unmap && !flag.munmap'"
        wrapper:
            "v1.3.2/bio/samtools/view"

    rule bwa_samse_initial:
        input:
            fastq = "{tmpdir}/seqdata/{sample}_0.fastp.fastq.gz",
            sai = "{tmpdir}/seqdata/{sample}_0.sai",
            idx = "{tmpdir}/refgenome/refgenome_ext.amb"
        output:
            pipe("{tmpdir}/seqdata/{sample}_MTonly_0.bam")
        message: "Generate BAM file for the single-end data: {wildcards.sample}"
        resources:
            mem = 4
        params:
            sort = "none",
            extra = lambda wildcards: "-r " + r'@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina'.format(sample = wildcards.sample)
        threads: 1
        wrapper:
            "v1.3.2/bio/bwa/samse"

    rule bwa_samse_initiial_filtering:
        input:
            "{tmpdir}/seqdata/{sample}_MTonly_0.bam"
        output:
            temp("{tmpdir}/seqdata/{sample}_MTonly_0.aligned.bam")
        message: "Only keep aligned single reads: {wildcards.sample}"
        resources:
            mem_gb = 2
        params:
            extra = "-uh -e '!flag.paired && !flag.unmap'"
        wrapper:
            "v1.3.2/bio/samtools/view"

    rule bwa_aln_initial_merge:
        input:
            pe = "{tmpdir}/seqdata/{sample}_MTonly_12.aligned.bam",
            se = "{tmpdir}/seqdata/{sample}_MTonly_0.aligned.bam"
        output:
            "{tmpdir}/seqdata/{sample}_MTonly.aligned.bam"
        message: "Merge the aligned paired and single reads: {wildcards.sample}"
        resources:
            mem_gb = 4
        params:
            extra = "-nu"
        wrapper:
            "v1.3.2/bio/samtools/merge"

    #rule bam_rewrap:
        ## TODO: get on conda and wrap into a container
        #input:
            #"{tmpdir}/seqdata/{sample}_MTonly.aligned.bam"
        #output:
            #"{tmpdir}/seqdata/{sample}_MTonly.rewrap.bam"
        #message: "Revert the circularisation of the MT genome: {wildcards.sample}"
        #resources:
            #mem_gb = 4
        #shell:
            #"""
            #samtools view -bh {input} | \
            #bam-rewrap MT:16569 > {output}
            #"""
            ##samtools view -h - | \
            ##bioawk -t '{{
                ##if ($2 == "0" || $2 == "16") {{gsub("16569","0",$8)}}; print
            ##}}' - | \
            ##samtools view -Sb - > {output}

elif config['aligner'] == "bwamem": 

    rule bwa_mem_summary:
        input:
            "{tmpdir}/seqdata/{sample}_MTonly.aligned.bam"
        output:
            touch("{tmpdir}/seqdata/{sample}.initial_alignment")
        shell:
            "touch {output}"

    rule bwa_mem_initial_pe:
        input:
            adapterremoval = "{tmpdir}/seqdata/{sample}.adapterremoval",
            prepare_refgenome = "{tmpdir}/refgenome/prepare_refgenome.done",
            reads = ["{tmpdir}/seqdata/{sample}_1.fastp.fastq.gz", "{tmpdir}/seqdata/{sample}_2.fastp.fastq.gz"],
            idx = "{tmpdir}/refgenome/refgenome.amb"
        output:
            pipe("{tmpdir}/seqdata/{sample}_MTonly_12.aligned.sam")
        message: "Align with BWA mem: {wildcards.sample}"
        resources:
            mem_gb = 8
        params:
                extra=r"-R '@RG\tID:{sample}\tSM:{sample}'"
        threads: 8
        wrapper:
            "v1.3.2/bio/bwa/mem"

    rule bwa_mem_initial_pe_bam:
        input:
            "{tmpdir}/seqdata/{sample}_MTonly_12.aligned.sam"
        output:
            temp("{tmpdir}/seqdata/{sample}_MTonly_12.aligned.bam")
        message: "Convert BWA mem PE output to BAM format: {wildcards.sample}"
        params:
            extra = "-Su"
        wrapper:
            "v1.3.2/bio/samtools/view"

    rule bwa_mem_initial_se:
        input:
            adapterremoval = "{tmpdir}/seqdata/{sample}.adapterremoval",
            prepare_refgenome = "{tmpdir}/refgenome/prepare_refgenome.done",
            reads = ["{tmpdir}/seqdata/{sample}_0.fastp.fastq.gz"],
            idx = "{tmpdir}/refgenome/refgenome.amb"
        output:
            "{tmpdir}/seqdata/{sample}_MTonly_0.aligned.sam"
        message: "Align with BWA mem: {wildcards.sample}"
        resources:
            mem_gb = 8
        params:
                extra=r"-R '@RG\tID:{sample}\tSM:{sample}'"
        threads: 8
        wrapper:
            "v1.3.2/bio/bwa/mem"

    rule bwa_mem_initial_se_bam:
        input:
            "{tmpdir}/seqdata/{sample}_MTonly_0.aligned.sam"
        output:
            temp("{tmpdir}/seqdata/{sample}_MTonly_0.aligned.bam")
        message: "Convert BWA mem SE output to BAM format: {wildcards.sample}"
        params:
            extra = "-Su"
        wrapper:
            "v1.3.2/bio/samtools/view"

    rule bwa_mem_initial_merge:
        input:
            pe = "{tmpdir}/seqdata/{sample}_MTonly_12.aligned.bam",
            se = "{tmpdir}/seqdata/{sample}_MTonly_0.aligned.bam"
        output:
            "{tmpdir}/seqdata/{sample}_MTonly.aligned.bam"
        message: "Merge the aligned paired and single reads: {wildcards.sample}"
        resources:
            mem_gb = 4
        params:
            extra = "-nu"
        wrapper:
            "v1.3.2/bio/samtools/merge"
