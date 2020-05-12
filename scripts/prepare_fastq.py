#!/usr/bin/env python
# Prepare FastQ files from input files
from glob import glob
import os
import subprocess

if snakemake.config['seqdatatype'] == "fastq":
    for fqfn in glob(f"{snakemake.config['seqdatadir']}/{snakemake.wildcards.sample}*.{snakemake.config['seqdatasuffix']}"):
        if any([f"_{i}" in os.path.basename(fqfn) for i in range(3)]):
            os.symlink(f"{fqfn}", f"{snakemake.params.outdir}/"
                        f"{os.path.basename(fqfn).replace('.' + snakemake.config['seqdatasuffix'], '')}"
                        ".raw_fastq.gz")
        else:
            os.symlink(f"{fqfn}", f"{snakemake.params.outdir}/"
                        f"{os.path.basename(fqfn).replace('.' + snakemake.config['seqdatasuffix'], '')}" 
                        "_1.raw_fastq.gz")
elif snakemake.config['seqdatatype'] == "bam":
    subprocess.run(f"samtools fastq -0 {snakemake.params.pe0} -1 {output[0]} "
                    f"-2 {snakemake.params.pe2} {snakemake.config['seqdatadir']}/{snakemake.wildcards.sample}.{snakemake.params.suffix}",
                    shell=True)
