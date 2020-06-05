#!/usr/bin/env python
# Prepare FastQ files from input files
from glob import glob
import os
import re
import subprocess

if snakemake.config['seqdatatype'] == "fastq":
    fqfns = glob(f"{snakemake.config['seqdatadir']}/{snakemake.wildcards.sample}*.{snakemake.config['seqdatasuffix']}")
    fqfns = [fn for fn in fqfns if re.search(rf"{snakemake.wildcards.sample}(_[0-2])*.{snakemake.config['seqdatasuffix']}",
                                             os.path.basename(fn))]
    for fqfn in fqfns:
        if any([f"_{i}" in os.path.basename(fqfn) for i in range(3)]):
            if not os.path.islink(f"{os.path.basename(fqfn).replace('.' + snakemake.config['seqdatasuffix'], '')}"
                                      ".raw_fastq.gz"):
                os.symlink(f"{fqfn}", f"{snakemake.params.outdir}/"
                           f"{os.path.basename(fqfn).replace('.' + snakemake.config['seqdatasuffix'], '')}"
                           ".raw_fastq.gz")
        else:
            if not os.path.islink(f"{os.path.basename(fqfn).replace('.' + snakemake.config['seqdatasuffix'], '')}"
                                      "_1.raw_fastq.gz"):
                os.symlink(f"{fqfn}", f"{snakemake.params.outdir}/"
                           f"{os.path.basename(fqfn).replace('.' + snakemake.config['seqdatasuffix'], '')}"
                           "_1.raw_fastq.gz")
elif snakemake.config['seqdatatype'] == "bam":
    subprocess.run(f"bam-fixpair -q -u {snakemake.config['seqdatadir']}/{snakemake.wildcards.sample}.{snakemake.params.suffix} | "
                   f"samtools fastq -0 {snakemake.params.pe0} -1 {snakemake.output[0]} "
                    f"-2 {snakemake.params.pe2} -",
                    shell=True)
