#!/usr/bin/env python
# Revert bam-wrap by removing multiple occurrences of reads
import pysam

if snakemake.params.state == "Pass":
    bamfile = pysam.AlignmentFile(snakemake.input[0])
    reads = []
    prevread = ""
    with pysam.AlignmentFile(snakemake.output[0], "wb", template=bamfile) as outfile:
        for i, read in enumerate(bamfile):
            if not read.is_unmapped and not read.mate_is_unmapped:
                if read.query_name != prevread and prevread != "":
                    if ((reads[0].is_paired and len(reads) > 2) or
                        (not reads[0].is_paired and len(reads) > 1)):  # one duplicated read
                        read_orient = [r.is_read1 for r in reads]
                        if sum(read_orient) > 1:  # forward read got duplicated
                            del reads[1]
                        else:  # reverse read got duplicated
                            if len(reads) > 2:  # for paired reads
                                del reads[2]
                            else:  # for single reads
                                del reads[1]
                    for r in reads:
                        outfile.write(r)
                    reads = []
                reads.append(read)
                prevread = read.query_name
else:
    with open(snakemake.output[0], "wt") as outfile:
        pass
