################################################################################
# Summary report for the MT pipeline
#
# Snakefile to create a summary report CSV that contains the typical columns
# Mark has usually asked for.
#
# Alex Huebner, 23/01/18
################################################################################

from snakemake.utils import R

workdir: config['projdir']

rule summary:
    input: "summary.csv",
           "summary.RData"

rule align_FastA:
    output: "consensus.fasta.gz"
    message: "Create multi-sequence alignment of all FastA sequences"
    shell:
        """ 
        mafft <(cat /mnt/genotyping/sk_pipelines/source/RSRS.fasta fasta/*.fasta) | gzip > {output}
        """ 

rule summary_report:
    input:  fas="consensus.fasta.gz"
    output: csv="summary.csv",
            RData="summary.RData"
    message: "Summarise the information into a table"
    run:
        R("""
        library(data.table)
        library(tidyverse)
        library(stringr)
        library(ape)

        # Processing of FastA files
        consensus_fas.fn <- gzfile("consensus.fasta.gz", open="r")
        consensus_fas <- read.dna(consensus_fas.fn, "fasta", as.character=T)
        rsrs_gaps <- table(consensus_fas[1,])[1]
        fastas <- do.call("rbind", apply(consensus_fas[2:nrow(consensus_fas),], 1, function(r) {{
                        alleles <- table(factor(r, levels=c("a", "c", "g", "t", "n", "-")))
                        data_frame(`consensus length [bp]` = sum(alleles[1:5]),
                                `indel [bp]` = (-1) * (alleles[6] - rsrs_gaps),
                                `N [bp]` = alleles[5])
                }})) %>%
                mutate(sample = dimnames(consensus_fas)[[1]][2:nrow(consensus_fas)])
        # Coverage
        print("Summarise the coverage")
        sequencing_depth <- fread("sequencing_depth.txt", header=T) %>%
                            gather("sample", "coverage", 3:ncol(.))
        sequencing_depth_summary <- sequencing_depth %>%
                                    group_by(sample) %>%
                                    summarise(`mean coverage` = mean(coverage),
                                            `min. coverage` = min(coverage),
                                            `max. coverage` = max(coverage))
        # Haplogroup
        print("Parse haplogroups")
        haplogrep <- fread("haplogroup_calls.hsd", select=c(1, 3, 4)) %>%
                     rename(haplogroup = Haplogroup)
        # Heteroplasmy
        print("Analyse minor alleles")
        heteroplasmy.fns <- list.files("heteroplasmy", pattern="heteroplasmy\\\.txt", full.names=T)
        heteroplasmy <- do.call("rbind", lapply(heteroplasmy.fns, function(fn) {{
                                sampleid <- str_match(basename(fn), "([A-Za-z0-9]+).heteroplasmy.txt")[,2]
                                print(sampleid)
                                hetero <- fread(fn) %>%
                                          left_join(filter(sequencing_depth, sample == sampleid) %>%
                                                    select(-chr, -sample),
                                                        by=c("position"="pos")) %>%
                                          filter(coverage > 500)
                                if (nrow(hetero) > 10) {{
                                    data_frame(sample = sampleid,
                                               `positions w/ minor alleles` = paste(nrow(hetero), "positions"))
                                }} else if (nrow(hetero) > 0) {{
                                    data_frame(sample = sampleid,
                                            `positions w/ minor alleles` = paste("positions", paste(hetero$position, collapse=",")))
                                }} else {{
                                    data_frame(sample = sampleid,
                                            `positions w/ minor alleles` = "")
                                }}
                            }}))
        # Contamination
        print("Analyse Mingkun's contamination test results")
        conta <- lapply(readLines("haplogroup_contamination.txt")[1:(length(heteroplasmy.fns) + 1)],
                        function(l) {{
                cols <- unlist(strsplit(l, "\t"))[1:2]                    
                if (length(cols) == 1) data_frame(sample = cols[1], `no. of HG defined by >= 5 heteroplasmies` = 0)
                else data_frame(sample = cols[1], `no. HG defined by >= 5 heteroplasmies` = cols[2])
                }})
        conta <- do.call("rbind", conta[2:length(conta)])  %>%
                 mutate(sample = factor(sample)) %>%
                 filter(`no. HG defined by >= 5 heteroplasmies` != "-") %>%
                 mutate(n = str_count(`no. HG defined by >= 5 heteroplasmies`, ",") + 1) %>%
                 select(sample, n) %>%
                 complete(sample, fill=list(n=0)) %>%
                 left_join(select(sequencing_depth_summary, sample, `mean coverage`), by="sample" ) %>%
                 mutate(n = ifelse(`mean coverage` >= 500, n, "not enough power")) %>%
                 select(-`mean coverage`) %>%
                 rename(`no. HG defined by >= 5 positions w/ minor allelies` = n)
        # MixEMT
        print("Analyse mixEMT results")
        mixemt_fns <- list.files("mixemt", pattern="mixemt\\\.log", full.names=T)
        mixemt <- do.call("rbind", lapply(mixemt_fns, function(fn) {{
                    if (file.info(fn)$size > 0) {{
                      mixemt_contribution <- fread(fn, header=F,
                                                   col.names=c("hgid", "hg", "proportion", "nReads")) %>%
                                             mutate(percent = round(proportion * 100, 1),
                                                    label = paste0(hg, ":", percent, "%")) %>%
                                             pull(label) %>%
                                             paste(., collapse=";")
                      data_frame(sample = str_replace(basename(fn), "\\\.mixemt\\\.log", ""),
                                 `mixEMT composition` = mixemt_contribution)
                    }} else {{
                      data_frame(sample = str_replace(basename(fn), "\\\.mixemt\\\.log", ""),
                                 `mixEMT composition` = "no inference possible")
                    }}
                  }}))
        # Join
        overview <- left_join(fastas, sequencing_depth_summary, by="sample") %>%
                    left_join(., haplogrep, by=c("sample"="SampleID")) %>%
                    left_join(., heteroplasmy, by="sample") %>%
                    left_join(., conta, by="sample")  %>%
                    left_join(., mixemt, by="sample")
        cat("haplogroup: haplogroup inferred by HaploGrep v2\n",
            file="{output.csv}")
        cat("indel: sum of insertions and deletions compared to RSRS (positive value: excess of insertions; negative value: excess of deletions)\n",
            file="{output.csv}", append=T)
        cat("positions w/ minor allele: number of positions, which carry a minor allele with a frequency of > 2%\n",
            file="{output.csv}", append=T)
        cat("no. HG defined by >= 5 positions w/ minor allelies: contamination indicator for contamination of an individual with a different haplogroup\n",
            file="{output.csv}", append=T)
        cat("mixEMT composition: contributing haplogroups identified from mixEMT\n\n",
            file="{output.csv}", append=T)
        write.table(overview[, c(4, 8, 5, 1, 3, 10, 11, 12)], file="{output.csv}", sep="\t",
                    row.names=F, quote=F, append=T)
        save(overview, file="{output.RData}")
        """)

