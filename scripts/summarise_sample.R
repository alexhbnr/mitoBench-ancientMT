# Load libaries
library(ape)
library(data.table)
library(tidyverse)

# Number of reads
nreads_MT <- tibble(sample = snakemake@wildcards[["sample"]],
                    `number of reads` = as.numeric(unlist(str_split(readLines(snakemake@input[["flagstat"]])[1], " "))[1]))
# Sequencing depth
seqdepth <- fread(snakemake@input[["seqdepth"]],
                  col.names = c("chr", "pos", "sample")) %>%
            as_tibble() %>%
            summarise(`mean coverage` = round(mean(sample), 1),
                      `std coverage` = round(sd(sample), 0),
                      `sites >= 5-fold coverage` = sum(sample >= 5),
                      `minimum coverage` = min(sample),
                      `maximum coverage` = max(sample)) %>%
            mutate(sample = snakemake@wildcards[["sample"]]) %>%
            select(sample, everything())
# Read length
readlength <- fread(snakemake@params[["readlengthdist"]]) %>%
              as_tibble() %>%
              group_by(Length) %>%
              summarise(n = sum(Occurrences)) %>%
              mutate(sumRL = Length * n) %>%
              arrange(desc(n), desc(sumRL)) %>%
              summarise(`mode of read length` = Length[1]) %>%
              mutate(sample = snakemake@wildcards[["sample"]])
# MixEMT
mixemt <- fread(snakemake@input[["mixemt"]],
                col.names = c("comp", "hg", "contr", "noReads")) %>%
          mutate(contr = str_c(contr * 100, "%"),
                 hg_summary = str_c(hg, " (", contr, ",", noReads,")")) %>%
          summarise(mixEMT = str_replace_all(toString(hg_summary), ", ", "; ")) %>%
          mutate(sample = snakemake@wildcards[["sample"]])
# contamMix
contamMix <- fread(snakemake@input[["contammix"]],
                   select = c(1:4),
                   col.names = c("error", "propAuth",
                                 "propAuth_lowerQuantile",
                                 "propAuth_upperQuantile")) %>%
             mutate(`proportion of contamination DNA (contamMix)` = str_c(round((1 - propAuth) * 100, 2), "% (",
                                                                          round((1 - propAuth_upperQuantile) * 100, 2), "-",
                                                                          round((1 - propAuth_lowerQuantile) * 100, 2), "%)"),
                    `flagged due to contamination (> 10%)` = 1 - propAuth >= 0.1,
                    sample = snakemake@wildcards[["sample"]]) %>%
             select(sample,
                    `proportion of contamination DNA (contamMix)`,
                    `flagged due to contamination (> 10%)`)
# snpAD
fas_snpAD <- read.dna(snakemake@input[["snpAD"]],
                      "fasta", as.character = T)
consensus_comparison <- tibble(noNs = sum(fas_snpAD[1,] == "n"),
                               `% of Ns in snpAD consensus sequence` = round(noNs * 100 / 16569, 1),
                               `flagged due to number of missing bases (> 1%)` = `% of Ns in snpAD consensus sequence` > 1) %>%
                        select(-noNs) %>%
                        mutate(sample = snakemake@wildcards[["sample"]])
#HaploGrep2
## Haplogroup assignment
count_no_polys <- function(s) {
  # Count the number of polymorphisms that are not caused by missing data
  if (is.na(s)) {
    rawstring <- unlist(str_split(s, " "))
    fltstring <- rawstring[!str_detect(rawstring, "N$")]
    return(length(fltstring))
  } else return(0)
}
haplogrep <- fread(snakemake@input[["haplogrep"]],
                   colClasses = list("character"=c("SampleID")), nrows =1) %>%
             mutate_at(vars(Not_Found_Polys, Remaining_Polys), count_no_polys) %>%
             select(haplogroup = Haplogroup,
                    `haplogroup quality` = Quality,
                    `no. of not-found polymorphisms` = Not_Found_Polys,
                    `no. of remaining polymorphisms` = Remaining_Polys) %>%
             mutate(sample = snakemake@wildcards[["sample"]])
# Summary
summary_table <- nreads_MT %>%
                 left_join(seqdepth, by = "sample") %>%
                 left_join(readlength, by = "sample") %>%
                 left_join(consensus_comparison, by = "sample") %>%
                 left_join(haplogrep, by = "sample") %>%
                 left_join(contamMix, by = "sample") %>%
                 left_join(mixemt, by = "sample") %>%
                 mutate(quality = ifelse(`flagged due to number of missing bases (> 1%)` | `flagged due to contamination (> 10%)`,
                                         "low", "high"))
fwrite(summary_table,
       sep = "\t",
       file = snakemake@output[[1]])
