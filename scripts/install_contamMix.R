#!/usr/bin/env Rscript
install.packages(c("coda", "getopt"),
                 repos = "https://cloud.r-project.org",
                lib = .libPaths()[length(.libPaths())])
install.packages(snakemake@params[["tarball"]],
                 type = "source",
                 lib = .libPaths()[length(.libPaths())])
