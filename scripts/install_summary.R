#!/usr/bin/env Rscript
install.packages(c("ape", "data.table", "tidyverse"),
                 repos = "https://cloud.r-project.org",
                lib = .libPaths()[length(.libPaths())])
