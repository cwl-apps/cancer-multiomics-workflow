source("qc_mrna.R")

args = commandArgs(trailingOnly=TRUE)

qc_mrna(mRNA_path = args[1])


