source("global_protein_qc.R")

args = commandArgs(trailingOnly=TRUE)

qlobal_protein_qc(global_protein_path = args[1])
