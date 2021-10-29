source("library.R")

qc_mrna(mrna_quant_path = "all_genes_no_intersect.rds", 
                   data_subtype_path = "mogsa_meta.tsv",
                   is_normalized = TRUE)

