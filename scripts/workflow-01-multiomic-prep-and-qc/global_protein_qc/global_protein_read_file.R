process_global_protein <- function(global_protein_path) {
  #Reads in global protein abundances and removes extra columns
  
  # Takes one argument: path of global protein
  
  global_protein = read.delim(global_protein_path, sep = "\t") %>% 
    filter(!(Gene %in% c("Mean", "Median", "StdDev"))) %>% 
    select(-NCBIGeneID, -Authority, -Description, -Organism, -Chromosome, -Locus) %>% 
    pivot_longer(-Gene) %>% 
    glimpse() %>%
    rename(sample = name, log_value = value) %>% 
    filter(!(stringr::str_detect(.$sample, "X263d3f"))) %>% 
    filter(!(stringr::str_detect(.$sample, "blcdb9"))) %>%
    filter(!(stringr::str_detect(.$sample, "c4155b"))) %>%
    filter(!(stringr::str_detect(.$sample, "Unshared"))) %>% 
    pivot_wider(values_from = log_value, names_from = sample) %>% 
    glimpse()
}