

global_protein_qc <- function(global_protein_path) {
  
  Prot.in = read.delim(name.files[i], sep="\t", header=T, check.names=F);
  Prot.in <- Prot.in[order ( Prot.in[,1]),]
  Prot.in <- Prot.in[!duplicated(Prot.in[1]),]
  ## rmv NA
  Prot.in <- na.omit(Prot.in)
  row.names(Prot.in ) <- Prot.in [,1]
  Prot.in <- Prot.in[ , -1];
  list.venn [[length(MOGSA.Input) + 1]] <- colnames( Prot.in)
  Prot.in <- Prot.in[, which ( colnames(Prot.in )%in% row.names(dat.sub))]
  pdf(paste0(outputFileFolder, "histogram_of_global_Protein.pdf"));
  hist(as.matrix(Prot.in ),main="global_Protein Expression",xlab=paste ("N(genes):", dim(Prot.in)[1],"; N(samples):", dim(Prot.in)[2]), breaks=100)
  dev.off();
  write.table(Prot.in, paste0( outputFileFolder,"Global_Protein_matrix_unique_GenesNo_NA.txt"), sep="\t", quote=F, row.names=T);
  MOGSA.Input[[length(MOGSA.Input) + 1]] <- Prot.in
  
  
}  