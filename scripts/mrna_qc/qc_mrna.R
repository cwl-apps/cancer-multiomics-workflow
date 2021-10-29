qc_mrna <- function(mRNA_path, is_raw) {
  
  #makes histogram of RNA_Log2_Normalized values
  
  #makes Venn diagram inputs.
  
  library(tidyverse) 
  library(edgeR)
  
  mRNA_quant = read_rds(mRNA_path) %>%
    glimpse()
  
  # make histrogram of RNA_Log2_Normailized_quantile
  
  if (is_raw){
    dat            <- log2(RNA.in.f.expressed + 1);
    write.table(dat, paste0(outputFileFolder, "RNA_logTransformed.txt"), sep="\t", quote=F, row.names=T);
    dat.n <- normalizeBetweenArrays(dat, method="quantile");
    write.table(dat.n, paste0( outputFileFolder,"RNA_Log2_Normalized_quantile.txt"), sep="\t", quote=F, row.names=T);
    pdf(paste0( outputFileFolder,"histogram_of_RNA_Log2_Normalized_quantile.pdf"));
    hist(as.matrix(dat.n),main="RNA normalized counts",xlab=paste ("N(genes):", dim(dat.n)[1],"; N(samples):", dim(dat.n)[2]),breaks=100)
    dev.off();
    dat.rna <- dat.n
  }else{
    dat.l          <- log2(RNA.in.f.expressed + 1);
    write.table(dat.l, paste0(outputFileFolder, "RNA_logTransformed.txt"), sep="\t", quote=F, row.names=T);
    pdf(paste0( outputFileFolder, "histogram of ", "RNA_Log2_Transformed.pdf"));
    hist(as.matrix(dat.l),main="Log2 RNA normalized counts",xlab=paste ("N(genes):", dim(dat.l)[1],"; N(samples):", dim(dat.l)[2]) ,breaks=100)
    dev.off();
    dat.rna <- dat.l
  }
  
  #make venn object outputs
  
  #make MOGSA compatible output
  
  
  
}