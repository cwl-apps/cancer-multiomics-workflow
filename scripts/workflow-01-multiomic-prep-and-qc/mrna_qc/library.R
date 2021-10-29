qc_mrna <- function(mrna_quant_path, data_subtype_path, is_normalized) {
  
  #mrna_quant_path = "all_genes_no_intersect.rds"
  #data_subtype_path = "mogsa_meta.tsv"
  #is_normalized = TRUE

  
  library(edgeR)
  dat.sub <- read.delim  (data_subtype_path,sep="\t", check.names=F)
  row.names(dat.sub) <- dat.sub [,1]
 
    RNA.in <-  readRDS(mrna_quant_path)
        RNA.in [1:3,1:4]
        RNA.in <- RNA.in[order ( RNA.in[,1]),]
        RNA.in <- RNA.in[!duplicated(RNA.in[1]),]
        dim(RNA.in)
        
        RNA.in = na.omit(RNA.in)
        
        row.names(RNA.in ) <- RNA.in[,1]
        RNA.in <- RNA.in[,-1]
        RNA.in.f <- RNA.in[, which ( colnames(RNA.in) %in% row.names(dat.sub))]
        dim(RNA.in.f)
        RNA.in.f <- na.omit(RNA.in.f)
        dim(RNA.in.f)
        IsExpr <- rowSums(cpm(RNA.in.f) > 1) >= nrow(dat.sub)/2;
        RNA.in.f.expressed <- RNA.in.f[IsExpr, ];
        RNA.in.f.expressed <- as.matrix(RNA.in.f.expressed)
        dim(RNA.in.f.expressed)
        ###### need to look at this.
        
        if (is_normalized){
          dat            <- log2(RNA.in.f.expressed + 1);
          write.table(dat, "RNA_logTransformed.txt", sep="\t", quote=F, row.names=T);
          
          #error here
          dat.n <- normalizeBetweenArrays(dat, method="quantile");
          write.table(dat.n, paste0( outputFileFolder,"RNA_Log2_Normalized_quantile.txt"), sep="\t", quote=F, row.names=T);
          pdf(paste0( outputFileFolder,"histogram_of_RNA_Log2_Normalized_quantile.pdf"));
          hist(as.matrix(dat.n),main="RNA normalized counts",xlab=paste ("N(genes):", dim(dat.n)[1],"; N(samples):", dim(dat.n)[2]),breaks=100)
          dev.off();
          dat.rna <- dat.n
        }
        
        if(!(is_normalized)){
          print ("it goes hese--  normalized count")
          dat.l          <- log2(RNA.in.f.expressed + 1);
          write.table(dat.l, paste0(outputFileFolder, "RNA_logTransformed.txt"), sep="\t", quote=F, row.names=T);
          pdf(paste0( outputFileFolder, "histogram of ", "RNA_Log2_Transformed.pdf"));
          hist(as.matrix(dat.l),main="Log2 RNA normalized counts",xlab=paste ("N(genes):", dim(dat.l)[1],"; N(samples):", dim(dat.l)[2]) ,breaks=100)
          dev.off();
          dat.rna <- dat.l
        }
  
  }
  
  
  