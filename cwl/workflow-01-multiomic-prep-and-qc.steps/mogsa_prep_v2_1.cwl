cwlVersion: v1.2
class: CommandLineTool
label: Multiple Data Type Quality Control
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: cgc-images.sbgenomics.com/david.roberson/mogsa:21.06.16
- class: InitialWorkDirRequirement
  listing:
  - entryname: functions.R
    writable: false
    entry: |
      
      
      mogsa_multi_data_type_processing <- function(mrna_quant, protein_quant, phospho_quant, cnv_quant, subtype_file) {
      ###### This is the preparation of input for MOGSA multiOmics. 
      ## output will be QC graphs and MOGSA input for next step. User should review before going to next step to make sure that what they expect.
      
      ### we are expecting user will input different combination # of datatypes. at the back end we should be able to get the description of the input files.
      ### Description of file and list of files will be in order. 
      
      descrition_files <- c("mrna_quant", "protein_quant", "phospho_quant","cnv_quant")
      
      
        
      }
      
      
      qc_mrna <- function() {
          
          
      }
      
      qc_global_protein <- function() {
          
          
      }
      
      qc_phospho_protein <- function() {
          
          
      }
      
      create_venn_diagram <- function() {
          
          
      }
      
      create_mogsa_compatible_object <- function() {
          
          
      }
      
      ## way to get the description of whether it is raw count or not.
      ## please change it TRUE if raw count is used.
      raw.count=FALSE;
      ## obtain the files Names from user.
      Index.files.in <- which ( descrition.files.user ==TRUE)
      Datatype.list <- NULL;
      name.files <- NULL;
      
      for ( i in 1:length(Index.files.in)){
        Datatype.list <- append(Datatype.list, descrition.files[Index.files.in[i]])
        name.files <- append(name.files, List_Files[Index.files.in[i]])
      }
      
      
      
      ######### RNA stuff
      library(edgeR);
      
      list.venn <- NULL;
      
      if (Datatype.list[i]== "rna_expression" ){
            RNA.in <-  read.delim(name.files[i], sep="\t", header=T, check.names=F);
            RNA.in [1:3,1:4]
            RNA.in <- RNA.in[order ( RNA.in[,1]),]
            RNA.in <- RNA.in[!duplicated(RNA.in[1]),]
            dim(RNA.in)
            row.names(RNA.in ) <- RNA.in[,1]
            RNA.in <- RNA.in[,-1]
            list.venn [[length(MOGSA.Input) + 1]] <- colnames( RNA.in)
            RNA.in.f <- RNA.in[, which ( colnames(RNA.in) %in% row.names(dat.sub))]
            dim(RNA.in.f)
            RNA.in.f <- na.omit(RNA.in.f)
            dim(RNA.in.f)
            IsExpr <- rowSums(cpm(RNA.in.f) > 1) >= nrow(dat.sub)/2;
            RNA.in.f.expressed <- RNA.in.f[IsExpr, ];
            RNA.in.f.expressed <- as.matrix(RNA.in.f.expressed)
            dim(RNA.in.f.expressed)
            ###### need to look at this.
            
            if ( raw.count){
              dat            <- log2(RNA.in.f.expressed + 1);
              write.table(dat, paste0(outputFileFolder, "RNA_logTransformed.txt"), sep="\t", quote=F, row.names=T);
              dat.n <- normalizeBetweenArrays(dat, method="quantile");
              write.table(dat.n, paste0( outputFileFolder,"RNA_Log2_Normalized_quantile.txt"), sep="\t", quote=F, row.names=T);
              pdf(paste0( outputFileFolder,"histogram_of_RNA_Log2_Normalized_quantile.pdf"));
              hist(as.matrix(dat.n),main="RNA normalized counts",xlab=paste ("N(genes):", dim(dat.n)[1],"; N(samples):", dim(dat.n)[2]),breaks=100)
              dev.off();
              dat.rna <- dat.n
            }else{
              print ("it goes hese--  normalized count")
              dat.l          <- log2(RNA.in.f.expressed + 1);
              write.table(dat.l, paste0(outputFileFolder, "RNA_logTransformed.txt"), sep="\t", quote=F, row.names=T);
              pdf(paste0( outputFileFolder, "histogram of ", "RNA_Log2_Transformed.pdf"));
              hist(as.matrix(dat.l),main="Log2 RNA normalized counts",xlab=paste ("N(genes):", dim(dat.l)[1],"; N(samples):", dim(dat.l)[2]) ,breaks=100)
              dev.off();
              dat.rna <- dat.l
            }
            ############# 
            
      
      #### Global protein code
      
      if (Datatype.list[i]== "protein_expression" ){
            print ( "go here- prot.in")
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
            
            
      #### phospoh protein code      
      if( Datatype.list[i]== "phospho_Expression" ){
            print ( "go phospho.in")
            phospho.in <- read.delim(name.files[i], sep="\t", header=T, check.names=F);
            phospho.in <- phospho.in[order ( phospho.in [,1]),]
            phospho.in<- phospho.in[!duplicated(phospho.in[1]),]
            row.names(phospho.in ) <- phospho.in [,1]
            phospho.in <- phospho.in[ , -1]
            ##
            list.venn [[length(MOGSA.Input) + 1]] <- colnames( phospho.in)
            ##############
            phospho.in <- phospho.in[, which ( colnames(phospho.in) %in% row.names(dat.sub))]
            pdf(paste0( outputFileFolder,"histogram_of_Phospho_Protein.pdf"));
            hist(as.matrix(phospho.in ),main="Phospho_Protein Expression",xlab=paste ("N(genes):", dim(phospho.in)[1],"; N(samples):", dim(phospho.in)[2]),breaks=100)
            dev.off();
            write.table(phospho.in, paste0( outputFileFolder,"Phospho_Protein_matrix_unique_GenesNo_NA.txt"), sep="\t", quote=F, row.names=T);
            MOGSA.Input[[length(MOGSA.Input) + 1]] <- phospho.in
          }else{
            cnv.in <- read.delim(name.files[i], sep="\t", header=T, check.names=F)
            cnv.in <- cnv.in[order ( cnv.in[,1]),]
            cnv.in <- cnv.in[!duplicated(cnv.in[1]),]
            cnv.in[1:2,1:3]
            row.names(cnv.in)  <- cnv.in[,1]
            cnv.in <- cnv.in[,-1]
            list.venn [[length(MOGSA.Input) + 1]] <- colnames( cnv.in)
            cnv.in <- cnv.in[, which ( colnames(cnv.in) %in% row.names(dat.sub))]
            iscnv.expressed <- which(rowSums(as.matrix(cnv.in))>0)
            cnv.in.expressed <- cnv.in[iscnv.expressed ,]
            pdf(paste0( outputFileFolder,"Histogram_of_CNA_Expression.pdf"));
            hist(as.matrix(cnv.in.expressed ),main="CNA Expression",xlab=paste ("N(genes):", dim(cnv.in.expressed)[1],"; N(samples):", dim(cnv.in.expressed)[2]),breaks=100)
            dev.off();
            write.table(phospho.in, paste0( outputFileFolder,"CNA_matrix_unique_GenesNo_NA_filtered_rowSumGreaterthanZero.txt"), sep="\t", quote=F, row.names=T);
          }
        }
        
      names( MOGSA.Input) <- Datatype.list
      names(list.venn) <-  Datatype.list
      ### venn diagram
      list.venn [[length(MOGSA.Input) + 1]] <- row.names( dat.sub)
      names(list.venn)[length(Datatype.list)+1] <- "subtypes"
        vp <- venn.diagram(list.venn, fill = 2:(length(list.venn)+1), alpha = 0.3, filename = NULL, cex= 2);
        pdf (paste0(outputFileFolder,"Venn_Diagram_datatypesWsubtypes.pdf"));
        grid.draw(vp)
        dev.off();
        sample_common <- Reduce(intersect,  list.venn)
        ########## adjust colnames the same order.
        for ( i in 1:length(MOGSA.Input )){
          MOGSA.Input[[i]] <- MOGSA.Input[[i]][,sample_common]
          print ( identical ( colnames(MOGSA.Input[[i]]), sample_common))
        }
        ##### double check inside
        for ( i in 2:length(MOGSA.Input )){
          
          print ( identical ( colnames(MOGSA.Input[[i]]), colnames(MOGSA.Input[[1]])))
        }
        write.table(dat.sub[sample_common,], paste0(outputFileFolder, "subtypes_updated_common_samples.txt"), sep="\t", quote=F, row.names=F);
        saveRDS(MOGSA.Input, file = paste0(outputFileFolder,"MOGSA_all_expression.rds"))
        ###########
        print ( "QC graphs available for review ")
      }
      
      
      
  - entryname: main.sh
    writable: false
    entry: |-
      #!/bin/bash
      #Rscript -e 'source(file("stdin"),echo=TRUE,spaced=TRUE)' <<RSCRIPT
      # *** Add Rscript below ***
      # functionalized R scripts should be pasted into functions.R
      source("functions.R") 
      library(tidyverse)
      
      ${if(inputs.mrna_quant){return 'qc_mrna("'+inputs.mrna.path+'")'}else{return "#no mRNA input"}}
      
      ${if(inputs.global_protein){return 'qc_global_protein("'+inputs.global_protein.path+'")'}}
      
      ${if(inputs.phospho_protein){return 'qc_phospho_protein("'+inputs.phospho_protein.path+'")'}}
      
      #create_venn_diagram()
      
      #create_mogsa_compatible_object()
      
      
      
        # *** Add Rscript above ***  
      RSCRIPT
      # *** Add any additional shell scripting here.
- class: InlineJavascriptRequirement

inputs:
- id: global_protein
  type: File
  inputBinding:
    position: 0
    shellQuote: false
- id: phospho_protein
  type: File
  inputBinding:
    position: 0
    shellQuote: false
- id: mrna_quant
  type: File?
  inputBinding:
    position: 0
    shellQuote: false

outputs: []

baseCommand:
- bash
- main.sh

hints:
- class: sbg:SaveLogs
  value: '*.sh'
- class: sbg:SaveLogs
  value: '*.R'
- class: sbg:SaveLogs
  value: standard.out
id: david.roberson/build-mogsa/mogsa-prep-v2-1/3
sbg:appVersion:
- v1.2
sbg:content_hash: a5742fdc57800063c0c1a3ebccb7537a0ee835e920c7bec193c8bfba1f0ee7c1e
sbg:contributors:
- david.roberson
sbg:createdBy: david.roberson
sbg:createdOn: 1634830262
sbg:id: david.roberson/build-mogsa/mogsa-prep-v2-1/3
sbg:image_url:
sbg:latestRevision: 3
sbg:modifiedBy: david.roberson
sbg:modifiedOn: 1634843063
sbg:project: david.roberson/build-mogsa
sbg:projectName: BUILD MOGSA
sbg:publisher: sbg
sbg:revision: 3
sbg:revisionNotes: ''
sbg:revisionsInfo:
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634830262
  sbg:revision: 0
  sbg:revisionNotes:
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634830305
  sbg:revision: 1
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634842588
  sbg:revision: 2
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634843063
  sbg:revision: 3
  sbg:revisionNotes: ''
sbg:sbgMaintained: false
sbg:validationErrors: []
