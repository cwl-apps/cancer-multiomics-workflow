#FROM bioconductor/bioconductor_docker:RELEASE_3_12
      
#RUN Rscript -e 'BiocManager::install("mogsa")'
#RUN Rscript -e 'install.packages("rmarkdown")'
#RUN Rscript -e 'BiocManager::install("BiocStyle")'
#RUN Rscript -e 'install.packages("magick")'
#RUN Rscript -e 'BiocManager::install("org.Hs.eg.db")'

#FROM cgc-images.sbgenomics.com/david.roberson/mogsa:1.24.0

#RUN Rscript -e 'install.packages("tidyverse")'
#RUN Rscript -e 'BiocManager::install("EnsDb.Hsapiens.v86", update = F, ask = F)'

#FROM cgc-images.sbgenomics.com/david.roberson/mogsa:210521
#RUN Rscript -e 'install.packages("VennDiagram")'
#RUN Rscript -e 'install.packages("venn")'

#FROM cgc-images.sbgenomics.com/david.roberson/mogsa:21.06.16
#RUN Rscript -e 'BiocManager::install("mixOmics", update = F, ask = F)'

FROM cgc-images.sbgenomics.com/david.roberson/mogsa:21.07.28
RUN Rscript -e 'BiocManager::install("edgeR", updated = F, ask = F)'
RUN Rscript -e 'BiocManager::install("TCGAbiolinks", update = F, ask = F)'
RUN Rscript -e 'BiocManager::install("biomaRt", update = F, ask = F)'
