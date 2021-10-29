cwlVersion: v1.2
class: CommandLineTool
label: process global protein
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: DockerRequirement
  dockerPull: cgc-images.sbgenomics.com/david.roberson/mogsa:21.06.16
- class: InitialWorkDirRequirement
  listing:
  - entryname: functions.R
    writable: false
    entry: |-
      process_global_protein <- function(global_protein_path) {
        #Reads in global protein abundances and removes extra columns
      
        # Takes one argument: path of global protein
        
        global_protein = read.delim(global_protein_path, sep = "\t") %>% 
          filter(!(Gene %in% c("Mean", "Median", "StdDev"))) %>% 
          select(-NCBIGeneID, -Authority, -Description, -Organism, -Chromosome, -Locus) %>% 
          pivot_longer(-Gene) %>% 
          rename(sample = name, log_value = value) %>% 
          filter(!(stringr::str_detect(.$sample, "X263d3f"))) %>% 
          filter(!(stringr::str_detect(.$sample, "blcdb9"))) %>%
          filter(!(stringr::str_detect(.$sample, "c4155b"))) %>%
          filter(!(stringr::str_detect(.$sample, "Unshared"))) %>% 
          pivot_wider(values_from = log_value, names_from = sample) %>% 
          glimpse()
        }
  - entryname: main.sh
    writable: false
    entry: |-
      #!/bin/bash
      Rscript -e 'source(file("stdin"),echo=TRUE,spaced=TRUE)' <<RSCRIPT
      # *** Add Rscript below ***
      # fucntionalized R scripts should be pasted into functions.R
      source("functions.R") 
      library(tidyverse)
      
      global_protein = process_global_protein(
          global_protein_path = "$(inputs.global_protein.path)"
          ) %>% 
        write_csv(paste0("$(inputs.global_protein.nameroot)", "_global_protein_log_share_ratio.csv"))
        # *** Add Rscript above ***  
      RSCRIPT
      # *** Add any additional shell scripting here.
- class: InlineJavascriptRequirement
  expressionLib:
  - |2-

    var setMetadata = function(file, metadata) {
        if (!('metadata' in file)) {
            file['metadata'] = {}
        }
        for (var key in metadata) {
            file['metadata'][key] = metadata[key];
        }
        return file
    };
    var inheritMetadata = function(o1, o2) {
        var commonMetadata = {};
        if (!o2) {
            return o1;
        };
        if (!Array.isArray(o2)) {
            o2 = [o2]
        }
        for (var i = 0; i < o2.length; i++) {
            var example = o2[i]['metadata'];
            for (var key in example) {
                if (i == 0)
                    commonMetadata[key] = example[key];
                else {
                    if (!(commonMetadata[key] == example[key])) {
                        delete commonMetadata[key]
                    }
                }
            }
            for (var key in commonMetadata) {
                if (!(key in example)) {
                    delete commonMetadata[key]
                }
            }
        }
        if (!Array.isArray(o1)) {
            o1 = setMetadata(o1, commonMetadata)
            if (o1.secondaryFiles) {
                o1.secondaryFiles = inheritMetadata(o1.secondaryFiles, o2)
            }
        } else {
            for (var i = 0; i < o1.length; i++) {
                o1[i] = setMetadata(o1[i], commonMetadata)
                if (o1[i].secondaryFiles) {
                    o1[i].secondaryFiles = inheritMetadata(o1[i].secondaryFiles, o2)
                }
            }
        }
        return o1;
    };

inputs:
- id: global_protein
  type: File
  sbg:fileTypes: TSV

outputs:
- id: processed_global_protein
  type: File?
  outputBinding:
    glob: '*.csv'
    outputEval: $(inheritMetadata(self, inputs.global_protein))
stdout: standard.out

baseCommand:
- bash
- main.sh

hints:
- class: sbg:SaveLogs
  value: '*.R'
- class: sbg:SaveLogs
  value: standard.out
- class: sbg:SaveLogs
  value: '*.sh'
- class: sbg:SaveLogs
  value: '*.py'
id: david.roberson/build-mogsa/process-global-protein/10
sbg:appVersion:
- v1.2
sbg:content_hash: a093bed6fa6ae23362141ec9e7ccced5b4c43ee4360cd09aa7deed4afd369f097
sbg:contributors:
- david.roberson
sbg:createdBy: david.roberson
sbg:createdOn: 1634306926
sbg:id: david.roberson/build-mogsa/process-global-protein/10
sbg:image_url:
sbg:latestRevision: 10
sbg:modifiedBy: david.roberson
sbg:modifiedOn: 1634313473
sbg:project: david.roberson/build-mogsa
sbg:projectName: BUILD MOGSA
sbg:publisher: sbg
sbg:revision: 10
sbg:revisionNotes: ''
sbg:revisionsInfo:
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634306926
  sbg:revision: 0
  sbg:revisionNotes:
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634307512
  sbg:revision: 1
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634307535
  sbg:revision: 2
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634310746
  sbg:revision: 3
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634310922
  sbg:revision: 4
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634311383
  sbg:revision: 5
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634311435
  sbg:revision: 6
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634312227
  sbg:revision: 7
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634312457
  sbg:revision: 8
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634313120
  sbg:revision: 9
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634313473
  sbg:revision: 10
  sbg:revisionNotes: ''
sbg:sbgMaintained: false
sbg:validationErrors: []
