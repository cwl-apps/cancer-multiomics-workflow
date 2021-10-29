cwlVersion: v1.2
class: CommandLineTool
label: Process Phospho Protein
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
      process_phospho_protein <- function(phospho_protein_path) {
        
      #reads in phospho protein abundances and removes extra columns - write new csv
      phospho_protein = read.delim(phospho_protein_path) %>% 
        select(-Peptide, -Organism) %>% 
        pivot_longer(c(-Gene, -Phosphosite)) %>% 
        filter(!(name %in% c("X263d3f.I.Log.Ratio", "blcdb9.I.Log.Ratio", "c4155b.C.Log.Ratio"))) %>% 
        pivot_wider(values_from = value, names_from = name) %>% 
        glimpse()

      }
  - entryname: main.sh
    writable: false
    entry: |-
      #!/bin/bash
      Rscript -e 'source(file("stdin"),echo=TRUE,spaced=TRUE)' <<RSCRIPT
      # *** Add Rscript below ***
      # functionalized R scripts should be pasted into functions.R
      source("functions.R") 
      library(tidyverse)
      process_phospho_protein(
          phospho_protein_path = "$(inputs.phospho_protein.path)") %>%
          glimpse() %>%
        write_csv(paste0("$(inputs.phospho_protein.nameroot)", "_phospho_protein_log_share_ratio.csv"))
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
- id: phospho_protein
  type: File
  sbg:fileTypes: TSV

outputs:
- id: processed_phospho_protein_file
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
id: david.roberson/build-mogsa/process-phospho-protein/5
sbg:appVersion:
- v1.2
sbg:content_hash: a65392982ca17007a9a3c2a785f13db6f5275600521084db41395c6e1c2eb2106
sbg:contributors:
- david.roberson
sbg:createdBy: david.roberson
sbg:createdOn: 1634314067
sbg:id: david.roberson/build-mogsa/process-phospho-protein/5
sbg:image_url:
sbg:latestRevision: 5
sbg:modifiedBy: david.roberson
sbg:modifiedOn: 1634315492
sbg:project: david.roberson/build-mogsa
sbg:projectName: BUILD MOGSA
sbg:publisher: sbg
sbg:revision: 5
sbg:revisionNotes: ''
sbg:revisionsInfo:
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634314067
  sbg:revision: 0
  sbg:revisionNotes: Copy of david.roberson/build-mogsa/tool-template-rscript/1
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634314442
  sbg:revision: 1
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634314514
  sbg:revision: 2
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634314716
  sbg:revision: 3
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634314995
  sbg:revision: 4
  sbg:revisionNotes: ''
- sbg:modifiedBy: david.roberson
  sbg:modifiedOn: 1634315492
  sbg:revision: 5
  sbg:revisionNotes: ''
sbg:sbgMaintained: false
sbg:validationErrors: []
