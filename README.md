# Chemoproteogenomics 
### Table of Contents: 

- [Generating sample-specifc variant contiaining FASTA databases](https://github.com/hdesai17/chemoproteogenomics#custom-database-generation) 
+ [MSFragger command line pipelines for 2-stage searches](https://github.com/hdesai17/chemoproteogenomics#msfragger-command-line-2-stage-search)

:exclamation: _The updated GUI is recommended over command-line scripts provided here._ 

The [FragPipe GUI](https://github.com/Nesvilab/FragPipe) is now compatible with 2-stage searches. Instructions on running are located in bioXiv manuscript: [Multi-omic stratification of the missense variant cysteinome](https://doi.org/10.1101/2023.08.12.553095) in supplementary information.

## Custom Database Generation

Generate sample-matched peptide variant-containing databases from VCFs. Both minimal Uniprot ID FASTA headers or detailed headers can be used in searching mass-spectromtery based proteomics data

#### Before Running
 1. Download or clone the repo

    `git clone https://github.com/hdesai17/chemoproteogenomics.git`
   
 2. Move VCF file into root directory (/chemoproteogenomics) or make sure the working directory contains VCF, Annotations/Tools folders and GenerateBD.sh script
    
 3. Download Genocode v28 protein coding translations and GTF annotation files as well as common SNPs missense changes [from this link](https://drive.google.com/drive/folders/1w1EaQC7q5uVudEMCGo-zREVJhK-YOC13?usp=sharing) and move to Annotations directory 
    
    `mv *gencode /path/to/working/directory/Annotations/` \
    `mv *common /path/to/working/directory/Annotations/`

#### Running

`./GenerateDB.sh` or `sh GenerateDB.sh`

:warning: Scripts depend on several R packages including _VariantAnnotation, GenomicFeatures, Biostrings, BSgenome.Hsapiens.UCSC.hg38, stringr, svMisc, and pbapply._

#### Outputs

In the Custom_Databases folder, there are variations of FASTA databases:
_- 2TS = two tryptic sites flanking variant sites; otherwise, they are whole protein sequences_
_- simple = only Uniprot ID (minimal) headers_
_- rev = contains reverse sequences specified as REV_
_- dedup = redundant peptide sequences are removed, regardless of transcript ID_

:warning: For minimal FASTA headers, additional post-processing is required to obtain variant IDs after using in FragPipe searches.

## MSFragger command-line 2-stage search
:exclamation: _The updated GUI is recommended over command-line scripts._ 

Process .raw MS files with an MSFragger pipeline using Philospher and Peptide Prophet for post-processing with optional IonQuant quantitation. 

#### Before Running

1. Download or clone the repo

  `git clone https://github.com/hdesai17/chemoproteogenomics.git`
  
:warning: Several files require path updates (see individual helper scripts)

#### Running
   
`./2-stage-run.sh` or `sh 2-stage-run.sh`
 


