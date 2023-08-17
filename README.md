# Chemoproteogenomics 
This resource contains scripts for generating sample-specifc variant contiaining FASTA databases for use in searching mass-spectromtery based proteomics data and MSFragger command line pipelines for 2-stage searches. 

The [FragPipe GUI](https://github.com/Nesvilab/FragPipe) is now compatible with 2-stage searches. \
Instructions on running are located in bioXiv manuscript: \
[Multi-omic stratification of the missense variant cysteinome](https://doi.org/10.1101/2023.08.12.553095) in supplementary information.

:exclamation: _The updated GUI is recommended over command-line scripts provided here._

## Custom Database Generation

Generate sample-matched peptide variant-containing databases with both minimal Uniprot ID FASTA headers or detailed headers.

### Before Running

 1. Move contents of CustomDB_Generation folder to working directory containing VCF files
    
    `cp -r CustomDB_Generation/* /path/to/your/working/directory/`
    
 3. Download Genocode v28 protein coding translations and GTF annotation files as well as common SNPs missense changes \
    [from this link](https://drive.google.com/drive/folders/1w1EaQC7q5uVudEMCGo-zREVJhK-YOC13?usp=sharing) and move to Annotations directory 
    
    `mv *gencode /path/to/your/working/directory/Annotations/` \
    `mv *common /path/to/your/working/directory/Annotations/`

### Running

`./GenerateDB.sh` or `sh GenerateDB.sh`

:warning: Scripts require several R packages including _VariantAnnotation, GenomicFeatures, Biostrings, BSgenome.Hsapiens.UCSC.hg38, stringr, svMisc, and pbapply._

For minimal FASTA headers, additional post-processing is required to obtain variant IDs after using in FragPipe searches.

## MSFragger command-line 2-stage search
:exclamation: _The updated GUI is recommended over command-line scripts._ 

Process .raw MS files with an MSFragger pipeline using Philospher and Peptide Prophet for post-processing with optional IonQuant quantitation. 

### Running

`./2-stage-run.sh` or `sh 2-stage-run.sh`
 
:warning: Several files require path updates (see individual helper scripts)

