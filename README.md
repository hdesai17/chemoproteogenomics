# Chemoproteogenomics 
This resource contains (1) scripts for generating sample-specifc variant contiaining FASTA databases for use in searching mass-spectromtery based proteomics data and (2) MSFragger command line pipelines for 2-stage searches. \
The [FragPipe GUI](https://github.com/Nesvilab/FragPipe) is now compatible with 2-stage searches. 
Instructions on running are located in bioXiv publication \
[Multi-omic stratification of the missense variant cysteinome](https://doi.org/10.1101/2023.08.12.553095) supplementary information.

_The updated GUI is recommended over command-line scripts._

## 1. Custom Database Generation

Generate sample-matched peptide variant-containing databases with both simple Uniprot ID FASTA headers or detailed headers.

### Before Running

 1. Move contents of CustomDB_Generation folder to working directory containing VCF files
    `cp -r CustomDB_Generation/* /path/to/your/working/directory/`
 3. Download Genocode v28 protein coding translations and GTF annotation files as well as common SNPs missense changes here \
    and move to Annotations directory
    `mv *gencode /path/to/your/working/directory/Annotations/` `mv *common /path/to/your/working/directory/Annotations/`

### Running

`sh GenerateDB.sh`

## 2. MSFragger command-line 2-stage search

Process .raw MS files with an MSFragger pipeline using Philospher and Peptide Prophet for post-processing with optional IonQuant quantitation. _The updated GUI is recommended over command-line scripts._

### Running

`sh 2-stage-run.sh`
 
_Note: Several files require path updates (see individual helper scripts)_

