# Chemoproteogenomics

FragPipe GUI is compatible with 2-stage search. 
Instructions on running are located in bioXiv publication \
[Multi-omic stratification of the missense variant cysteinome](https://doi.org/10.1101/2023.08.12.553095) supplementary information.

_The updated GUI is recommended over command-line scripts._

## 1. Custom Database Generation

Generate sample-matched peptide variant-containing databases with both simple Uniprot ID FASTA headers or detailed headers.

### Running

`sh customDB-run.sh`

## 2. MSFragger command-line 2-stage search

Process .raw MS files with an MSFragger pipeline using Philospher and Peptide Prophet for post-processing with optional IonQuant quantitation. _The updated GUI is recommended over command-line scripts._

### Running

`sh 2-stage-run.sh`
 
_Note: Several files require path updates (see individual helper scripts)_

