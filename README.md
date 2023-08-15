# Chemoproteogenomics

FragPipe GUI is compatible with 2-stage search 

link paper

_The updated GUI is recommended over command-line scripts._



## 1. Custom Database Generation

Generate sample-matched peptide custom databases with both simple Uniprot ID FASTA headers or detailed headers.

### Running

`sh customDB-run.sh`

## 2. MSFragger command-line 2-stage search

_The updated GUI is recommended over command-line scripts._

Process .raw MS files with an MSFragger pipeline using Philospher and Peptide Prophet for post-processing with optional IonQuant quantitation.

### Running

`sh 2-stage-run.sh`
 
_Note: Several files require path updates (see individual helper scripts)_

