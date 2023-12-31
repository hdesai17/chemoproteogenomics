#!/bin/bash

############################Before running#####################################
#Change ionquant to TRUE if quant workflow line 27
#Change H/L masses here 26-39
#Change path to reference database line 42
#Change path to variant database line 73
#Change filepaths to MSFragger jar, philosopher file, and fragger.params lines 18-20
#Change file path to Pipeline.sh lines 49 & 82
#Change file path to get_WT_scans.R 61
#Start in directory containing raw files
###############################################################################

echo "Set variables"
path=`pwd` #directory with raw files
name=basename `pwd`

MSFragger_jar=/path/to/Tools/FragPipe-18.0/fragpipe/tools/MSFragger-3.5/MSFragger-3.5.jar
philosopher=/path/to/Tools/FragPipe-18.0/fragpipe/tools/philosopher_v4.2.2_linux_amd64/philosopher
fragger_params=/path/to/fragger-Biotin_OTOT.params

sed -i 's/excluded_scan_list_file=.*//' "$fragger_params"

#IonQuant parameters:

#SET H/L masses for quant
ionquant=FALSE

#H/L IPIAA
#light=C463.2366
#heavy=C467.2529

#H/L Biotin
light=C463.2366
heavy=C469.2742

#H/L TEV tags
#light=C521.30741 
#heavy=C527.32122

#First pass search with WT database,edit fragger params to contain WT reference
sed -i 's/database_name =.*/database_name = \/path\/to\/_ip2_ip2_data_kbackus_database__Uniprot_Human_18432CCDS_canonical_contaminant_01-01-2020_reversed.fasta/' "$fragger_params"
sed -i 's/decoy_prefix =.*/decoy_prefix = Reverse_/' "$fragger_params"

#Match variables to fragger params file
Reference_Database=$(grep -oP "database_name = \K.*" $fragger_params)
prefix=$(grep -oP "decoy_prefix = \K.*" $fragger_params)

source /path/to/2-stage-pipeline.sh  

##Get WT scan list#######
echo "make output dir"
pwd
mkdir outputs
FILES="*.raw"
for i in $FILES; do
mv ${i%.raw} outputs
done
cd outputs

Rscript /path/to/get_WT_scans.R
mv * ../
########################


cd ../..
mkdir ${path}_variant 
mv ${path}/*.raw ${path}_variant
cd ${path}_variant


#2nd pass search, edit fragger.params to variant database parameters
sed -i 's/database_name =.*/database_name = \/path\/to\/Custom_Databases\/Final\/simple_H358-rm15_upper.fa/' "$fragger_params"
sed -i 's/decoy_prefix =.*/decoy_prefix = REV_/' "$fragger_params"
echo excluded_scan_list_file=${path}/scans.txt | cat - "$fragger_params" > temp && mv temp "$fragger_params"


#Match variables fragger params file
Reference_Database=$(grep -oP "database_name = \K.*" $fragger_params)
prefix=$(grep -oP "decoy_prefix = \K.*" $fragger_params)

source /path/to/2-stage-pipeline.sh  

sed -i 's/excluded_scan_list_file=.*//' "$fragger_params"

