#!/bin/bash
##This script is command line version of the FragPipe 2-stage FDR analysis
##Start in a folder containing all .raw files
##NOTE: line 93 contains paths to IonQuant jar files

echo "Create File lists"
path=`pwd` 
FILES="*.raw"
readlink -f *.raw | sed -e 's!.*/!!'> temp.txt 
readlink -f *.raw | sed -e 's/.raw/\//g' > temp2.txt
readlink -f *.raw > temp3.txt

sed -e 's/.raw/.pep.xml/g' temp.txt|sed 's/^/interact-/'> temp4.txt
paste -d'\0' temp2.txt temp4.txt > temp9.txt
awk '$1=$1' ORS=' ' temp9.txt > filelist_proteinprophet.txt
sed -e 's/.raw/\/psm.tsv/g' temp3.txt|sed 's/^/--psm /' > temp4.txt

pwd | xargs echo '--specdir' | tee -a temp4.txt

awk '$1=$1' ORS=' ' temp4.txt > filelist_ionquant.txt #fp18
rm temp*.txt
chmod 777 filelist*

echo "Here is list of raw files"
echo $FILES

echo "Clean and initialize working directory"
$philosopher workspace --clean
$philosopher workspace --init

echo "MSFragger search" 
java -Xmx100g -jar $MSFragger_jar $fragger_params $FILES
echo java -Xmx100g -jar $MSFragger_jar $fragger_params $FILES
echo "Loop through folders for philosopher analysis"
n=0
for i in $FILES; do
n=$((n+1))
echo "This is the file"
echo $i
echo "This is n value"
echo $n
cd $path
mkdir ${i%.raw}
chmod 777 ${i%.raw}
echo "Move pepXML outputs"
mv ${i%.raw}.pepXML ${i%.raw}
cd ${i%.raw}
echo "This is the working directory"
pwd
echo "Clean and initialize working directory"
$philosopher workspace --clean
$philosopher workspace --init
echo "Annotate DB" #must annotate database in the philosopher workspace before running peptide and protein prophet
echo "these are the variables"
echo $prefix
echo $Reference_Database
$philosopher database --annotate $Reference_Database --prefix $prefix
echo "Peptide Prophet"
$philosopher peptideprophet --decoyprobs --ppm --accmass --nonparam --expectscore --decoy $prefix --database $Reference_Database ${i%.raw}.pepXML
done
echo "This is the working directory"
echo $path
cd $path
n=0
for i in $FILES; do
  cd $path
  cd ${i%.raw}
  n=$((n+1))
  if [ $n -eq 1 ]; then
  echo "the first iteration $n"
  echo "Protein Prophet"
  $philosopher proteinprophet --maxppmdiff 2000000 --output combined $(< ${path}/filelist_proteinprophet.txt)
  mv combined* ../
  echo "philosopher Filter"
  $philosopher filter --sequential --prot 0.01 --tag $prefix --pepxml interact*pep.xml --protxml ${path}/combined.prot.xml --razor 
  mv .meta/razor.bin ${path}/.meta
  else
  echo "the second iteration $n"
  $philosopher filter --sequential --prot 0.01 --tag $prefix --pepxml interact*pep.xml --protxml ${path}/combined.prot.xml --razorbin "${path}/.meta/razor.bin"
  fi
  echo "philosopher Report"
  $philosopher report
  $philosopher workspace --clean --nocheck
done

cd ..
$philosopher workspace --clean
$philosopher workspace --init

echo "Run Ion Quant" 
if [ $ionquant = "TRUE" ]; then

java -Xmx90G -Dlibs.bruker.dir=/path/to/Tools/FragPipe-18.0/fragpipe/tools/MSFragger-3.5/ext/bruker -Dlibs.thermo.dir=/path/to/Tools/FragPipe-18.0/fragpipe/tools/MSFragger-3.5/ext/thermo -cp /path/to/Tools/FragPipe-18.0/fragpipe/tools/jfreechart-1.5.3.jar:/home/user/Tools/FragPipe-18.0/fragpipe/tools/batmass-io-1.25.5.jar:/path/to/Tools/FragPipe-18.0/fragpipe/tools/IonQuant-1.8.0.jar ionquant.IonQuant --threads 40 --ionmobility 0 --minexps 1 --mbr 0 --maxlfq 0 --requantify 1 --mztol 10 --imtol 0.05 --rttol 0.4 --mbrmincorr 0 --mbrrttol 1 --mbrimtol 0.05 --mbrtoprun 10 --ionfdr 0.01 --proteinfdr 1 --peptidefdr 1 --normalization 0 --minisotopes 1 --minscans 1 --writeindex 0 --light $light --heavy $heavy --tp 3 --minfreq 0.5 --minions 1 --locprob 0.75 --uniqueness 0 --multidir . $(< filelist_ionquant.txt)
else echo "IonQunat OFF"
fi
echo "Done"


