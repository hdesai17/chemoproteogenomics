#Get WT scans to pipe to MSFragger search
#Loops through psms files in 'outputs' folder and outputs list of scan numbers to file scans.txt
getwd()
fnames<-list.files()
files = lapply(paste(fnames,"psm.tsv", sep= "/"), read.csv, sep = "\t") ###change to correct number of files
scans<-list()
for (i in 1:length(files)){
  files[[i]][,1]<-gsub("*\\.0*",".",files[[i]][,1])  ##this should remove leading zeros
  scans<-c(files[[i]][,1],scans)
}
scans<-unlist(scans)
write.table(scans, row.names = F,col.names = F,quote = F, file = "scans.txt")
