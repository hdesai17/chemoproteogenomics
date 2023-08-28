##Edited to include multi-variant SAAVs (REFAA and VARAA must be same length)

all_missense<-all_missense[-grep("\\*",all_missense$VARAA)]
all_missense<-all_missense[-grep("\\*",all_missense$REFAA)]

if (combos){all_missense<-all_missense[nchar(all_missense$REFAA) == 1 & nchar(all_missense$VARAA) == 1]}


temppep=pep[which(names(pep)=="23")]

for (i in 1:length(all_missense)){
  print("Iteration")
  print(i)
  temp=pep[which(names(pep)%in%as.character(all_missense$TXID[i]))]
  proloc=all_missense$PROTEINLOC[i][[1]]
  if (length(temp) == 0){
    print("not in the list")
    next
  }
  if (length(temp[[1]]) >= as.numeric(proloc)[1]){
    if (length(proloc) > 1){
      seq=seq(as.numeric(proloc[[1]]),as.numeric(proloc[[2]]))
    }else{seq=proloc}
    ref=as.character(temp[[1]][seq])
    ref1=temp[[1]][seq]
    ref2=as.character(all_missense$REFAA[i][[1]])
    alt=as.character(all_missense$VARAA[i][[1]])
    alt1=all_missense$VARAA[i][[1]]
    
    if (ref == ref2 & nchar(ref) ==nchar(alt)){
      pos=c(mapply(function(x, y) which(x != y), strsplit(ref, ""), strsplit(alt, ""))) 
      temp[[1]][seq][pos]<-as(tolower(as.character(alt1[pos])),"AAString")
      names(temp)<-paste(names(temp),all_missense$paramRangeID[i],paste(ref2,paste0(proloc,collapse = ":"),alt,sep=""), sep = " ")
      temppep<-c(temppep,temp)
    }else{
      print("no match")
      print(i)
    }
  }else{
    print("not the same length")
    print(i)
    print(all_missense$PROTEINLOC[i])
  }
} 


rna_exome<-paste(str_split_fixed(names(temppep),"[_]+",2)[,1], str_split_fixed(names(temppep),"[ ]+",3)[,3], sep = " ")    

for (i in 1:length(temppep)){
  
  dups=temppep[grep(rna_exome[i], rna_exome)]
  if (length(dups)>1 & nchar(as.character(names(dups)[1])) < 40){
    print(i)
    names(temppep)[grep(rna_exome[i], rna_exome)]<-paste(names(dups),collapse = " and ")
  }
}


temppep_2<- temppep[!duplicated(names(temppep))] #removed the extra pairs of similar names from above loop

all_same_txids <- c(temppep_2[duplicated(str_split_fixed(names(temppep_2),"[_]+",2)[,1])],temppep_2[duplicated(str_split_fixed(names(temppep_2),"[_]+",2)[,1], fromLast = TRUE)])


print("Finished")
#save(temppep, file="temppep.RData")
system("mkdir Temp")
save(temppep_2, file="Temp/temppep_2.RData")
#save(rna_exome, file="rna_exome.RData")
save(all_same_txids, file="Temp/all_same_txids.RData")
print("Saved")
