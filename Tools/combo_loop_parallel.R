library(stringr)
library(foreach)
library(doParallel)
library(VariantAnnotation)
library(GenomicFeatures)
library(svMisc)
library(doSNOW)

#setup parallel backend to use many processors #pkill -fe RSOCKnode.R to kill
load("Temp/all_same_txids.RData")
load("Temp/temppep_2.RData")
cores=detectCores()
cl <- makeCluster(cores-2) 
registerDoSNOW(cl)
same_txids<-unique(all_same_txids[grep(sample_name,names(all_same_txids))]) #removes duplicate protein sequences
iterations <- length(unique(str_split_fixed(names(same_txids),"[_]+",2)[,1]))
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)



same_txids<-all_same_txids[grep(sample_name,names(all_same_txids))]  ##grep cell line
same_txids<-unique(same_txids) #removes duplicate protein sequences ignoring transcript ID
string=AAStringSet()
dups= AAStringSet() 
modified_peps= same_txids[which(names(same_txids)== 0)]
all_combo_peps <- list()
comb<-list()


source('Tools/helper_function.R')  

close_pos <-NULL
new_pos <- NULL
loop<-unique(str_split_fixed(names(same_txids),"[_]+",2)[,1])

final<-foreach(i=1:iterations, .combine=c,
               .packages = c('stringr','VariantAnnotation','GenomicFeatures','svMisc','doParallel','foreach'), .options.snow = opts) %dopar% {
  out_peps<-list()
  out_peps2<-list()
  out_peps3<-list()
  print(i)
  string=same_txids[which(str_split_fixed(names(same_txids),"[_]+",2)[,1]== str_split_fixed(loop[i],"[_]+",2)[1])] ###string has single variant peptides for each TXID
  loc_var<-NULL
  adapted_list<- Biostrings::AAStringSet()
  if(length(string)>1 & length(string) < 16){
  for (w in 1:length(string)){
    for(x in 1:(length(string[[w]]))){
      if(!grepl("^[[:upper:]]+$", string[[w]][x]) & as.character(string[[w]][x])!= "*"){ ## if the letter is not uppercase and the letter is not *
        loc_var<-c(loc_var,x)                                                             ### x is the position in the string
        x<-loc_var
        new<-list()
        for (q in 1:length(loc_var)){
        close_pos<-as.data.frame(which(abs(x-loc_var[q]) <= 30))
        new<-append(new,close_pos)
        } 
        new<-unique(new)
      }
    }
  }
  for (z in 1:length(new)){
          close_pos<-new[[z]]
          if (length(close_pos) >= 2){  
          adapted_list<-string[close_pos]
  if(length(adapted_list) > 2){
      print(paste("adapted list is length", length(adapted_list)))
      print(paste("multiple peptides combo replacement for", length(adapted_list), "peptides"))
        print(loop[i])
	
        for(j in 1:length(adapted_list)){
          if(j > 1){
            dimsofcomb=combn(length(adapted_list), j)
            ch_peps=all_multicombs_of_var_peps(list_of_peps = adapted_list, dims_of_combs = dimsofcomb)
            out_peps= c(out_peps, ch_peps) 
           }
          }
  }else if (length(adapted_list) == 2){
   print("duplicate")
    for (k in 1:length(adapted_list[[1]])){
      if (adapted_list[[1]][k]!= adapted_list[[2]][k]){
        if (tolower(adapted_list[[1]][k]) == as.character(adapted_list[[1]][k])){
          adapted_list[[2]][k] <-as.character(adapted_list[[1]][k]) 
          out_peps2<-adapted_list[2]
          names(out_peps2)<-paste(names(adapted_list[1]),names(adapted_list[2]),sep = ",")
          out_peps3<-c(out_peps2,out_peps3)
     }
    }
   }
  }else{
    print(paste("do nothing because",length(adapted_list), sep= " "))
    print(loop[i])
  }
    
   }
  } 
    all_combo_peps[[i]]= c(out_peps,out_peps3)
  }else{
    print(paste("too many variants in", loop[i], sep = " "))
  }
     }

print("Finished")


all_combo_peps<-final ##output from loop above
library(stringr)
temppep_cellline<-temppep_2[grep(sample_name,str_split_fixed(names(temppep_2),"[_ ]+",3)[,2])]#do for all cell lines to add back in single variant peptides
modpep=unlist(all_combo_peps) 
all_modpeps=modpep[[1]]
for(i in 2:length(modpep)){
  all_modpeps=c(all_modpeps,modpep[[i]])
}
print("converted to stringset")
modppep_cellline<-all_modpeps[grep(sample_name,str_split_fixed(names(all_modpeps),"[_ ]+",3)[,2])]#do for all cell lines to add back in single variant peptides
#add back temppep
all_custom_peps <-c(all_modpeps, temppep_cellline)
print("combining temppep and all_modpeps")

stopCluster(cl)
