#Load Packages
library(VariantAnnotation)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)
library(svMisc)
library(pbapply)
library(Biostrings)

#Load objects and inputs
load("Tools/seqinfo1.RData")
load("Tools/ccds.RData")
load("Tools/common_snps_missense.RData")
gtf_file <- "Tools/gencode.v28.annotation.gtf"
prot<-readAAStringSet("Tools/gencode.v28.pc_translations.fa", format = "fasta", use.names = T)
directory<-getwd()                               
db_directory<-paste(getwd(),"/Custom_Databases/",sep="")

createTxDB <- function(x, y) {
  txdb <- makeTxDbFromGFF(x,
                          format = c("gtf"),
                          dataSource = "gencodev28",
                          organism = NA,
                          taxonomyId = NA,
                          circ_seqs = DEFAULT_CIRC_SEQS,
                          chrominfo = y,
                          miRBaseBuild = NA,
                          metadata = NULL,
                          dbxrefTag = "tx_id")
return(txdb)
}

print("Make TxDB object and get Gene/Transcipt IDs") 
txdb<- createTxDB(gtf_file, seqinfo1)
txid = keys(txdb,"TXID")
df = biomaRt::select(txdb, txid, columns= c("TXNAME","GENEID"),"TXID")
new_txids<-df[which(df$TXNAME%in%ccds$nucleotide_ID),1]
new_df<-df[which(df$TXNAME%in%ccds$nucleotide_ID),]
new_df<-cbind(new_df,ccds[match(new_df$TXNAME,ccds$nucleotide_ID),10])
colnames(new_df)[4]<-"Uniprot_ID"
print("Done")

print("Get protein sequences and subset")
prot_id<-str_split_fixed(names(prot), "\\|",3)[,2]
fasta_headers<-str_split_fixed(names(prot), "\\|",8)   ###Gencode FASTA headers
table(as.character(ccds$nucleotide_ID)%in%unique(prot_id))
names(prot)<-prot_id ##make the names txname
pep<-prot[which(names(prot)%in%ccds$nucleotide_ID),] ##the peps that correspond to the ccds CCDS transcript id's
pep_txids<-biomaRt::select(txdb, names(pep), columns= c("TXID","GENEID"),"TXNAME")
names(pep)<-pep_txids[,3]
subseq(pep, start =width(pep)+1)<-"*"   ###puts an asterisk at the end of every AA sequence
rm(txdb)
print("Done")

print("Loop through VCFs and get predicted coding changes")

fnames<-list.files(pattern=".vcf*",full.names = T)##change name _filtered or _Aligned
fnames2<-gsub("_Aligned\\..*","",fnames)
fnames2<-gsub("./","",fnames2)
list<-list()

for (i in 1:length(fnames)){
  library(VariantAnnotation)
  fl <- fnames[i]
  vcf <- VariantAnnotation::readVcf(fl, "hg38")
  txdb <- createTxDB(gtf_file, seqinfo(vcf))
  genome <- getBSgenome("hg38")
  coding <- VariantAnnotation::predictCoding(vcf, txdb, seqSource = genome)
  all_missense <- coding[coding$CONSEQUENCE == "nonsynonymous", ]
  all_missense <- all_missense[all_missense$TXID %in% new_txids, ]
  snps_to_keep <- paste(all_missense$TXID, all_missense$REFAA, all_missense$VARAA)
  snps_common <- paste(common_snps$TXID, common_snps$REFAA, common_snps$VARAA)
  all_missense_snps <- all_missense[snps_to_keep %in% snps_common, ]
  all_missense_rare <- all_missense[!snps_to_keep %in% snps_common, ]
  param_range_id <- toupper(paste0(sub("_RNA.*", "", fnames2[i]), "_RARE_RNA"))
  all_missense_rare$paramRangeID <- param_range_id
  param_range_id <- toupper(paste0(sub("_RNA.*", "", fnames2[i]), "_SNP_RNA"))
  all_missense_snps$paramRangeID <- param_range_id
  all_missense <- c(all_missense_rare, all_missense_snps)
  all_missense$paramRangeID <- as.factor(all_missense$paramRangeID)
  list <- c(list, assign(paste("all_missense", fnames2[i], sep = "_"), all_missense))
}
print("Done")

#RNA<-do.call(c,RNA_missense_all)#Exome<-do.call(c,Exome_missense_all)#all_missense<-c(RNA,Exome)

print("Get a dataframe of predicted coding changes with gene names and transcript IDs")
changelist <- data.frame(all_missense, row.names = NULL)
changelist$ALT <- unlist(lapply(changelist$ALT, function(x) as.character(unlist(x))))
add <- fasta_headers[, 7][match(changelist$GENEID, fasta_headers[, 3])]
changelist$add <- add
identifiers <- paste(changelist$add, paste("p.", changelist$REFAA, changelist$PROTEINLOC, changelist$VARAA, sep = ""), sep = "_")
changelist$identifiers <- identifiers
changelist2 <- data.frame(cell_line = gsub("_.*", "", changelist$paramRangeID), gene = changelist$add)
ref_length <- nchar(changelist$REF)
var_length <- nchar(changelist$ALT)
changelist <- cbind(changelist, ref_length, var_length, changelist2$cell_line)
changelist <- changelist[changelist$ref_length == 1 & changelist$var_length == 1, ]
ensembl_txids<-fasta_headers[,2][match(changelist$GENEID, fasta_headers[,3])] ##get ensemble transcript ids
changelist<-cbind(changelist,ensembl_txids)
print("Done")

print("Replace Reference AA and combine exome and RNA variant proteins")
setwd(directory)
source("Tools/replaceAA.R")
print("Done with loop")


all_custom_peps<-temppep_2
print("Formatting headers")
names(all_custom_peps)<- paste("sp","|",fasta_headers[match(df[match(str_split_fixed(names(all_custom_peps)," ",3)[,1],df$TXID),3],fasta_headers[,2]),3], "|",fasta_headers[match(df[match(str_split_fixed(names(all_custom_peps)," ",3)[,1],df$TXID),3],fasta_headers[,2]),7],"|", df[match(str_split_fixed(names(all_custom_peps)," ",3)[,1],df$TXID),3], "|", names(all_custom_peps), "|", sep = "")
print("Done formatting headers")


library(pbapply)
findlower<-function(amino_string){
  amino_string=str_split(amino_string,"")
  amino_string=unlist(amino_string)
  amino_string=as.character(amino_string)
  ret=sapply(amino_string, function(x) tolower(x) == x & tolower(x) != "*")
  ret=which(ret)
  return(ret)
}

print("Put variant locations in headers")
matrix<-as.data.frame(all_custom_peps)
mat<-as.matrix(matrix)
test<-pbsapply(mat[,1],function(y) findlower(y))
test<-pbsapply(test, function(w) paste(w,names(w),sep = ""))
test=pbsapply(test, function(z) paste(z, collapse = ""))
test=unlist(test)
names=unname(test)
print("Done putting locations")

print("Put variant locations in headers 2")
names(all_custom_peps)<-paste(paste(names(all_custom_peps),test,sep = ""),"|", sep = "")
list<-str_split_fixed(names(all_custom_peps),"\\|",7)
list2<-str_split_fixed(list[,6],",",7) 

order<-apply(list2, 1, function(x){cf2<-as.numeric(gsub("[^[:digit:]]", "", substr(x,stop = nchar(x),start = nchar(x)-5)))
names(cf2) <- seq_along(cf2)
gsub("[^[:digit:]]", "", substr(x,stop = nchar(x),start = nchar(x)-5))[as.numeric(names(sort(cf2)))]
x[as.numeric(names(sort(cf2)))]})
ordered<-lapply(order,function(x) paste(x,collapse = ","))
bind<-do.call(rbind,ordered)
list[,6]<-bind
names<-apply(list,1,function(x) paste(x,collapse = "|"))

names(all_custom_peps)<-names
subseq(all_custom_peps, start =width(all_custom_peps))<-""   ##takes the asterisk out
print("Done putting locations")

print("Add in Uniprot ID")
names(all_custom_peps)<-paste(str_split_fixed(names(all_custom_peps),"\\|",2)[,1],ccds[match(str_split_fixed(names(all_custom_peps),"\\|",7)[,4], ccds$nucleotide_ID),10],str_split_fixed(names(all_custom_peps),"\\|",2)[,2],sep = "|")
all_custom_peps<-all_custom_peps[!duplicated(names(all_custom_peps)),]
names(all_custom_peps)<-gsub(" ","-",names(all_custom_peps))

print("Write XStringSet")
setwd(directory)
system("mkdir Custom_Databases")
Biostrings::writeXStringSet(all_custom_peps, filepath = "Custom_Databases/CellLine_custom.fa")
print("Done")
setwd(db_directory)
system("python3 ../Tools/Automated_DB_Processor.py")
system("python3 ../Tools/Automated_Header_Processor.py")

fasta<-readAAStringSet(file="2TS_CellLine_edited.fa") ##2TS_blank_ goes here
names(fasta)<-sub("sp\\|","", names(fasta))
names(fasta)<-sub("-.*","", names(fasta))
names(fasta)<-sub("\\|.*","", names(fasta))
fasta<-fasta[!duplicated(fasta)] #remove same tryptic peptides derived from different transcripts
rev_fasta<-Biostrings::reverse(fasta)
names(rev_fasta)<-paste0("rev_",names(rev_fasta))
final_fasta<-c(fasta,rev_fasta)
writeXStringSet(final_fasta, filepath = "2TS_CellLine_edited_simple_dedup_rev.fa", format = "fasta", width=800)

fasta2<-readAAStringSet(file="2TS_CellLine_edited.fa") ##2TS_blank_ goes here
fasta2<-fasta2[!duplicated(fasta2)]
rev_fasta2<-Biostrings::reverse(fasta2)
names(rev_fasta2)<-paste0("rev_",names(rev_fasta2))
final_fasta2<-c(fasta2,rev_fasta2)
writeXStringSet(final_fasta2, filepath = "2TS_CellLine_edited_dedup_rev.fa", format = "fasta", width=800)

setwd(db_directory)
system("dd if=2TS_CellLine_edited_simple_dedup_rev.fa of=2TS_CellLine_edited_simple_dedup_rev_upper.fa conv=ucase")
system("dd if=2TS_CellLine_edited_dedup_rev.fa of=2TS_CellLine_edited_dedup_rev_upper.fa conv=ucase")