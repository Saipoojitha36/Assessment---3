
library("seqinr")
library("R.utils")
library("rBLAST")
library("ape")
library("ORFik")
library("Biostrings")

#Question 1 - Download whole set of E.coli gene DNA sequences

download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",
             destfile = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")

# decompressing the file
R.utils::gunzip("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz", overwrite = TRUE)
# Creating a blast database
makeblastdb("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa", dbtype = "nucl", "-parse_seqids")
#number of sequences in E.coli is 4140

#Question 2 - Download sample fasta sequences
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa",
              destfile = "sample.fa")
# read into r script
E.c <- read.fasta("sample.fa")
myseq <- E.c [[61]]
str(myseq)
#calculating the GC 
seqinr::GC(myseq)
#length of myseq in basepairs
seqinr::getLength(myseq)


#Question 3 - 
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/mutblast_functions.R",
              destfile = "mutblast_functions.R")
source("mutblast.R")
mutb <- myblastn_tab (myseq = myseq, db ="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
head(mutb)
str(mutb)
mysseqid <- as.character(mutb$sseqid)
mysseqid
hits <- as.character(mutb$sseqid[1:3])
hits
hits <- c(hits, "AAC76307")
hits

#Question 4 - Find number of mismatches  between original and mutated seq?

mutator <- function (myseq, nmut)
myseq_mt <- mutator(myseq, 50)
myblastn_tab(myseq_mt, db = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")


#Question 5 -

Result <- function(myseq, nmut) {
  seqmut <- mutator(myseq=myseq, nmut=nmut)
res <- myblastn_tab(myseq=seqmut, db = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
if(is.null(res)) {
  myres=0} else {myres=1
  }
return(myres)
}
replicate(n = 100, expr = Result(myseq, 100) )
mean(replicate(n=100,expr = Result(myseq, 100) ) )
str(mean)
  finalres <- function (n) {
    mean(replicate(n, expr = Result(myseq, 100) ) )
}
n <- c(10, 100, 200, 300)  
sapply(n, finalres)
  


  
