
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
Ec <- read.fasta("sample.fa")
Myseq <- Ec [[61]]
str(Myseq)
#calculating the GC 
seqinr::GC(Myseq)
#length of myseq in basepairs
seqinr::getLength(Myseq)


#Question 3 - Create Blast function with provided R function and do blast searches
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/mutblast_functions.R",
              destfile = "mutblast_functions.R")
source("mutblast.R")
# Use Blast to identify ecoli sequense
mutb <- myblastn_tab (myseq = Myseq, db ="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
head(mutb)
str(mutb)
mysseqid <- as.character(mutb$sseqid)
mysseqid
# top 3 hit genes with percent identity, E-values and bit scores
hits <- as.character(mutb$sseqid[1:3])
hits
hits <- c(hits, "AAC76307")
hits

#Question 4 - Find number of mismatches  between original and mutated seq?
mutator(Myseq, 50)
Myseq_mut <- mutator(Myseq, 50)
myblastn_tab( Myseq_mut, db = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")


#Question 5 - Perform mutating and blasting sequence to determine number of proprtion of sites need to be altered

Result <- function (Myseq, nmut) {
  mutseq <- mutator( Myseq, nmut)
  res <- myblastn_tab(myseq = mutseq, db = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa")
if ( is.null(res) ) {
  myres=0 } else { myres=1
  }
return(myres)
}
replicate(n = 100, expr = Result(Myseq, 100) )
mean(replicate(n=100,expr = Result(Myseq, 100) ) )
str(mean)
  finalres <- function (n) {
    mean(replicate(n, expr = Result(Myseq, 100) ) )
}
n <- c(10, 50,100, 150, 200, 250, 300)  
sapply(n, finalres)

#Ques 6 -  Plot a chart with increasing proprotion of mutated bases

plot(n, sapply(n,finalres), type = 'o')
  


  
