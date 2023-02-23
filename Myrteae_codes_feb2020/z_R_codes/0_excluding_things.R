# renaming and excluding sequences from alignments:

# rm(list=ls())

# TV: renomear e excluir sequencias

library(phylotools)
setwd("~/Desktop/Myrteae_codes/0_alignments/")

rename<-read.table("RENAME.txt", h=T)
remove<-as.character(read.csv("REMOVE.txt")[,1])

rename.fasta(infile = "concat_final_11.02.2020.fasta", rename, 
             outfile = "SupTree_Nuc_Plastid.fasta")

rm.sequence.fasta("SupTree_Nuc_Plastid.fasta", 
                  outfile = "SupTree_Nuc_Plastid_cleaned.fasta", to.rm = remove)

