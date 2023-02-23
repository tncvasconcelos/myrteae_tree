#######################
# rm(list=ls())

library(phylotools)
library(ape)
library(phangorn)

#### Subset Myrteae tree and Treespace ####
# Setting working directory #
setwd("~/Desktop/Myrteae_codes/0_alignments/")

# Loading full alignment #
full_alignment <- read.dna("SupTree_Nuc_Plastid_cleaned.fasta", format="fasta", as.matrix=T)
as.alignment(full_alignment)->tips
tips$nam->tips

# Extracting each region #

ETS <- full_alignment[,1:705]
ITS <- full_alignment[,706:1754]
matK <- full_alignment[,1755:2653]
ndhf <- full_alignment[,2654:3545]
psba <- full_alignment[,3546:5035]
rpl16 <- full_alignment[,5036:6332]
rpl32 <- full_alignment[,6333:7804]
trnlF <- full_alignment[,7805:8942]
trnq <- full_alignment[,8943:10966]

####################################################
# Extracting three blocks for Results section 3.2. #

# Only tips with all regions covered (52)
all_52 <- list(rpl16, rpl32, trnq, psba, ETS, ITS, matK, trnlF, ndhf)
# Tips with matrix complete for five markers used in Eugeninae focused studies (329)
all_Eugeniinae <- list(rpl16, rpl32, trnq, psba, ITS) # Eugeniinae
# Tips with matrix complete for five markers used in Myrciinae focused studies (195)
all_Myrciinae <- list(trnq, psba, ITS, trnlF, ndhf) # Myrciinae 

# Creating new aligments with selected sequences #
# ALL 52 #
names<-list()
for(i in 1:length(all_52)){
  del.rowgapsonly(all_52[[i]])->a0
  as.alignment(a0)->a0
  a0$nam->n0
  names[[i]]<-n0
}
all0 <-Reduce(intersect,names)
'%!in%' <- function(x,y)!('%in%'(x,y))
remove <- tips %!in% all0 
r<-subset(tips, remove)
rm.sequence.fasta("SupTree_Nuc_Plastid_cleaned.fasta", outfile = "all_52.fasta", to.rm = r)

# Eugeniinae #
names<-list()
for(i in 1:length(all_Eugeniinae)){
  del.rowgapsonly(all_Eugeniinae[[i]])->a0
  as.alignment(a0)->a0
  a0$nam->n0
  names[[i]]<-n0
}
all0 <-Reduce(intersect,names)
'%!in%' <- function(x,y)!('%in%'(x,y))
remove <- tips %!in% all0 
r<-subset(tips, remove)
rm.sequence.fasta("SupTree_Nuc_Plastid_cleaned.fasta", outfile = "all_eugeniinae.fasta", to.rm = r)

# Myrciinae #
names<-list()
for(i in 1:length(all_Myrciinae)){
  del.rowgapsonly(all_Myrciinae[[i]])->a0
  as.alignment(a0)->a0
  a0$nam->n0
  names[[i]]<-n0
}
all0 <-Reduce(intersect,names)
'%!in%' <- function(x,y)!('%in%'(x,y))
remove <- tips %!in% all0 
r<-subset(tips, remove)
rm.sequence.fasta("SupTree_Nuc_Plastid_cleaned.fasta", outfile = "all_myrciinae.fasta", to.rm = r)


#########
# saving separate phylip files for tree inference #

all_52 <- read.dna("all_52.fasta", format="fasta", as.matrix=T)
all_myrciinae <- read.dna("all_myrciinae.fasta", format="fasta", as.matrix=T)
all_eugeniinae <- read.dna("all_eugeniinae.fasta", format="fasta", as.matrix=T)

# All 52

ETS <- 1:705
ITS <- 706:1754
matK <- 1755:2653
ndhf <- 2654:3545
psba <- 3546:5035
rpl16 <- 5036:6332
rpl32 <- 6333:7804
trnlF <- 7805:8942
trnq <- 8943:10966
nuclear <- 1:1754
plastid <-  1755:10966


setwd("~/Desktop/Myrteae_codes/1_subsets/")

write.phyDat(all_52[,matK], file="all52_matK.txt")
write.phyDat(all_52[,trnlF], file="all52_trnlF.txt")
write.phyDat(all_52[,ndhf], file="all52_ndhf.txt")
write.phyDat(all_52[,rpl16], file="all52_rpl16.txt")
write.phyDat(all_52[,rpl32], file="all52_rpl32.txt")
write.phyDat(all_52[,trnq], file="all52_trnq.txt")
write.phyDat(all_52[,psba], file="all52_psba.txt")
write.phyDat(all_52[,ETS], file="all52_ETS.txt")
write.phyDat(all_52[,ITS], file="all52_ITS.txt")
write.phyDat(all_52[,nuclear], file="all52_nuclear.txt")
write.phyDat(all_52[,plastid], file="all52_plastid.txt")
write.phyDat(all_52, file="all52_full.txt")

###
write.phyDat(all_eugeniinae[,rpl16], file="eugeniinae_rpl16.txt")
write.phyDat(all_eugeniinae[,rpl32], file="eugeniinae_rpl32.txt")
write.phyDat(all_eugeniinae[,trnq], file="eugeniinae_trnq.txt")
write.phyDat(all_eugeniinae[,psba], file="eugeniinae_psba.txt")
write.phyDat(all_eugeniinae[,ITS], file="eugeniinae_ITS.txt")
write.phyDat(all_eugeniinae[,ITS], file="eugeniinae_nuclear.txt")
write.phyDat(all_eugeniinae[,c(rpl16,rpl32,trnq,psba)], file="eugeniinae_plastid.txt")
write.phyDat(all_eugeniinae[,c(rpl16,rpl32,trnq,psba,ITS)], file="eugeniinae_full.txt")

###
write.phyDat(all_myrciinae[,trnlF], file="myrciinae_trnlF.txt")
write.phyDat(all_myrciinae[,ndhf], file="myrciinae_ndhf.txt")
write.phyDat(all_myrciinae[,trnq], file="myrciinae_trnq.txt")
write.phyDat(all_myrciinae[,psba], file="myrciinae_psba.txt")
write.phyDat(all_myrciinae[,ITS], file="myrciinae_ITS.txt")
write.phyDat(all_myrciinae[,ITS], file="myrciinae_nuclear.txt")
write.phyDat(all_myrciinae[,c(trnlF,ndhf,trnq,psba)], file="myrciinae_plastid.txt")
write.phyDat(all_myrciinae[,c(trnlF,ndhf,trnq,psba,ITS)], file="myrciinae_full.txt")

