#######################
# rm(list=ls())

library(phylotools)
library(ape)
library(phangorn)

#### Subset Myrteae tree and Treespace ####
# Setting working directory #
setwd("C:/Users/VGS/Desktop/Papers_Van/andamento/SuperTree/Thaís_codes/8april2023")

# Loading full alignment #
full_alignment <- read.dna("all_concatenated_8abril2023.fasta", format="fasta", as.matrix=T)
as.alignment(full_alignment)->tips
tips$nam->tips

# Extracting each region #

ETS <- full_alignment[,1:641]
ITS <- full_alignment[,642:1613]
matK <- full_alignment[,1614:2520]
ndhf <- full_alignment[,2521:3310]
psba <- full_alignment[,3311:4750]
rpl16 <- full_alignment[,4751:6068]
rpl32 <- full_alignment[,6069:7521]
trnlF <- full_alignment[,7522:8485]
trnq <- full_alignment[,8486:10312]

####################################################
# Extracting three blocks for Results section 3.2. #

# Only tips with all regions covered (52)
all_53 <- list(rpl16, rpl32, trnq, psba, ETS, ITS, matK, trnlF, ndhf)

# Creating new aligments with selected sequences #
# ALL 53 #
names<-list()
for(i in 1:length(all_53)){
  del.rowgapsonly(all_53[[i]])->a0
  as.alignment(a0)->a0
  a0$nam->n0
  names[[i]]<-n0
}
all0 <-Reduce(intersect,names)
'%!in%' <- function(x,y)!('%in%'(x,y))
remove <- tips %!in% all0 
r<-subset(tips, remove)
rm.sequence.fasta("all_concatenated_8abril2023.fasta", outfile = "all_53.fasta", to.rm = r)

#######################################################################
# rm(list=ls())

# Creating Eugeninae alignment #
remove<-as.character(read.csv("REMOVE_Eugeninaefile.txt")[,1])

rm.sequence.fasta("all_concatenated_8abril2023.fasta", 
                  outfile = "Eugeninae_9markers.fasta", to.rm = remove)

Eugeninae_alignment <- read.dna("Eugeninae_9markers.fasta", format="fasta", as.matrix=T)
as.alignment(Eugeninae_alignment)->tips
tips$nam->tips

# Extracting each region #

ETS <- Eugeninae_alignment[,1:641]
ITS <- Eugeninae_alignment[,642:1613]
matK <- Eugeninae_alignment[,1614:2520]
ndhf <- Eugeninae_alignment[,2521:3310]
psba <- Eugeninae_alignment[,3311:4750]
rpl16 <- Eugeninae_alignment[,4751:6068]
rpl32 <- Eugeninae_alignment[,6069:7521]
trnlF <- Eugeninae_alignment[,7522:8485]
trnq <- Eugeninae_alignment[,8486:10312]

####################################################
# Extracting three blocks for Results section 3.2. #

# Only tips with five regions covered (193 Eugeninae + 3 outgroup)
all_Eugeninae <- list(rpl16, rpl32, trnq, psba, ITS)

# Creating new aligments with selected sequences #
names<-list()
for(i in 1:length(all_Eugeninae)){
  del.rowgapsonly(all_Eugeninae[[i]])->a0
  as.alignment(a0)->a0
  a0$nam->n0
  names[[i]]<-n0
}
all0 <-Reduce(intersect,names)
'%!in%' <- function(x,y)!('%in%'(x,y))
remove <- tips %!in% all0 
r<-subset(tips, remove)
rm.sequence.fasta("Eugeninae_9markers.fasta", outfile = "Eugeninae_5markers_complete.fasta", to.rm = r)



#######################################################################
# rm(list=ls())

# Creating Myrciinae alignment #
remove<-as.character(read.csv("REMOVE_Myrciinaefile.txt")[,1])

rm.sequence.fasta("all_concatenated_8abril2023.fasta", 
                  outfile = "Myrciinae_9markers.fasta", to.rm = remove)

Myrciinae_alignment <- read.dna("Myrciinae_9markers.fasta", format="fasta", as.matrix=T)
as.alignment(Myrciinae_alignment)->tips
tips$nam->tips

# Extracting each region #

ETS <- Myrciinae_alignment[,1:641]
ITS <- Myrciinae_alignment[,642:1613]
matK <- Myrciinae_alignment[,1614:2520]
ndhf <- Myrciinae_alignment[,2521:3310]
psba <- Myrciinae_alignment[,3311:4750]
rpl16 <- Myrciinae_alignment[,4751:6068]
rpl32 <- Myrciinae_alignment[,6069:7521]
trnlF <- Myrciinae_alignment[,7522:8485]
trnq <- Myrciinae_alignment[,8486:10312]

####################################################
# Extracting three blocks for Results section 3.2. #

# Only tips with five regions covered (188 Myrciinae + 3 outgroup)
all_Myrciinae <- list(trnlF, ndhf, trnq, psba, ITS)

# Creating new aligments with selected sequences #
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
rm.sequence.fasta("Myrciinae_9markers.fasta", outfile = "Myrciinae_5markers_complete.fasta", to.rm = r)





#########
# saving separate phylip files for tree inference #

all_53 <- read.dna("all_53.fasta", format="fasta", as.matrix=T)
all_myrciinae <- read.dna("Myrciinae_5markers_complete.fasta", format="fasta", as.matrix=T)
all_eugeniinae <- read.dna("Eugeninae_5markers_complete.fasta", format="fasta", as.matrix=T)

# All 53
ETS <- 1:641
ITS <- 642:1613
matK <- 1614:2520
ndhf <- 2521:3310
psba <- 3311:4750
rpl16 <- 4751:6068
rpl32 <- 6069:7521
trnlF <- 7522:8485
trnq <- 8486:10312
nuclear <- 1:1613
plastid <-  1614:10312


setwd("C:/Users/VGS/Desktop/Papers_Van/andamento/SuperTree/Thaís_codes/8april2023/1_subsets/")

write.phyDat(all_53[,matK], file="all53_matK.txt")
write.phyDat(all_53[,trnlF], file="all53_trnlF.txt")
write.phyDat(all_53[,ndhf], file="all53_ndhf.txt")
write.phyDat(all_53[,rpl16], file="all53_rpl16.txt")
write.phyDat(all_53[,rpl32], file="all53_rpl32.txt")
write.phyDat(all_53[,trnq], file="all53_trnq.txt")
write.phyDat(all_53[,psba], file="all53_psba.txt")
write.phyDat(all_53[,ETS], file="all53_ETS.txt")
write.phyDat(all_53[,ITS], file="all53_ITS.txt")
write.phyDat(all_53[,nuclear], file="all53_nuclear.txt")
write.phyDat(all_53[,plastid], file="all53_plastid.txt")
write.phyDat(all_53, file="all53_full.txt")

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

