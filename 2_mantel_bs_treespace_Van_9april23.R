############################################################
# Mantel tests, bootstrap medians and treespace in Myrteae #
############################################################

# rm(list=ls()) # uncomment to clean environment

# load packages
library(treespace)
library(tidyverse)
library(stringr)
library(vegan)
library(ape)
library(robustbase)
library(scales)
library(phyloch)


### wd <- "C:/Users/VGS/Desktop/Papers_Van/andamento/SuperTree/Myrteae_codes/1_subsets"
setwd("C:/Users/VGS/Desktop/Papers_Van/andamento/SuperTree/Thaís_codes/8april2023/2_Mantel_bs")


#############################
###### 1. MANTEL TESTS ######
# set working directory
#setwd(paste(wd, "mantel", sep=""))

######################

# tests similarities among topologies of different markers #

#### ALL 53 ####
list.files(pattern="_all53_tree") -> tree
sapply(tree, read.tree, simplify=F) -> trees
sub("_all53_tree", "", tree) -> labs

vcv(trees[[which(labs=="full")]])->tree_full
tree_full<-tree_full[order(rownames(tree_full)), order(colnames(tree_full))]

result_mantel<-list()
for(i in 1:length(trees)){
  vcv(trees[[i]])->t0
  tree[i]->name_tree
  t0<-t0[order(rownames(t0)), order(colnames(t0))]
  result<-mantel(t0, tree_full, method="pearson", permutations=1000)
  result_mantel[[i]]<-c(name_tree, result)
}
results_final_all53<-c()
n0 <- c()
for(u in 1:length(result_mantel)){
  x0<- round(as.numeric(result_mantel[[u]][4]) ,2)
  results_final_all53[u] <- paste(str_split(unlist(result_mantel[[u]][1]), "_")[[1]][1], 
                                  x0, sep=" ")
  n0[u] <- str_split(unlist(result_mantel[[u]][1]), "_")[[1]][1]
}
names(results_final_all53) <- n0

#### Eugeninae ####

list.files(pattern="_Eugeninae_tree") -> tree
sapply(tree, read.tree, simplify=F) -> trees
sub("_Eugeninae_tree", "", tree) -> labs

vcv(trees[[which(labs=="full")]])->tree_full
tree_full<-tree_full[order(rownames(tree_full)), order(colnames(tree_full))]

result_mantel<-list()
for(i in 1:length(trees)){
  vcv(trees[[i]])->t0
  tree[i]->name_tree
  t0<-t0[order(rownames(t0)), order(colnames(t0))]
  result<-mantel(t0, tree_full, method="pearson", permutations=1000)
  result_mantel[[i]]<-c(name_tree, result)
}

results_final_eug<-c()
n0 <- c()
for(u in 1:length(result_mantel)){
  x0<- round(as.numeric(result_mantel[[u]][4]) ,2)
  results_final_eug[u] <- paste(str_split(unlist(result_mantel[[u]][1]), "_")[[1]][1], 
                                x0, sep=" ")
  n0[u] <- str_split(unlist(result_mantel[[u]][1]), "_")[[1]][1]
}
names(results_final_eug) <- n0

#### Myrciinae ####
list.files(pattern="_Myrciinae_tree") -> tree
sapply(tree, read.tree, simplify=F) -> trees
sub("_Myrciinae_tree", "", tree) -> labs

vcv(trees[[which(labs=="full")]])->tree_full
tree_full<-tree_full[order(rownames(tree_full)), order(colnames(tree_full))]

result_mantel<-list()
for(i in 1:length(trees)){
  vcv(trees[[i]])->t0
  tree[i]->name_tree
  t0<-t0[order(rownames(t0)), order(colnames(t0))]
  result<-mantel(t0, tree_full, method="pearson", permutations=1000)
  result_mantel[[i]]<-c(name_tree, result)
}

results_final_myr<-c()
n0 <- c()
for(u in 1:length(result_mantel)){
  x0<- round(as.numeric(result_mantel[[u]][4]) ,2)
  results_final_myr[u] <- paste(str_split(unlist(result_mantel[[u]][1]), "_")[[1]][1], 
                                x0, sep=" ")
  n0[u] <- str_split(unlist(result_mantel[[u]][1]), "_")[[1]][1]
}
names(results_final_myr) <- n0


#####################################
## Boxplot with bootstrap values ####
# shows which trees have higher resolution #

# Defining color code #
markers <- c("trnlF","rpl16","psba",
             "matK","ndhf","rpl32",
             "ITS","ETS","trnq",
             "nuclear","plastid","full")
colors <- c("#c7e9b4","#edf8b1","#1d91c0",
            "#081d58","#7fcdbb","#225ea8",
            "#762a83","#c2a5cf","#253494",
            "#6e016b","#41b6c4","grey")

names(colors) <- markers

# its   #762a83
# ets   #c2a5cf

# rpl16    #edf8b1
# trnl       #c7e9b4
# ndhf    #7fcdbb
# psba   #1d91c0
# rpl32   #225ea8
# trnQ    #253494
# matk   #081d58

#41b6c4  - all plastid
#88419d or #6e016b - all nuclear

# adicionais se precisar

#ffffd9
#045a8d
#016c59



#### ALL 53 ####
# bootstrap 
list.files(pattern="_all53_bs.txt") -> bs
sapply(bs, readChar, nchars=1e8, simplify=F) -> bs_all
r <-gsub("[\\(\\)]", "", regmatches(bs_all[[1]], gregexpr("\\[.*?\\]", bs_all[[1]]))[[1]])

results_median_all53 <- data.frame(matrix(NA, nrow = length(r), ncol = length(bs_all)))
for(i in 1:length(bs_all)){
  r0<-gsub("[\\(\\)]", "", regmatches(bs_all[[i]], gregexpr("\\[.*?\\]", bs_all[[i]]))[[1]])
  r0<-as.numeric(str_extract_all(r0, "[0-9]+"))
  results_median_all53[,i]<-r0
  colnames(results_median_all53)[i]<-paste(bs[i])
}

mns <- colMedians(as.matrix(results_median_all53))
results_median_all53 <- results_median_all53[,order(mns)]
names(results_median_all53) <- sub("_all53_bs.txt", "", names(results_median_all53))
boxplot(results_median_all53, col= colors, las=1)

#### EUGENIINAE ####
list.files(pattern="_Eugeninae_bs.txt") -> bs
sapply(bs, readChar, nchars=1e8, simplify=F) -> bs_all
r <-gsub("[\\(\\)]", "", regmatches(bs_all[[1]], gregexpr("\\[.*?\\]", bs_all[[1]]))[[1]])

results_median_eug <- data.frame(matrix(NA, nrow = length(r), ncol = length(bs_all)))
for(i in 1:length(bs_all)){
  r0<-gsub("[\\(\\)]", "", regmatches(bs_all[[i]], gregexpr("\\[.*?\\]", bs_all[[i]]))[[1]])
  r0<-as.numeric(str_extract_all(r0, "[0-9]+"))
  results_median_eug[,i]<-r0
  colnames(results_median_eug)[i]<-paste(bs[i])
}

mns <- colMedians(as.matrix(results_median_eug))
results_median_eug <- results_median_eug[,order(mns)]
names(results_median_eug) <- sub("_Eugeninae_bs.txt", "", names(results_median_eug))
boxplot(results_median_eug, col= colors[which(names(colors) %in% names(results_median_eug))], las=1)

#### MYRCIINAE ####
list.files(pattern="_Myrciinae_bs.txt") -> bs
sapply(bs, readChar, nchars=1e8, simplify=F) -> bs_all
r <-gsub("[\\(\\)]", "", regmatches(bs_all[[1]], gregexpr("\\[.*?\\]", bs_all[[1]]))[[1]])

results_median_myr <- data.frame(matrix(NA, nrow = length(r), ncol = length(bs_all)))
for(i in 1:length(bs_all)){
  r0<-gsub("[\\(\\)]", "", regmatches(bs_all[[i]], gregexpr("\\[.*?\\]", bs_all[[i]]))[[1]])
  r0<-as.numeric(str_extract_all(r0, "[0-9]+"))
  results_median_myr[,i]<-r0
  colnames(results_median_myr)[i]<-paste(bs[i])
}

mns <- colMedians(as.matrix(results_median_myr))
results_median_myr <- results_median_myr[,order(mns)]
names(results_median_myr) <- sub("_Myrciinae_bs.txt", "", names(results_median_myr))
boxplot(results_median_myr, col= colors[which(names(colors) %in% names(results_median_myr))], las=1)


######################
##### Treespace ######
# infers differences in both resolution and topology #

setwd("C:/Users/VGS/Desktop/Papers_Van/andamento/SuperTree/Thaís_codes/8april2023/3_treespace")

#### ALL 53 ####
#### treespace 1 - markers ####
list.files(pattern="_all53_trees.txt") -> files
sub("_all53_trees.txt", "", files[-c(2,6,7)]) -> labs
sapply(files[-c(2,6,7)], read.tree, simplify = F) -> trees
# selecting 100 trees
n=100
for (i in 1:length(trees)) {
  trees[[i]] -> t0
  sample(t0, n, replace = F) -> t0
  t0 -> trees[[i]]
}
## name trees
for (i in 1:length(trees)) {
  trees[[i]] -> t0
  paste(labs[i], 1:length(t0), sep="-") -> names(t0)
  t0 -> trees[[i]]
}
.compressTipLabel(unlist(trees, recursive=F)) -> trees
unlist(lapply(strsplit(names(trees), "\\."), "[", 3)) -> names(trees)
trees -> all.trees
.compressTipLabel(all.trees) -> all.trees
all.trees
names(trees) -> names(all.trees)
names(all.trees) -> qualifier_all53_markers
unlist(lapply(strsplit(qualifier_all53_markers, "-"), "[", 1)) -> qualifier_all53_markers
# Root trees
outgroup = "Eucalyptus.perriniana_Lucas.283"
lapply(all.trees, root, outgroup, resolve.root=T) -> all.trees
.compressTipLabel(all.trees) -> all.trees
# Run treespace
res_all53_markers <- treespace(all.trees, nf=3)

#### treespace 2 - partitions ####
list.files(pattern="_all53_trees.txt") -> files
sub("_all53_trees.txt", "", files[c(2,6,7)]) -> labs
sapply(files[c(2,6,7)], read.tree, simplify = F) -> trees
# selecting 100 trees
n=100
for (i in 1:length(trees)) {
  trees[[i]] -> t0
  sample(t0, n, replace = F) -> t0
  t0 -> trees[[i]]
}
# name trees
for (i in 1:length(trees)) {
  trees[[i]] -> t0
  paste(labs[i], 1:length(t0), sep="-") -> names(t0)
  t0 -> trees[[i]]
}
.compressTipLabel(unlist(trees, recursive=F)) -> trees
unlist(lapply(strsplit(names(trees), "\\."), "[", 3)) -> names(trees)
trees -> all.trees
.compressTipLabel(all.trees) -> all.trees
all.trees
names(trees) -> names(all.trees)
names(all.trees) -> qualifier_all53_part
unlist(lapply(strsplit(qualifier_all53_part, "-"), "[", 1)) -> qualifier_all53_part
# Root trees
outgroup = "Eucalyptus.perriniana_Lucas.283"
lapply(all.trees, root, outgroup, resolve.root=T) -> all.trees
.compressTipLabel(all.trees) -> all.trees
# Run treespace
res_all53_part <- treespace(all.trees, nf=3)


#### EUGENINAE ####
#### treespace 1 - markers ####
list.files(pattern="_Eugeninae_trees.txt") -> files
sub("_Eugeninae_trees.txt", "", files[-c(1,3,4)]) -> labs
sapply(files[-c(1,3,4)], read.tree, simplify = F) -> trees
# selecting 100 trees
n=100
for (i in 1:length(trees)) {
  trees[[i]] -> t0
  sample(t0, n, replace = F) -> t0
  t0 -> trees[[i]]
}
# name trees
for (i in 1:length(trees)) {
  trees[[i]] -> t0
  paste(labs[i], 1:length(t0), sep="-") -> names(t0)
  t0 -> trees[[i]]
}
.compressTipLabel(unlist(trees, recursive=F)) -> trees
unlist(lapply(strsplit(names(trees), "\\."), "[", 3)) -> names(trees)
trees -> all.trees
.compressTipLabel(all.trees) -> all.trees
all.trees
names(trees) -> names(all.trees)
names(all.trees) -> qualifier_eug_markers
unlist(lapply(strsplit(qualifier_eug_markers, "-"), "[", 1)) -> qualifier_eug_markers
# Root trees
outgroup = "Eucalyptus.perriniana_Lucas.283"
lapply(all.trees, root, outgroup, resolve.root=T) -> all.trees
.compressTipLabel(all.trees) -> all.trees
# Run treespace
res_eug_markers <- treespace(all.trees, nf=3)

#### treespace  2 - partitions ####
list.files(pattern="_Eugeninae_trees.txt") -> files
sub("_Eugeninae_trees.txt", "", files[c(1,3,4)]) -> labs
sapply(files[c(1,3,4)], read.tree, simplify = F) -> trees
# selecting 100 trees
n=100
for (i in 1:length(trees)) {
  trees[[i]] -> t0
  sample(t0, n, replace = F) -> t0
  t0 -> trees[[i]]
}
# name trees
for (i in 1:length(trees)) {
  trees[[i]] -> t0
  paste(labs[i], 1:length(t0), sep="-") -> names(t0)
  t0 -> trees[[i]]
}
.compressTipLabel(unlist(trees, recursive=F)) -> trees
unlist(lapply(strsplit(names(trees), "\\."), "[", 3)) -> names(trees)
trees -> all.trees
.compressTipLabel(all.trees) -> all.trees
all.trees
names(trees) -> names(all.trees)
names(all.trees) -> qualifier_eug_part
unlist(lapply(strsplit(qualifier_eug_part, "-"), "[", 1)) -> qualifier_eug_part
# Root trees
outgroup = "Eucalyptus.perriniana_Lucas.283"
lapply(all.trees, root, outgroup, resolve.root=T) -> all.trees
.compressTipLabel(all.trees) -> all.trees
# Run treespace
res_eug_part <- treespace(all.trees, nf=3)


#### Myrciinae ####
#### treespace 1 - markers ####
list.files(pattern="_Myrciinae_trees.txt") -> files
sub("_Myrciinae_trees.txt", "", files[-c(1,4,5)]) -> labs
sapply(files[-c(1,4,5)], read.tree, simplify = F) -> trees
# selecting 100 trees
n=100
for (i in 1:length(trees)) {
  trees[[i]] -> t0
  sample(t0, n, replace = F) -> t0
  t0 -> trees[[i]]
}
# name trees
for (i in 1:length(trees)) {
  trees[[i]] -> t0
  paste(labs[i], 1:length(t0), sep="-") -> names(t0)
  t0 -> trees[[i]]
}
.compressTipLabel(unlist(trees, recursive=F)) -> trees
unlist(lapply(strsplit(names(trees), "\\."), "[", 3)) -> names(trees)
trees -> all.trees
.compressTipLabel(all.trees) -> all.trees
all.trees
names(trees) -> names(all.trees)
names(all.trees) -> qualifier_myr_markers
unlist(lapply(strsplit(qualifier_myr_markers, "-"), "[", 1)) -> qualifier_myr_markers
# Root trees
outgroup = "Eucalyptus.perriniana_Lucas.283"
lapply(all.trees, root, outgroup, resolve.root=T) -> all.trees
.compressTipLabel(all.trees) -> all.trees
# Run treespace
res_myr_markers <- treespace(all.trees, nf=3)

#### treespace 2 - partitions ####
list.files(pattern="_Myrciinae_trees.txt") -> files
sub("_Myrciinae_trees.txt", "", files[c(1,4,5)]) -> labs
sapply(files[c(1,4,5)], read.tree, simplify = F) -> trees
# selecting 100 trees
n=100


for (i in 1:length(trees)) {
  trees[[i]] -> t0
  sample(t0, n, replace = F) -> t0
  t0 -> trees[[i]]
}
# name trees
for (i in 1:length(trees)) {
  trees[[i]] -> t0
  paste(labs[i], 1:length(t0), sep="-") -> names(t0)
  t0 -> trees[[i]]
}
.compressTipLabel(unlist(trees, recursive=F)) -> trees
unlist(lapply(strsplit(names(trees), "\\."), "[", 3)) -> names(trees)
trees -> all.trees
.compressTipLabel(all.trees) -> all.trees
all.trees
names(trees) -> names(all.trees)
names(all.trees) -> qualifier_myr_part
unlist(lapply(strsplit(qualifier_myr_part, "-"), "[", 1)) -> qualifier_myr_part

# Root trees
outgroup = "Eucalyptus.perriniana_Lucas.283"
lapply(all.trees, root, outgroup, resolve.root=T) -> all.trees
.compressTipLabel(all.trees) -> all.trees
# Run treespace
res_myr_part <- treespace(all.trees, nf=3)


######## PLOT ########
#setwd(paste(wd, "final/", sep=""))

pdf(file="teste.pdf", width = 16, height = 7)
par(mfrow=c(3,3), margin(1,1,1,1))
# all 53
boxplot(results_median_all53,col=colors[order(match(names(colors),names(results_median_all53)))],
        las=1, cex.axis=0.75, cex.lab=0.75,
        xlab=paste(results_final_all53[order(match(names(results_final_all53),
                                                 names(results_median_all53)))],collapse = "; "))
title(main="All markers (53 tips)")

s.class(res_all53_markers$pco$li, fac=as.factor(qualifier_all53_markers), cellipse = 0, 
        cstar=1, label=NULL, pch=19,cgrid=1, cpoint = 0.5, 
        col=colors[levels(as.factor(qualifier_all53_markers))])
s.chull(res_all53_markers$pco$li, fac=as.factor(qualifier_all53_markers), add.plot = T, 
        clabel=1.2, col=colors[levels(as.factor(qualifier_all53_markers))], 
        label=names(colors[levels(as.factor(qualifier_all53_markers))]), 
        optchull=1)

s.class(res_all53_part$pco$li, fac=as.factor(qualifier_all53_part), cellipse = 0, cstar=1, 
        label=NULL, pch=19,cgrid=1, cpoint = 0.5, 
        col=colors[levels(as.factor(qualifier_all53_part))])
s.chull(res_all53_part$pco$li, fac=as.factor(qualifier_all53_part), 
        add.plot = T, clabel=1.2, col=colors[levels(as.factor(qualifier_all53_part))], 
        label=names(colors[levels(as.factor(qualifier_all53_part))]), 
        optchull=1)

# Myrciinae
boxplot(results_median_myr,col=colors[order(match(names(colors),names(results_median_myr)))],
        las=1, cex.axis=0.75,  cex.lab=0.75,
        xlab=paste(results_final_myr[order(match(names(results_final_myr),
                                                                       names(results_median_myr)))],collapse = "; "))
title(main="Myrciinae markers (169 tips)")
s.class(res_myr_markers$pco$li, fac=as.factor(qualifier_myr_markers), cellipse = 0, cstar=1, 
        label=NULL, pch=19, cgrid=1, cpoint = 0.5, 
        col=colors[levels(as.factor(qualifier_myr_markers))])
s.chull(res_myr_markers$pco$li, fac=as.factor(qualifier_myr_markers), add.plot = T, 
        clabel=1.2, col=colors[levels(as.factor(qualifier_myr_markers))], 
        label=names(colors[levels(as.factor(qualifier_myr_markers))]), 
        optchull=1)

s.class(res_myr_part$pco$li, fac=as.factor(qualifier_myr_part), cellipse = 0, 
        cstar=1, label=NULL, pch=19, 
        cgrid=1, cpoint = 0.5, col=colors[levels(as.factor(qualifier_myr_part))])
s.chull(res_myr_part$pco$li, fac=as.factor(qualifier_myr_part), add.plot = T, 
        clabel=1.2, col=colors[levels(as.factor(qualifier_myr_part))], 
        label=names(colors[levels(as.factor(qualifier_myr_part))]), 
        optchull=1)

# Eugeniinae
boxplot(results_median_eug, col=colors[order(match(names(colors),names(results_median_eug)))],
        las=1, cex.axis=0.75, cex.lab=0.75, 
        xlab=paste(results_final_eug[order(match(names(results_final_eug), 
                                                 names(results_median_eug)))],collapse = "; "))
title(main="Eugeniinae markers (197 tips)")
s.class(res_eug_markers$pco$li, fac=as.factor(qualifier_eug_markers), cellipse = 0, cstar=1, 
        label=NULL, pch=19, cgrid=1, cpoint = 0.5, 
        col=colors[levels(as.factor(qualifier_eug_markers))])
s.chull(res_eug_markers$pco$li, fac=as.factor(qualifier_eug_markers), add.plot = T, 
        clabel=1.2, col=colors[levels(as.factor(qualifier_eug_markers))], 
        label=names(colors[levels(as.factor(qualifier_eug_markers))]), 
        optchull=1)

s.class(res_eug_part$pco$li, fac=as.factor(qualifier_eug_part), cellipse = 0, 
        cstar=1, label=NULL, pch=19, 
        cgrid=1, cpoint = 0.5, col=colors[levels(as.factor(qualifier_eug_part))])
s.chull(res_eug_part$pco$li, fac=as.factor(qualifier_eug_part), add.plot = T, 
        clabel=1.2, col=colors[levels(as.factor(qualifier_eug_part))], 
        label=names(colors[levels(as.factor(qualifier_eug_part))]), 
        optchull=1)

dev.off()

