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
#devtools::install_github("fmichonneau/phyloch")

setwd("~/Desktop/myrteae_tree")


#-------------------------------------------
#############################
###### 1. MANTEL TESTS ######
############################

# tests similarities among topologies of different markers #
#### ALL 53 ####
tree <- list.files("2_Mantel_bs","_all53_tree", full.names = T)
trees <- sapply(tree, read.tree, simplify=F)
labels <- sub(paste0(c("2_Mantel_bs/","_all53_tree"), collapse="|"), "", tree) 

tree_full <- vcv(trees[[grep("full",labels)]])
tree_full <- tree_full[order(rownames(tree_full)), order(colnames(tree_full))]

result_mantel <- list()
for(i in 1:length(trees)){
  t0 <- vcv(trees[[i]])
  name_tree <- tree[i]
  t0 <- t0[order(rownames(t0)), order(colnames(t0))]
  result <- mantel(t0, tree_full, method="pearson", permutations=1000)
  result_mantel[[i]]<-c(name_tree, result)
}

results_final_all53 <- c()
for(u in 1:length(result_mantel)){
  results_final_all53[u] <- round(as.numeric(result_mantel[[u]][4]) ,2)
}
names(results_final_all53) <- labels

#----------------------------
#### Eugeninae ####

tree <- list.files("2_Mantel_bs","_Eugeninae_tree", full.names = T)
trees <- sapply(tree, read.tree, simplify=F) 
labels <- sub(paste0(c("2_Mantel_bs/","_Eugeninae_tree"), collapse="|"), "", tree) 

tree_full <- vcv(trees[[grep("full",labels)]])
tree_full <- tree_full[order(rownames(tree_full)), order(colnames(tree_full))]

result_mantel <- list()
for(i in 1:length(trees)){
  t0 <- vcv(trees[[i]])
  name_tree <- tree[i]
  t0 <- t0[order(rownames(t0)), order(colnames(t0))]
  result <- mantel(t0, tree_full, method="pearson", permutations=1000)
  result_mantel[[i]] <- c(name_tree, result)
}

results_final_eug <- c()
for(u in 1:length(result_mantel)){
  results_final_eug[u] <- round(as.numeric(result_mantel[[u]][4]) ,2)
}
names(results_final_eug) <- labels

#----------------------------
#### Myrciinae ####
tree <- list.files("2_Mantel_bs","_Myrciinae_tree", full.names = T)
trees <- sapply(tree, read.tree, simplify=F) 
labels <- sub(paste0(c("2_Mantel_bs/","_Myrciinae_tree"), collapse="|"), "", tree) 

tree_full <- vcv(trees[[grep("full",labels)]])
tree_full <- tree_full[order(rownames(tree_full)), order(colnames(tree_full))]

result_mantel <- list()
for(i in 1:length(trees)){
  t0 <- vcv(trees[[i]])
  name_tree <- tree[i]
  t0 <- t0[order(rownames(t0)), order(colnames(t0))]
  result <- mantel(t0, tree_full, method="pearson", permutations=1000)
  result_mantel[[i]] <- c(name_tree, result)
}

results_final_myr <- c()
for(u in 1:length(result_mantel)){
  results_final_myr[u] <- round(as.numeric(result_mantel[[u]][4]) ,2)
}
names(results_final_myr) <- labels

results_final_all53 <- as.data.frame(results_final_all53)
results_final_eug <- as.data.frame(results_final_eug)
results_final_myr <- as.data.frame(results_final_myr)
colnames(results_final_all53) <- colnames(results_final_eug) <- colnames(results_final_myr) <- "mantel"
write.csv(rbind(results_final_all53, results_final_eug, results_final_myr), file="all_mantel_results.csv")

#####################################
## Boxplot with bootstrap values ####
# shows which trees have higher resolution #

# Defining color code #
#markers <- c("trnlF","rpl16","psba",
             # "matK","ndhf","rpl32",
             # "ITS","ETS","trnq",
             # "nuclear","plastid","full")


#----------------------------
#### ALL 53 ####
# bootstrap 

bs <- list.files(path="2_Mantel_bs", pattern="_all53_bs", full.names = T)
bs_all <- sapply(bs, readChar, nchars=1e8, simplify=F) 
r <- gsub("[\\(\\)]", "", regmatches(bs_all[[1]], gregexpr("\\[.*?\\]", bs_all[[1]]))[[1]])

results_median_all53 <- data.frame(matrix(NA, nrow = length(r), ncol = length(bs_all)))
for(i in 1:length(bs_all)){
  r0 <- gsub("[\\(\\)]", "", regmatches(bs_all[[i]], gregexpr("\\[.*?\\]", bs_all[[i]]))[[1]])
  r0 <- as.numeric(str_extract_all(r0, "[0-9]+"))
  results_median_all53[,i] <- r0
  colnames(results_median_all53)[i] <- paste(bs[i])
}

mns <- colMedians(as.matrix(results_median_all53))
results_median_all53 <- results_median_all53[,order(mns)]
names(results_median_all53) <- gsub(paste0(c("2_Mantel_bs/","_all53_bs.txt"), collapse="|"), "", names(results_median_all53))

#boxplot(results_median_all53, col= colors, las=1)

results_median_all53 <- reshape2::melt(results_median_all53)
colors <- hcl.colors(length(unique(results_median_all53$variable)), palette = "Viridis", alpha = 0.75)

all_53_bs_boxplot <- ggplot(results_median_all53, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_manual(values=colors) +
  theme_classic() +
  coord_flip() +
  theme(plot.title=element_text(size=10))

pdf("plots/bootstrap_comparison_all53.pdf", width=6, height=6)
all_53_bs_boxplot
dev.off()

#### EUGENIINAE ####
bs <- list.files(path="2_Mantel_bs", pattern="_Eugeninae_bs.txt", full.names = T)
bs_all <- sapply(bs, readChar, nchars=1e8, simplify=F) 
r <- gsub("[\\(\\)]", "", regmatches(bs_all[[1]], gregexpr("\\[.*?\\]", bs_all[[1]]))[[1]])

results_median_eug <- data.frame(matrix(NA, nrow = length(r), ncol = length(bs_all)))
for(i in 1:length(bs_all)){
  r0<-gsub("[\\(\\)]", "", regmatches(bs_all[[i]], gregexpr("\\[.*?\\]", bs_all[[i]]))[[1]])
  r0<-as.numeric(str_extract_all(r0, "[0-9]+"))
  results_median_eug[,i]<-r0
  colnames(results_median_eug)[i]<-paste(bs[i])
}

mns <- colMedians(as.matrix(results_median_eug))
results_median_eug <- results_median_eug[,order(mns)]
names(results_median_eug) <- gsub(paste0(c("2_Mantel_bs/","_Eugeninae_bs.txt"), collapse="|"), "", names(results_median_eug))
#boxplot(results_median_eug, col= colors[which(names(colors) %in% names(results_median_eug))], las=1)

results_median_eug <- reshape2::melt(results_median_eug)
colors <- hcl.colors(length(unique(results_median_eug$variable)), palette = "Viridis", alpha = 0.75)

eug_bs_boxplot <- ggplot(results_median_eug, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_manual(values=colors) +
  theme_classic() +
  coord_flip() +
  theme(plot.title=element_text(size=10))

pdf("plots/bootstrap_comparison_eug.pdf", width=6, height=6)
eug_bs_boxplot
dev.off()


#### MYRCIINAE ####
bs <- list.files(path="2_Mantel_bs", pattern="_Myrciinae_bs.txt", full.names = T)
bs_all <- sapply(bs, readChar, nchars=1e8, simplify=F) 
r <- gsub("[\\(\\)]", "", regmatches(bs_all[[1]], gregexpr("\\[.*?\\]", bs_all[[1]]))[[1]])

results_median_myr <- data.frame(matrix(NA, nrow = length(r), ncol = length(bs_all)))
for(i in 1:length(bs_all)){
  r0<-gsub("[\\(\\)]", "", regmatches(bs_all[[i]], gregexpr("\\[.*?\\]", bs_all[[i]]))[[1]])
  r0<-as.numeric(str_extract_all(r0, "[0-9]+"))
  results_median_myr[,i]<-r0
  colnames(results_median_myr)[i]<-paste(bs[i])
}

mns <- colMedians(as.matrix(results_median_myr))
results_median_myr <- results_median_myr[,order(mns)]
names(results_median_myr) <- gsub(paste0(c("2_Mantel_bs/","_Myrciinae_bs.txt"), collapse="|"),"", names(results_median_myr))
#boxplot(results_median_myr, col= colors[which(names(colors) %in% names(results_median_myr))], las=1)

results_median_myr <- reshape2::melt(results_median_myr)
colors <- hcl.colors(length(unique(results_median_myr$variable)), palette = "Viridis", alpha = 0.75)

myr_bs_boxplot <- ggplot(results_median_myr, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_manual(values=colors) +
  theme_classic() +
  coord_flip() +
  theme(plot.title=element_text(size=10))

pdf("plots/bootstrap_comparison_myr.pdf", width=6, height=6)
myr_bs_boxplot
dev.off()

