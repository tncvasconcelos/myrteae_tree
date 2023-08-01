library(ape)
library(treespace)

setwd("~/Desktop/myrteae_tree/3_treespace/")

#####################################
## Boxplot with bootstrap values ####
# shows which trees have higher resolution #

# Defining color code #
markers <- c("trnlF","rpl16","psba",
             "matK","ndhf","rpl32",
             "ITS","ETS","trnq",
             "nuclear","plastid","full")

palettes_hcl <- hcl.pals()
colors <- hcl.colors(length(markers), palette = "Viridis", alpha = 0.75)
names(colors) <- markers


######################
##### Treespace ######
# infers differences in both resolution and topology #

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

#---------------------------
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

#-----------------------------------------
######## PLOT ########
#setwd(paste(wd, "final/", sep=""))

partitions <- c("full","plastid","nuclear")
palettes_hcl <- hcl.pals()
colors <- hcl.colors(12, palette = "SunsetDark", alpha = 1)[c(1,6,10)]
names(colors) <- partitions


# pdf(file="teste.pdf", width = 16, height = 7)
# par(mfrow=c(3,3), margin(1,1,1,1))
# all 53
pdf("../plots/treespace_partitions_all.pdf", width=6, height=6)
s.class(res_all53_part$pco$li, fac=as.factor(qualifier_all53_part), cellipse = 0, cstar=1, 
        label=NULL, pch=19,cgrid=1, cpoint = 0.5, 
        col=colors[levels(as.factor(qualifier_all53_part))])
s.chull(res_all53_part$pco$li, fac=as.factor(qualifier_all53_part), 
        add.plot = T, clabel=1.2, col=colors[levels(as.factor(qualifier_all53_part))], 
        label=names(colors[levels(as.factor(qualifier_all53_part))]), 
        optchull=1)
dev.off()

pdf("../plots/treespace_partitions_eug.pdf", width=6, height=6)
s.class(res_eug_part$pco$li, fac=as.factor(qualifier_eug_part), cellipse = 0, 
        cstar=1, label=NULL, pch=19, 
        cgrid=1, cpoint = 0.5, col=colors[levels(as.factor(qualifier_eug_part))])
s.chull(res_eug_part$pco$li, fac=as.factor(qualifier_eug_part), add.plot = T, 
        clabel=1.2, col=colors[levels(as.factor(qualifier_eug_part))], 
        label=names(colors[levels(as.factor(qualifier_eug_part))]), 
        optchull=1)
dev.off()

pdf("../plots/treespace_partitions_myr.pdf", width=6, height=6)
s.class(res_myr_part$pco$li, fac=as.factor(qualifier_myr_part), cellipse = 0, 
        cstar=1, label=NULL, pch=19, 
        cgrid=1, cpoint = 0.5, col=colors[levels(as.factor(qualifier_myr_part))])
s.chull(res_myr_part$pco$li, fac=as.factor(qualifier_myr_part), add.plot = T, 
        clabel=1.2, col=colors[levels(as.factor(qualifier_myr_part))], 
        label=names(colors[levels(as.factor(qualifier_myr_part))]), 
        optchull=1)
dev.off()

# To measure proportion of variance explained
tree.D <- vegdist(res_all53_markers$pco$li, "euclidean")
pcoa_result <- ape::pcoa(tree.D)
# Extract the eigenvalues
eigenvalues <- pcoa_result$values
# Calculate the proportion of variance explained
variance_explained <- eigenvalues / sum(eigenvalues)
# Print the proportion of variance explained by each axis
print(variance_explained)

dev.off()
#-----------------------------------------
# Myrciinae


# To measure proportion of variance explained
tree.D <- vegdist(res_myr_part$pco$li, "euclidean")
pcoa_result <- ape::pcoa(tree.D)
# Extract the eigenvalues
eigenvalues <- pcoa_result$values
# Calculate the proportion of variance explained
variance_explained <- eigenvalues / sum(eigenvalues)
# Print the proportion of variance explained by each axis
print(variance_explained[1:3,])

#-----------------------------------------
# Eugeniinae

# To measure proportion of variance explained
tree.D <- vegdist(res_eug_part$pco$li, "euclidean")
pcoa_result <- ape::pcoa(tree.D)
# Extract the eigenvalues
eigenvalues <- pcoa_result$values
# Calculate the proportion of variance explained
variance_explained <- eigenvalues / sum(eigenvalues)
# Print the proportion of variance explained by each axis
print(variance_explained)
#dev.off()

#-----------------------------
# for individual markers
s.class(res_all53_markers$pco$li, fac=as.factor(qualifier_all53_markers), cellipse = 0, 
        cstar=1, label=NULL, pch=19,cgrid=1, cpoint = 0.5, 
        col=colors[levels(as.factor(qualifier_all53_markers))])
s.chull(res_all53_markers$pco$li, fac=as.factor(qualifier_all53_markers), add.plot = T, 
        clabel=1.2, col=colors[levels(as.factor(qualifier_all53_markers))], 
        label=names(colors[levels(as.factor(qualifier_all53_markers))]), 
        optchull=1)

s.class(res_myr_markers$pco$li, fac=as.factor(qualifier_myr_markers), cellipse = 0, cstar=1, 
        label=NULL, pch=19, cgrid=1, cpoint = 0.5, 
        col=colors[levels(as.factor(qualifier_myr_markers))])
s.chull(res_myr_markers$pco$li, fac=as.factor(qualifier_myr_markers), add.plot = T, 
        clabel=1.2, col=colors[levels(as.factor(qualifier_myr_markers))], 
        label=names(colors[levels(as.factor(qualifier_myr_markers))]), 
        optchull=1)
s.class(res_eug_markers$pco$li, fac=as.factor(qualifier_eug_markers), cellipse = 0, cstar=1, 
        label=NULL, pch=19, cgrid=1, cpoint = 0.5, 
        col=colors[levels(as.factor(qualifier_eug_markers))])
s.chull(res_eug_markers$pco$li, fac=as.factor(qualifier_eug_markers), add.plot = T, 
        clabel=1.2, col=colors[levels(as.factor(qualifier_eug_markers))], 
        label=names(colors[levels(as.factor(qualifier_eug_markers))]), 
        optchull=1)
