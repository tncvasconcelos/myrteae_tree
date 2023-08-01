#rm(list=ls())
# setwd("~/Desktop/myrteae_traits")
library(RColorBrewer)
library(maptools)
library(raster)
data("wrld_simpl")

#----------------
dist_sample <- read.table("WCVP/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("WCVP/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")

myrteae_genera <- read.csv("neotropical_myrteae_genera_WCVP.csv")
tree <- read.tree("mmc_target_common_June27_pruned_no_out.tre")

pimentinae <- c("Acca","Amomyrtella","Amomyrtus","Campomanesia","Curitiba",
                "Feijoa","Legrandia","Mosiera","Myrrhinium","Pimenta","Psidium")
ugniinae <- c("Ugni","Myrteola","Lenwebbia","Lophomyrtus")
myrciinae <- c("Myrcia")
myrtinae <- c("Accara","Chamguava","Myrtus","Calycolpus")
pliniinae <- c("Plinia","Myrciaria","Neomitranthes","Siphoneugena","Algrizea")
eugeninae <- c("Calycorectes","Eugenia","Myrcianthes")
luminae <- c("Luma","Myrceugenia","Temu")
blepharocalycinae <- "Blepharocalyx"

all_vars <- subset(all_vars, all_vars$family %in% "Myrtaceae")
all_vars <- subset(all_vars, all_vars$taxon_status %in% "Accepted")
all_vars <- subset(all_vars, all_vars$species != "")
all_vars <- subset(all_vars, all_vars$genus %in% myrteae_genera$neotropical_myrteae_genera)
all_species <- unique(all_vars$taxon_name)

all_subtribes <- list(pimentinae,ugniinae,myrciinae,myrtinae,pliniinae,eugeninae,luminae,blepharocalycinae)
names(all_subtribes) <- c("pimentinae","ugniinae","myrciinae","myrtinae","pliniinae","eugeninae","luminae","blepharocalycinae")
coverage <- matrix(nrow=length(all_subtribes), ncol=3)
for(i in 1:length(all_subtribes)) {
  one_subtribe <- all_subtribes[[i]]
  total_length <- length(grep(paste(one_subtribe, collapse="|"), all_species))
  sampled <- length(grep(paste(one_subtribe, collapse="|"), tree$tip.label))
  prop_coverage <- sampled/total_length
  coverage[i,c(1:3)] <- c(total_length, sampled, prop_coverage)
}
rownames(coverage) <- names(all_subtribes)
colnames(coverage) <- c("total_richness","sampled","coverage")

write.csv(coverage,file="sampling_coverage.csv",row.names=T)

sum(coverage[,2]) / sum(coverage[,1]) 

###
# Plot for figure 
tree_data <- read.csv("tree_data.csv")
#write.csv(tree_data, "tree_data.csv", row.names=F)

library(ggplot2)

tree_data <- subset(tree_data, tree_data$Group != "Myrteae")

plot1 <- ggplot(tree_data, aes(x = reorder(Group, -order), y = age_mean)) + 
  geom_point(size=2) +
  scale_y_reverse() +
  coord_flip(ylim = c(45, 0)) +
  theme_bw() +
  geom_errorbar(aes(ymin = age_max, ymax = age_min), width =  0.75) 

plot2 <- ggplot() + 
  geom_bar(aes(y = taxonomic_coverage, x = reorder(Group, -order)), data = tree_data,
                          stat="identity") +
  coord_flip(ylim = c(0, 1))  + theme_bw() 

pdf("age_and_prop.pdf")
plot1
plot2
dev.off()

