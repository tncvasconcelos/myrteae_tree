
#### TAXONOMIC COVERAGE ####

### 
setwd("~/Desktop/Myrteae_codes/3_taxonomic_gaps/")

# TV: nao fiz muita coisa nesse, mas se conseguir a lista com os valores certos eu facço o resto
taxa <- read.csv("taxonomic.csv", h=F)[-1,]
taxa = t(taxa[,2:4])
total <- as.numeric(taxa[2,])
phylo <- as.numeric(taxa[3,])
unsampled <- total - phylo
proportion <- phylo / total
taxa <- rbind(taxa, unsampled, proportion)
taxa <- taxa[,order(as.numeric(taxa[2, ]), decreasing=T)]
colnames(taxa) <- taxa[1,]
barplot(as.matrix(taxa[-c(1,2,5),]), cex.names = 0.4, las = 2)

