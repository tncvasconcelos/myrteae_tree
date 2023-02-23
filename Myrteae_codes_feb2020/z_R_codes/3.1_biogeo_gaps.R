#### Biogeographical COVERAGE ####

####################################


setwd("~/Desktop/Myrteae_codes/3_species_distribution/")

# rm(list=ls())

library(stringr)
library(dismo) 
library(CoordinateCleaner)
library(dplyr)
library(monographaR)
library(schoolmath)
library(RColorBrewer)
data(wrld_simpl)

# TV: primeiro puxa os pontos das spp da filo

###### phylo #####
phylo <- as.character(read.csv("phylo_list.csv")[,1])
phylo <- gsub("_.*", "", phylo)
bla <- c("aff.","sp.","sp+[[:digit:]]","\\<sp\\>","_", "cf.", "var.")
grep(paste(bla, collapse="|"), phylo) -> bla
phylo[-bla]->phylo
unique(phylo)->phylo

### GBIF ###
start_time <- Sys.time()
datalist_phylo = list()
for (z in 1:length(phylo)){
  tt <- strsplit(phylo,".",fixed = TRUE)[[z]]
  tt[1] -> x
  tt[2] -> y
  print(c(x,y)) # to check names
  dat <- gbif(x, y) # download points from gbif
  datalist_phylo[[z]] <- dat 
}
big_data_phylo <- data.table::rbindlist(datalist_phylo, fill=T)

end_time <- Sys.time()
end_time - start_time # Time difference of 1.776503 hours

big_data_phylo <- as.data.frame(big_data_phylo[,c("species","lon","lat")])

# TV: aqui puxa os pontos das spp da lista do VPA ("Vascular Plants of America"),
# eu acho que o ideal na verdade seria puxar uma lista feita pelo WCSPF do Kew, mas
# nao consegui um jeito facil de fazer isso... Essa lista do VPA tem um monte de nome bizarro,
# mas eu acredito que o GBIF corrija para os sinonimos na hora de puxar (mas não tenho certeza)
# enfim, tem que ver isso... o ideal seria conseguir uma lista melhor de todas as especies
# Neotropicais, eu acho que eé possivel conseguir isso pelo Kew (perguntar pra Eve)

### VPA ####

vpa <- as.character(read.csv("AdvancedSearchResults.csv")[,2])
first<- word(vpa, 1)
second<- word(vpa, 2)
vpa <- paste(first, second)

# some data cleaning
bla <- c("NA")
grep(paste(bla, collapse="|"), vpa) -> bla
vpa[-bla] -> vpa 
unique(vpa) -> vpa

#
vpa_cleaned<-c()
for (z in 1:length(vpa)){
  tt <- strsplit(vpa," ",fixed = TRUE)[[z]]
  tt[1] -> x
  tt[2] -> y
  cap <- substr(y, 1, 1)
  if(is.na(x)){
    print(c(x,y)) # 
    print("deu xabu - deletando")
  } 
  else{
    if(str_detect(cap, "^[:upper:]+$")){
      print(c(x,y)) # 
      print("deu xabu - deletando")
    }
    else {
      v0 <- paste(x, y)
      vpa_cleaned <- c(v0, vpa_cleaned)
    }
  }
}

### GBIF ###
start_time <- Sys.time()
datalist_VPA = list()
for (z in 1:length(vpa_cleaned2)){
  tt <- strsplit(vpa_cleaned2," ",fixed = TRUE)[[z]]
  tt[1] -> x
  tt[2] -> y
  print(c(x,y)) # to check names
  dat <- gbif(x, y) # download points from gbif
  datalist_VPA[[z]] <- dat 
}

big_data_vpa <- data.table::rbindlist(datalist_VPA, fill=T)

end_time <- Sys.time()
end_time - start_time # Time difference of 3.997529 hours

big_data_vpa <- as.data.frame(big_data_vpa[,c("species","lon","lat")])

# TV: aqui tem um pouco de data cleaning 

#################################
# first round of data cleaning #
################################

big_data_vpa <- read.csv("cleaned_points_vpa.csv")[,2:4] 
big_data_phylo <- read.csv("cleaned_points_phylo.csv")[,2:4]

# vpa
points <- subset(big_data_vpa, !is.na(lon) & !is.na(lat)) 
points <- cc_cen(points, lon="lon", lat="lat", species ="species")
points <- cc_cap(points, lon="lon", lat="lat", species ="species")
points <- cc_dupl(points, lon="lon", lat="lat", species ="species")
points <- cc_equ(points, lon="lon", lat="lat")
points <- cc_inst(points, lon="lon", lat="lat", species ="species")
points <- cc_val(points, lon="lon", lat="lat")
points <- cc_sea(points, lon="lon", lat="lat")
# restricting to the Neotropics 
points %>% filter(lon < -20, lon >- 130, lat < 40, lat >- 40) -> cleaned_points_vpa

# phylo
points <- subset(big_data_phylo, !is.na(lon) & !is.na(lat)) 
points <- cc_cen(points, lon="lon", lat="lat", species ="species")
points <- cc_cap(points, lon="lon", lat="lat", species ="species")
points <- cc_dupl(points, lon="lon", lat="lat", species ="species")
points <- cc_equ(points, lon="lon", lat="lat")
points <- cc_inst(points, lon="lon", lat="lat", species ="species")
points <- cc_val(points, lon="lon", lat="lat")
points <- cc_sea(points, lon="lon", lat="lat")
# restricting to the Neotropics 
points %>% filter(lon < -20, lon >- 130, lat < 40, lat >- 40) -> cleaned_points_phylo


#### combining tables for second round ####
col1 <- rep("vpa", length(big_data_vpa$Species))
big_data_vpa <- cbind(col1, big_data_vpa)
col1 <- rep("phylo", length(big_data_phylo$Species))
big_data_phylo <- cbind(col1, big_data_phylo)
myrteae_total <- rbind(big_data_vpa, big_data_phylo)

####
# TV: Aqui filtra uns pseudocentroides

# 3.3. Filtering centroids #
# centroids for: 
# Brazil, Bolivia, Peru, Colombia, Venezuela, Chile... 
# (should probably add others..)
loncent <- c(-52.87310,-64.435,-74.14,-72.8667,-65.9119,-70.8906,
             -64.9208, -79.016069, -78.7520357, -58.1689,	-56.0180680,
             -55.625, -52.9708, -58.7017,	-83.9519, -90.1792, -89.35083,
             -89.93333)
latcent <- c(-10.833900,-16.7261,-9.1839,3.8811,7.0758,-35.8156,
             -35.3869, 21.6229002, -1.4238195, -23.2025, -32.7995198,
             4.1006, 3.8558, 4.7364, 10.0114, 15.7422, 16.08111,
             13.81667)

# flaging and removing centroids
remove_centroid <- c()
for(i in 1:length(myrteae_total$Species)) {
  for(u in 1:length(latcent))
    if(myrteae_total[i,3] == loncent[u] && myrteae_total[i,4] == latcent[u]){
      print(paste(i, "is centroid!", sep=" "))
      remove_centroid <- c(remove_centroid, i)
    }
}
myrteae_total <- myrteae_total[-remove_centroid, ]

# TV: aqui tira tudo que é ponto que nao tem casa decimal (innacurate)

# 3.4. Removing inaccurate points #
# removing points without decimal cases #
keep_accurate <- c()
for(i in 1:length(myrteae_total$Latitude)) {
  if(is.decimal(myrteae_total[i,3]) && is.decimal(myrteae_total[i,4])){
    print(paste(i, "is accurate!", sep=" "))
    keep_accurate <- c(keep_accurate, i)
  }
}
myrteae_total <- myrteae_total[keep_accurate, ]

# TV: aqui tira outliers

# 3.5. Removing outliers in distribution #
names(table(myrteae_total$Species))->species
cleaned_point <- data.frame()
for(u in 1:length(species)){
  sp0 <- myrteae_total[myrteae_total$Species==species[u],]
  out_lat <- boxplot.stats(sp0$Latitude)$out
  out_lon <- boxplot.stats(sp0$Longitude)$out
  sp <- sp0[ ! sp0$Latitude %in% out_lat, ]
  sp <- sp[ ! sp$Longitude %in% out_lon, ]
  cleaned_point <- rbind(cleaned_point, sp)
}

###################

# writing tables
vpa <- cleaned_point[cleaned_point$col1=="vpa",]
vpa <- vpa[,-1]
phylo <- cleaned_point[cleaned_point$col1=="phylo",]
phylo <- phylo[,-1]

write.csv(vpa, file="cleaned_points_vpa.csv") 
write.csv(phylo, file="cleaned_points_phylo.csv")

###################################
######## contrasting maps #########
###################################

vpa<- read.csv("cleaned_points_vpa.csv")[,2:4] 
phylo<- read.csv("cleaned_points_phylo.csv")[,2:4] 

###

# checar a distribuição
nameCol <- c("sp", "Longitude","Latitude")
colnames(vpa) <- nameCol
xx1<-mapDiversity(vpa, resolution = 1, plot = TRUE, plot.with.grid = F, legend = T, export = F) #

colnames(phylo) <- nameCol
xx2<-mapDiversity(phylo, resolution = 1, plot = F, plot.with.grid = F, legend = T, export = F) #

# TV: Nao sei se eé a melhor maneira de plotar esses mapas e tambem naão troquei as cores...

#####

pdf(file="teste.pdf", width=12, height=10)
par(mfrow=c(2,2))

plot(crop(xx1, extent(-130, -20, -40, 40)),
     main = "Total Neotropical Myrtaceae (vpa)", zlim=c(1,220))
plot(wrld_simpl, add=T,  border="black", lwd=0.5)

plot(crop(xx2, extent(-130, -20, -40, 40)),
     main = "Total Neotropical Myrtaceae (phylo)", zlim=c(1,220))
plot(wrld_simpl, add=T,  border="black", lwd=0.5)

all <- lapply(list(xx1, xx2), crop, extent(-130, -20, -40, 40))
res1 <- Reduce("-", all, accumulate = TRUE)
l1<-length(res1)
plot(res1[[l1]], main = "Total Neotropical Myrtaceae (difference)", zlim=c(1,90))
plot(wrld_simpl, add=T)

all <- lapply(list(xx1, xx2), crop, extent(-130, -20, -40, 40))
res2 <- Reduce("/", all, accumulate = TRUE)
l1<-length(res2)
plot(res2[[l1]], main = "Total Neotropical Myrtaceae (proportion)", zlim=c(0,5))
plot(wrld_simpl, add=T, lwd=0.5)

dev.off()

####


