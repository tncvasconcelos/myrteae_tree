# rm(list=ls())
setwd("~/Desktop/myrteae_tree")
library(ape)

#-----------------------------
# Updating names
raxml_tree <- read.tree("trees_final/RAxML_result_June25.tre")
beast_tree <- read.tree("trees_final/mmc_target_common_June25.tre")

beast_tree$tip.label <- gsub(".*-","",beast_tree$tip.label) #removing label annotations

tips_to_update <- read.csv("name_updates_may23.csv")

# Arguments are the tree and a data.frame with one column with old names and a column with new names.
# Tips that are not to be changed should be left in blank ("") on second column.
UpdateTips<- function(one_tree, tips_to_update){
  tree_tips <- one_tree$tip.label
  tips_to_update <- subset(tips_to_update, tips_to_update[,2]!="")
  for(tip_index in 1:length(tree_tips)) {
    if(tree_tips[tip_index] %in% tips_to_update[,1]) {
      tree_tips[tip_index] <- tips_to_update[,2][which(tips_to_update[,1] == tree_tips[tip_index])]
    }
  }  
  one_tree$tip.label <- tree_tips
  return(one_tree)
}

raxml_tree_updated <- UpdateTips(raxml_tree, tips_to_update)
beast_tree_updated <- UpdateTips(beast_tree, tips_to_update)

write.tree(raxml_tree_updated, file="trees_final/RAxML_result_June27.tre")
write.tree(beast_tree_updated, file="trees_final/mmc_target_common_June27.tre")

# names should be updated in other lists too:
tips_to_remove <- read.csv("names_to_remove_1_sp_only.csv")
tips_no_voucher <- read.csv("names_no_vouchers.csv")
tips_to_update <- subset(tips_to_update, tips_to_update[,2]!="")
sort(tips_to_update$old_name)

tips_to_change <- tips_to_remove[,1][which(tips_to_remove[,1] %in% tips_to_update[,1])] 
if(length(tips_to_change)>0) {
  for(i in 1:length(tips_to_change)) {
    tips_to_remove[,1][which(tips_to_remove[,1]==tips_to_change[i])] <- tips_to_update[,2][which(tips_to_update[,1]==tips_to_change[i])]
  } 
}

tips_to_change <- tips_no_voucher[,1][which(tips_no_voucher[,1] %in% tips_to_update[,1])] 
if(length(tips_to_change)>0) {
   for(i in 1:length(tips_to_change)) {
     tips_no_voucher[,1][which(tips_no_voucher[,1]==tips_to_change[i])] <- tips_to_update[,2][which(tips_to_update[,1]==tips_to_change[i])]
   } 
}

write.csv(tips_to_remove, "names_to_remove_1_sp_only.csv", row.names = F)
write.csv(tips_no_voucher, "names_no_vouchers.csv", row.names = F)

#-----------------------------
# Keeping only one tip per species
# Load most up-to-date trees:
raxml_tree <- read.tree("trees_final/RAxML_result_June27.tre")
beast_tree <- read.tree("trees_final/mmc_target_common_June27.tre")
tips_to_remove <- read.csv("names_to_remove_1_sp_only.csv")
tips_no_voucher <- read.csv("names_no_vouchers.csv")

raxml_tree_pruned <- drop.tip(raxml_tree, tips_to_remove[,1])
beast_tree_pruned <- drop.tip(beast_tree, tips_to_remove[,1])

raxml_tree_pruned <- keep.tip(raxml_tree_pruned, tips_no_voucher[,1])
beast_tree_pruned <- keep.tip(beast_tree_pruned, tips_no_voucher[,1])

#-----------------------------
# Removing vouchers
# raxml_tree <- read.tree("trees_final/RAxML_result_June27_pruned.tre")
# beast_tree <- read.tree("trees_final/mmc_target_common_June27_pruned.tre")

raxml_tree_pruned <- UpdateTips(raxml_tree_pruned, tips_no_voucher)
beast_tree_pruned <- UpdateTips(beast_tree_pruned, tips_no_voucher)

raxml_tree_pruned <- ladderize(raxml_tree_pruned)
beast_tree_pruned <- ladderize(beast_tree_pruned)

write.tree(raxml_tree_pruned, file="trees_final/RAxML_result_June27_pruned.tre")
write.tree(beast_tree_pruned, file="trees_final/mmc_target_common_June27_pruned.tre")

# Eugenia.involucrata_Genome.OP650216
#write.csv(beast_tree_pruned$tip.label, file="prune_outgroups.csv",row.names = F)
#-----------------------------
# Removing ougroups
raxml_tree <- read.tree("trees_final/RAxML_result_June27_pruned.tre")
beast_tree <- read.tree("trees_final/mmc_target_common_June27_pruned.tre")

outgroups <- read.csv("prune_outgroups.csv")

raxml_tree_no_out <- drop.tip(raxml_tree, outgroups[,1][which(outgroups[,2]!="")])
beast_tree_no_out <- drop.tip(beast_tree, outgroups[,1][which(outgroups[,2]!="")])

write.tree(raxml_tree_no_out, file="trees_final/RAxML_result_June27_pruned_no_out.tre")
write.tree(beast_tree_no_out, file="trees_final/mmc_target_common_June27_pruned_no_out.tre")


