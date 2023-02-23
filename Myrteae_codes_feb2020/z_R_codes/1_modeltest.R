###### MODEL TEST ####

library(phangorn)

setwd("~/Desktop/Myrteae_codes/1_subsets/")

list.files(pattern=".txt") -> files
sapply(files, read.phyDat, simplify=F) -> align
sub(".txt", "", files) -> labs

results <- c()
for(i in 1:length(align)){
  Sys.time() -> start_time
  
  m0 <- modelTest(align[[i]], 
                  model = c("JC", "F81", "K80", "HKY", "SYM", "GTR"), 
                  G = TRUE, I = TRUE, multicore=TRUE, mc.cores=2)  
  r0 <- paste(labs[i], m0$Model[which.max(m0$AICw)], round(max(m0$AICw),3))
  results <- c(results, r0)
  
  Sys.time() -> end_time
  print(c(labs[i], "done!"))
  print(end_time-start_time)
}


write.csv(results, file="modeltest_results.txt") # TV: eu deixei mostrando soó o melhor modelo, 
# mas tem como puxar todos os outros na tabela



