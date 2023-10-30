library(CMScaller)
library(data.table)



GSE="GSE_ID" #submit GSE Id for dataset to be processed
setwd("path_to_working folder")
path_to_results<-paste0("path/", GSE, "/", GSE)
count_matrix <- fread(paste0(path_to_results, "_countmatrix.txt"), sep = "\t") #Assuming you have saved the counts as GSE_countmatrix.txt

Subtyping = function(count_matrix){
  #load count matrix
  count_matrix <- as.matrix(KeepMaxGene(count_matrix)) # to keep maximally expressed genes only
  png(paste0(path_to_results, "_CMSClassifierHeatmap.png") , width = 10, height = 10, units = "cm", res = 500)
  #CMS classifying
  CMSclassifier <- CMScaller(count_matrix, rowNames = "symbol", RNAseq = F)
  dev.off()
  CMSclassifier <- cbind(rownames(CMSclassifier), CMSclassifier)
  return(CMSclassifier)
  
}
  
fwrite(CMSclassifier, paste0(path_to_results, "_results.txt"), sep = "\t")

