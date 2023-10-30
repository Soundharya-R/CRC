library(viridis)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(robustbase)
library(rstatix)
library(data.table)
library(ggsci)
library(ggpubr)
library(ggnewscale)


setwd("path to folder containing the count matrix file")
     counts<-fread("./single-cell-count-matrix-file.txt", header = T) #column names = cell barcodes
     counts<-as.data.frame(counts); gc(reset = T)
     rownames(counts)<-counts$V1; counts<-counts[,-1]
     rownames_temp <- rownames(counts)
     counts <- apply(counts, 2, scale)
     rownames(counts) <- rownames_temp
     counts <- as.data.frame(counts[, apply(counts, 2, function(x) !any(is.na(x)))] + 1)
     entropy_table <- data.frame(row.names = colnames(counts))
     pathway_gene_list <- data.frame(fread("./tab-delimited-gmt-containg-relevant-signatures", fill = T), row.names = 1)
     for (i in 1:nrow(pathway_gene_list))
     {
       pathway_genes <- as.character(pathway_gene_list[i,])
       pathway_genes <- pathway_genes[which(pathway_genes != "")]
       count_matrix_forrun <- na.omit(counts[pathway_genes,])
       probs <- t(t(count_matrix_forrun)/colSums(count_matrix_forrun))
       entropy_cell <- -apply(probs * log(probs + 1e-10)/log(nrow(count_matrix_forrun)), 2, sum)
       entropy_table[,i] <- entropy_cell
       names(entropy_table)[i] <- rownames(pathway_gene_list)[i]
     }
     
     rm(counts, count_matrix_forrun, probs, pathway_gene_list, entropy_cell, pathway_genes, i); gc(reset = T)
     colnames(entropy_table)[c(1:nrow(entropy_table)] <- c("names of pathways chosen")    

     write.table(entropy_plot, "./entropy.txt", row.names = T, sep = '\t', quote = F)
