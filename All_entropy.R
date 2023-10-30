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
source("count_to_TPM_function.R")
source("UsefulFunctions.R")

setwd("D:/ProjectWork/CRC/Codes/GSE132257/")
#     counts<-fread("./GSE144735_countmatrix.txt", header = T) 
#     counts<-as.data.frame(counts); gc(reset = T)
#     rownames(counts)<-counts$V1
#     counts<-counts[,-1]
#     counts[1:5,1:5]
#     # count_matrix <- countToTpm(count_matrix) 
#     rownames_temp <- rownames(counts)
#     counts <- apply(counts, 2, scale)
#     rownames(counts) <- rownames_temp
#     counts <- as.data.frame(counts[, apply(counts, 2, function(x) !any(is.na(x)))] + 1)
#     entropy_table <- data.frame(row.names = colnames(counts))
#     pathway_gene_list <- data.frame(fread("../Entropy_signatures.txt", fill = T), row.names = 1)
#     # pathway_gene_list <- data.frame(fread("signatures_for_entropy2.txt", fill = T), row.names = 1)
#     for (i in 1:nrow(pathway_gene_list))
#     {
#       pathway_genes <- as.character(pathway_gene_list[i,])
#       pathway_genes <- pathway_genes[which(pathway_genes != "")]
#       count_matrix_forrun <- na.omit(counts[pathway_genes,])
#       probs <- t(t(count_matrix_forrun)/colSums(count_matrix_forrun))
#       entropy_cell <- -apply(probs * log(probs + 1e-10)/log(nrow(count_matrix_forrun)), 2, sum)
#       entropy_table[,i] <- entropy_cell
#       names(entropy_table)[i] <- rownames(pathway_gene_list)[i]
#     }
#     
#     rm(counts, count_matrix_forrun, probs, pathway_gene_list, entropy_cell, pathway_genes, i); gc(reset = T)
#     colnames(entropy_table)[c(1:7)] <- c("FAO","Glycolysis","OXPHOS","hEMT","Epi","Mes","PDL1")
#     
#     CMS_data <- read.table("./CMS.txt")
#     entropy_table <- cbind(CMS_data, entropy_table)
#     entropy_plot <- entropy_table[!entropy_table$CMScaller == "CMSU",]
#     entropy_plot <- entropy_plot[, -c(2,3)]
#     write.table(entropy_plot, "./entropy.txt", row.names = T, sep = '\t', quote = F)
# #___________    
#     entropy_plot <- entropy_plot[, -c(1,2,3)]  # entropy_plot <- entropy_plot[, -c(2,3)]
#     # entropy_plot$Sample <- rownames(entropy_plot)
#     entropy_plot$sample <- rownames(entropy_plot)
#     entropy_plot <- melt(entropy_plot, id.vars = "sample")
#     # entropy_plot <- entropy_plot[entropy_plot$variable != "Sample"]
#     entropy_plot$scaled <- scale(entropy_plot$value); names(entropy_plot)[2] <- "signature"
#     stat.test <- entropy_plot %>% t_test(value ~ signature) 
#     stat.test <- stat.test %>% add_xy_position(x = "CMScaller", step.increase = 0.5)
#     stat.test$y.position <- stat.test$y.position + 0.0004
#     stat.test
#     ggplot(entropy_plot, aes(y = value, x = signature)) +
#       geom_boxplot(aes(fill = signature),show.legend = F,  alpha = 0.9, outlier.shape = NA) + scale_fill_lancet() + 
#       # geom_jitter(shape=3, position = position_jitter(0.6), 
#                   # color = "black", fill = "white", size = 0.01) +
#       theme_classic() +  
#       labs(title = "", y = "Entropy", x = "") + 
#       # scale_y_continuous(expand = c(0,0)) +
#       scale_x_discrete(labels = c("FAO", "Glycolysis", "OXPHOS", "hEMT", "Epi", "Mes","PD-L1")) +
#       theme(panel.grid = element_blank(), 
#             panel.border = element_rect(colour = "black", fill = NA, size = 1.22),
#             axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10), 
#             legend.title = element_blank())+
#       coord_cartesian(ylim = c(0.83, 0.985))
#       stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, size = 2)
#     
#     ggsave("./all_entropy.pdf")
    
#---------------------------------------------------------------------------------------------------------
    MeanBaseCreation <- function(entropy_plot){
      entropy_line <- do.call(data.frame, aggregate(. ~ CMScaller, entropy_plot, function(x) c("mean" = mean(x), "sd" = sd(x))))
      return(entropy_line)
    }
    setwd("/home/csb/Desktop/Manas/CRC/GSE132465/")
    entropy_plot <- read.table("./entropy.txt", row.names = 1, header = T)
    entropy_plot <- entropy_plot[c(1, 6, 7)]
    entropy_line <- MeanBaseCreation(entropy_plot)
    entropy_line_mean <- entropy_line[,c(1,2,4)]; entropy_line_mean <- melt(entropy_line_mean, id.vars = "CMScaller"); 
    names(entropy_line_mean)[2:3] <- c("variable1","mean")
    entropy_line_sd <- entropy_line[,c(1,3,5)]; entropy_line_sd <- melt(entropy_line_sd, id.vars = "CMScaller")
    entropy_line <- cbind(entropy_line_mean, entropy_line_sd[,-1])
    ggplot(entropy_line, aes_string(x = "mean", y = "value", label = "variable1")) + geom_path(aes(group = variable1)) +
      geom_point(aes(color = CMScaller)) + scale_color_jama() + theme_classic() +  labs(x = "Mean", y = "Standard deviation") +
      theme(legend.title = element_blank(), 
            panel.grid = element_blank(), 
            panel.border = element_rect(colour = "black", fill = NA, size = 1.22),) +  
      annotate("text", x = 0.935, y = 0.003, label= "E entropy") +
      annotate("text", x = 0.916, y = 0.0072, label= "M entropy") 
    ggsave(paste0("./","EM_mean_sd",".png"), dpi = 1000, width = 9.8, height = 6.5, units = "cm")
    
#------------------------------   
    setwd("F:/CRC/Codes/GSE132465/")
    metrics <- c("FAO", "Glycolysis", "OXPHOS", "Epi", "Mes","PD-L1")
    for (i in 1:length(metrics)) 
    {
    entropy_plot <- read.table("./entropy.txt", row.names = 1, header = T)
    names(entropy_plot)[8] <- "PD-L1"
    # entropy_plot$sample <- rownames(entropy_plot)
    entropy_plot <- melt(entropy_plot, id.vars = "CMScaller")
    # entropy_plot <- entropy_plot[entropy_plot$variable != "Sample"]
    # entropy_plot$scaled <- scale(entropy_plot$value)
      entropy_plot <- entropy_plot[entropy_plot$variable == metrics[i],]
      entropy_plot <- entropy_plot[,-2]
      colnames(entropy_plot)[2] <- "value"
      stat.test <- entropy_plot %>%  t_test(value ~ CMScaller) 
      stat.test <- stat.test %>% add_xy_position(x = "CMScaller", scales = "free_y")
      # stat.test$y.position <- stat.test$y.position + 0.0004
      stat.test
      plot <- ggplot(entropy_plot, aes(y = value, x = CMScaller)) +
        geom_boxplot(aes(fill = CMScaller), show.legend = F,  alpha = 1, outlier.shape = NA) + 
        # scale_fill_lancet() + 
        scale_fill_igv() + 
        # geom_jitter(shape=3, position = position_jitter(0.6),
        # color = "black", fill = "white", size = 0.01) +
        theme_classic() +  
        labs(title = "", y = "Entropy", x = "") +
        # scale_y_continuous(expand = c(0.07,0.1)) +
        # scale_x_discrete(labels = c("FAO", "Glycolysis", "OXPHOS", "hEMT", "Epi", "Mes","PD-L1")) +
        theme(panel.grid = element_blank(), 
              panel.border = element_rect(colour = "black", fill = NA, size = 1.1),
              axis.text.x = element_text(angle = 0, hjust = 0.44, vjust = 0, size = 12, face = "bold", colour = "black"), 
              legend.title = element_blank(),
              plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + ggtitle(metrics[i])
      # coord_cartesian(ylim = c(0.83, 0.985))
      plot <- plot + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, size = 5, hide.ns = F, step.increase = 0.096)
      plot
      ggsave(paste0("./Replotted/",metrics[i],"_entropy.pdf"))
      ggsave(paste0("./Replotted/",metrics[i],"_entropy.png"), height = 5.44, width = 4.38)
    }
    
    