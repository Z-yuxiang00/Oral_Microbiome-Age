rm(list=ls())
graphics.off() 
setwd("pathway")
library(tidyverse)
library(microeco) 
library(magrittr)
library(ggtree)
library(dplyr)

#---------
asv_data <- read.csv('sample_asv.csv', row.names = 1)
asv_data_filtered <- asv_data[rowSums(asv_data > 0) > 0.01 * ncol(asv_data), ]
asv_data_relative <- sweep(asv_data_filtered, 2, colSums(asv_data_filtered), FUN = "/")
asv_data_relative <- sweep(asv_data, 2, colSums(asv_data), FUN = "/")
sample_table <- read.csv('sample_information.csv', row.names = 1)
tax_table <- read.csv('tax_table.csv', row.names = 1)
feature_table <- asv_data_relative
tax_table %<>% tidy_taxonomy 
#---------
dataset <- microtable$new(sample_table = sample_table,
                          otu_table = feature_table, 
                          tax_table = tax_table)

set.seed(8)
lefse <- trans_diff$new(dataset = dataset, 
                        method = "lefse", 
                        group = "Age_group", 
                        alpha = 0.05, 
                        lefse_subgroup = NULL, 
                        taxa_level = "all", 
                        lefse_min_subsam = 3
) 

#--------
lefse_results <- lefse$res_diff
filtered_features <- lefse_results[lefse_results$LDA > 2, ]
write.csv(filtered_features,file = "LDA2_filtered_features.csv")
use_number <- which(lefse$res_diff$LDA > 2)
bar_plot <- lefse$plot_diff_bar(use_number = use_number, 
                    width = 0.8, 
                    group_order = c("Young", "Middle-aged", "Old"))
bar_plot + 
    ggplot2::scale_fill_manual(
        values = c("Young" = "#F9766E", "Middle-aged" = "#429A60", "Old" = "#6DBDF2")
    )

#----------
lefse_results <- lefse$res_diff
filtered_features <- lefse_results[lefse_results$LDA > 2, ]
print(filtered_features)
print(nrow(filtered_features))

cladogram <- lefse$plot_diff_cladogram(use_taxa_num =1000, 
                                       use_feature_num = nrow(filtered_features), 
                                       clade_label_level = 5,
                                       group_order = c("Young", "Middle-aged", "Old"),
                                       color = c( "Middle-aged" = "#429A60", "Old" = "#6DBDF2"),
                                       filter_taxa = 0, 
                                       only_select_show = F 
                                       )

cladogram