remove(list = ls())
graphics.off()  

library(phyloseq)
library(vegan)
library(picante)

setwd("pathway")
data <- read.csv(file = "sample_asv.csv", header = T)

otu_table_data <- data[, -ncol(data)]
rownames(otu_table_data) <- otu_table_data$ID
otu_table_data <- otu_table_data[, -1]
otu_matrix <- as.matrix(otu_table_data)
taxonomy_data <- data[, c(1, ncol(data))]
colnames(taxonomy_data) <- c("OTU", "taxonomy")
taxonomy_split <- strsplit(as.character(taxonomy_data$taxonomy), ";")
taxonomy_matrix <- do.call(rbind, lapply(taxonomy_split, function(x) {
  x <- gsub("^\\s+|\\s+$", "", x) 
  x <- c(x, rep(NA, 7 - length(x))) 
}))
colnames(taxonomy_matrix) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(taxonomy_matrix) <- taxonomy_data$OTU

phy_tree <- read.tree("rooted_tree.tre") 
otu_table_obj <- otu_table(otu_matrix, taxa_are_rows = TRUE)
taxonomy_table_obj <- tax_table(as.matrix(taxonomy_matrix))
physeq <- phyloseq(otu_table_obj, taxonomy_table_obj, phy_tree(phy_tree))

otu_df <- as.data.frame(otu_table(physeq))
head(otu_df)

min_depth <- min(sample_sums(physeq)) 
ps_rarefied <- rarefy_even_depth(physeq, sample.size = min_depth, rngseed = 123)
alpha_diversity <- estimate_richness(ps_rarefied, measures = c("Shannon", "Simpson", "Chao1","Observed"))
print(alpha_diversity)
pd_values <- pd(t(otu_table(ps_rarefied)), phy_tree(ps_rarefied), include.root = TRUE)
alpha_diversity$PD_Whole_Tree <- pd_values$PD
print(alpha_diversity)

metadata <- read.csv(file = "sample_information.csv", header = TRUE, row.names = 1)
alpha_diversity$Sample <- rownames(alpha_diversity)
combined_data <- merge(alpha_diversity, metadata, by = "Sample")
write.csv(combined_data, file = "combined_data.csv", row.names = FALSE)


