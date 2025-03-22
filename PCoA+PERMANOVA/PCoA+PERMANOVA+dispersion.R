remove(list = ls())
graphics.off()
setwd("pathway")

#---------------------------------------------------------------------------------------------------
library(vegan)
bc_1 <- read.csv(file = "sample_asv.csv", header = T, row.names = 1)
bc_1 <- na.omit(bc_1) 
bc_1 <- sweep(bc_1, 2, colSums(bc_1), "/") 
bc_1 <- t(bc_1)
# bc_1 <- decostand(bc_1, method = "pa") 
bc_2 <- read.csv(file = "sample_information.csv", header = T, row.names = 1)
dune.env<-as.matrix(bc_2)
dune<-as.matrix(bc_1)

#--------------------------------------------------------------------------------------------------
library(vegan)
set.seed(8)
dune_dist <- vegdist(dune, method="bray", binary=F)
dune_pcoa <- cmdscale(dune_dist, k=3, eig=T)
dune_pcoa_points <- as.data.frame(dune_pcoa$points)
sum_eig <- sum(dune_pcoa$eig)
eig_percent <- round(dune_pcoa$eig/sum_eig*100,1)
colnames(dune_pcoa_points) <- paste0("PCoA", 1:3)
dune_pcoa_result <- cbind(dune_pcoa_points, dune.env)
head(dune_pcoa_result)

#-------------------------------------------------------------------------------

library(ggsci)
library(extrafont)
library(ggplot2)

color_palette <- c("Young" = "#a0bdd5", "Middle-aged" ="#ca8ba8" , "Old" = "#8583a9")
shape_palette <- c("Young" = 17, "Middle-aged" = 16, "Old" = 15)  

p_age <-ggplot(dune_pcoa_result, aes(x=PCoA1, y=PCoA2, color=factor(Age_group), shape=factor(Age_group))) + 
    geom_point(size=3.5, alpha=1) + 
    stat_ellipse(level=0.95, linetype=1, linewidth=0.8, aes(fill = factor(Age_group)),geom = "polygon", alpha = 0.3) +  
    scale_color_manual(values=color_palette) +  
    scale_shape_manual(values=shape_palette) +  
    scale_fill_manual(values = color_palette, guide = "none") + 
    labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
         y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
         color="Age Group", shape="Age Group") +  
    theme_classic(base_size = 15) + 
    theme(legend.position="right", 
          legend.title=element_text(size=10, family="Times New Roman"), 
          legend.text=element_text(size=10, family="Times New Roman"), 
          legend.key.size = unit(1.5, "lines"),  
          axis.title=element_text(size=10, family="Times New Roman"), 
          axis.text=element_text(size=10, family="Times New Roman"),
          plot.title = element_text(size = 16, family="Times New Roman", hjust = 0.5)) +  
    ggtitle("PCoA of Age")  
p_age

color_palette <- c("female" = "#E69A7A", "male" = '#7EB48F')  
shape_palette <- c("female" = 17, "male" = 16) 
p_Sex <-ggplot(dune_pcoa_result, aes(x=PCoA1, y=PCoA2, color=factor(Sex), shape=factor(Sex))) + 
    geom_point(size=3.5, alpha=1) +  
    stat_ellipse(level=0.95, linetype=1, linewidth=0.8, aes(fill = factor(Sex)),geom = "polygon", alpha = 0.3) +  
    scale_color_manual(values=color_palette) +  
    scale_shape_manual(values=shape_palette) +  
    scale_fill_manual(values = color_palette, guide = "none") +  
    labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
         y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
         color="Sex Group", shape="Sex Group") +  
    theme_classic(base_size = 15) +  
    theme(legend.position="right", 
          legend.title=element_text(size=10, family="Times New Roman"), 
          legend.text=element_text(size=10, family="Times New Roman"), 
          legend.key.size = unit(1.5, "lines"),  
          axis.title=element_text(size=10, family="Times New Roman"), 
          axis.text=element_text(size=10, family="Times New Roman"),
          plot.title = element_text(size = 16, family="Times New Roman", hjust = 0.5)) +  
    ggtitle("PCoA of Sex")  
p_Sex

#---------------------------------------------------------------------------------------------------------
library(vegan)
dune.env_df <- as.data.frame(dune.env)
set.seed(8) 
dune.div <- adonis2(dune ~ Age_group+Sex+Age_group+BMI_group+Smoking+Alcohol, data = dune.env_df, permutations = 999, method="bray", by="margin")
dune.div

#--------------------------------------------------------------------------------
library(ggsci)
library(extrafont)
library(vegan)
library(ggplot2)
r2_value <- round(dune.div$R2[1], 4)
p_value <- round(dune.div$`Pr(>F)`[1],4)
label_text <- sprintf("adonis R²: %.4f\nP-value:    %.4f", r2_value, p_value)
p_age_2 <- p_age +
    annotate("text", x = Inf, y = Inf, label = label_text, 
             hjust = 1.1, vjust = 2, size = 3.5, family = "Times New Roman")
p_age_2

r2_value <- round(dune.div$R2[2], 4)
p_value <- round(dune.div$`Pr(>F)`[2],4)
label_text <- sprintf("adonis R²: %.4f\nP-value:    %.4f", r2_value, p_value)
p_Sex_2 <- p_Sex +
    annotate("text", x = Inf, y = Inf, label = label_text, 
             hjust = 1.0, vjust = 2, size = 3.5, family = "Times New Roman")
p_Sex_2

#---------PERMANOVA------------
adonis_result <- data.frame(
    Factor = c("Age", "Sex", "BMI", "Smoking", "Alcohol"),
    R2 = c(dune.div$R2[1:5]),
    P = c(dune.div$`Pr(>F)`[1:5])
)

adonis_result$Significance <- cut(
    adonis_result$P,
    breaks = c(-Inf, 0.001, 0.01, 0.05, 1),
    labels = c("***", "**", "*", "")
)

adonis_result$Factor <- factor(adonis_result$Factor, levels = adonis_result$Factor[order(adonis_result$R2, decreasing = F)])

adonis_result$logP <- -log10(adonis_result$P)
library(ggplot2)

p_PERMANOVA <- ggplot(adonis_result, aes(x = R2, y = Factor)) +
    geom_bar(stat = "identity", aes(fill = logP), width = 0.5) +
    geom_text(aes(label = Significance), hjust = -0.2, color = "black", size = 5, family = "Times New Roman") +
    scale_fill_gradient(name = "-log10(P)", low = "lightblue", high = "darkblue") +
    labs(x = expression(R^2), y = "Factor", title = "PERMANOVA with Bray–curtis distance") +
    coord_cartesian(xlim = c(0, max(adonis_result$R2) )) +  
    theme_minimal(base_size = 15) +  
    theme(
        axis.text.x = element_text(size = 12, family = "Times New Roman"),
        axis.text.y = element_text(size = 12, family = "Times New Roman"),
        axis.title = element_text(size = 14, family = "Times New Roman"),
        plot.title = element_text(size = 16, family = "Times New Roman", hjust = 0.5),
        panel.grid.minor = element_blank(), 
    )
print(p_PERMANOVA)