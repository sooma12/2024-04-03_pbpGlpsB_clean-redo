# deseq_2024-04-04.R
# DESeq2

# Installed these packages using Packrat:
# install.packages("tidyverse")
# install.packages("vsn") ## not for this R version
# install.packages("ggrepel")
# install.packages("DEGreport") ## not for this R version
# install.packages("pheatmap")
# BiocManager::install(pattern="apeglm")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
#.libpaths()
#BiocManager::install("apeglm")


# Set up packrat (package manager)
library(packrat)
pr_dir <- '/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/software'
setwd(pr_dir)
packrat::init()


# Load libraries
library(DESeq2)
library(ggplot2)
#library(vsn)
library(dplyr)
library(tidyverse)
library(ggrepel)
# library(DEGreport)
library(pheatmap)
library(apeglm)

wd <- "/work/geisingerlab/Mark/rnaSeq/2024-04-03_pbpGlpsB_clean-redo"
setwd(wd)

## Read in raw data and prepare feature count table
feature_count <- read.table("./data/featurecounts/counts.txt", header=TRUE, row.names = 1)
# Grab the count data
data <- feature_count[,6:14]
# view(data)

## Set up column names in feature count data AND prep metadata table with strains and conditions
# Trim column names down to just the sample IDs
column_names <- colnames(data)
column_names <- sub("X.work.geisingerlab.Mark.rnaSeq.2024.04.03_pbpGlpsB_clean.redo.data.mapped.", "", column_names)
column_names <- sub("Aligned.sortedByCoord.out.bam", "", column_names)
column_names
colnames(data) <- column_names
# Use regex to get condition (mutation) from strain IDs
# ".+?(?=_)":   .+?  means "match any character (.) any number of times (+?)"
# (?=_): a positive lookahead: find where the character - is matched
conditions <- str_extract(column_names, ".+?(?=_)")
# Use column_names and conditions to make metadata table
meta <- data.frame(column_names, conditions)
colnames(meta) <- c('id', 'condition')
meta  # Verify that IDs and conditions are as expected
# Write out metadata table (save!!)
#meta_filepath <- '/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/data'
#meta_file <- file.path(meta_filepath, 'metadata.txt')
#write.table(meta, meta_file, row.names=FALSE)  # TODO: Come back to this. Ok to have "id" label, or needs to be blank as in examples?  May want to change colnames before this.

## Load data, pre-filter low-count genes, and relevel to set WT as the reference
des_data <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ condition)
smallestGroupSize <- 3  # should be size of smallest group; I did 3 replicates
keep <- rowSums(counts(des_data) >= 10) >= smallestGroupSize  # keep data where count is >10 in all 3 samples
# relevel dds$condition to set WT as the reference
des_data$condition <- relevel(des_data$condition, ref = "WT")

dds <- DESeq(des_data)

# LFC shrinkage helps visualize and rank genes...
resultsNames(dds)  # get names of coefficients to shrink
## "Intercept" "condition_DlpsB_vs_WT" "condition_DpbpG_vs_WT"
# Use default shrinkage, which is apeglm method (Zhu, Ibrahim, and Love 2018)
resLFC_pbpG <- lfcShrink(dds, coef="condition_DpbpG_vs_WT", type="apeglm")
resLFC_lpsB <- lfcShrink(dds, coef="condition_DlpsB_vs_WT", type="apeglm")

# Clean up row names
pbpG_rownames <- sub('gene-', '', row.names(resLFC_pbpG))
lpsB_rownames <- sub('gene-', '', row.names(resLFC_lpsB))
row.names(resLFC_pbpG) <- pbpG_rownames
row.names(resLFC_lpsB) <- lpsB_rownames

# Count significant hits
sum(resLFC_pbpG$padj < 0.1, na.rm=TRUE)
sum(resLFC_pbpG$padj < 0.1 & abs(resLFC_pbpG$log2FoldChange) > 1, na.rm=TRUE)
sum(resLFC_lpsB$padj < 0.1, na.rm=TRUE)
sum(resLFC_lpsB$padj < 0.1 & abs(resLFC_lpsB$log2FoldChange) > 1, na.rm=TRUE)

# Check largest foldchanges from significant hits
resSigPbpG <- resLFC_pbpG[which(resLFC_pbpG$padj < 0.1), ]

#downregulated
head(resSigPbpG[order(resSigPbpG$log2FoldChange),])  # Verify pbpG (RS16895) is top hit
#upregulated
tail(resSigPbpG[order(resSigPbpG$log2FoldChange),])

resSigLpsB <- resLFC_lpsB[which(resLFC_lpsB$padj < 0.1), ]
#downregulated
head(resSigLpsB[order(resSigLpsB$log2FoldChange),])  # Verify lpsB (RS15950) is top hit
#upregulated
tail(resSigLpsB[order(resSigLpsB$log2FoldChange),])  # Verify lpsB (RS15950) is top hit

# MA plot
plotMA(resLFC_pbpG, ylim=c(-1, 1))
plotMA(resLFC_lpsB, ylim=c(-1, 1))

# Check PCA plot
rld <- rlog(dds, blind=TRUE)
plotPCA(rld, intgroup="condition") + geom_text(aes(label=name))

# Plot the gene with the smallest value 
#plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
plotCounts(dds, gene="gene-ACX60_RS15100", intgroup="condition")

# Save data
#file_name = 'NAME/OF/FILE'
#write.table(<results/go/here>, file_name, quote=FALSE)
# Can also pre-filter for significant hits
#padj.cutoff <- 0.05 # False Discovery Rate cutoff
#significant_results <- res[which(res$padj < padj.cutoff),]
#write.table(significant_results, file_name, quote=FALSE)

#get session information, including package versions... can save to a variable, then write out.
#print(sessionInfo())


# Save current analyses
out_dir <- "/work/geisingerlab/Mark/rnaSeq/2024-04-03_pbpGlpsB_clean-redo"
setwd(out_dir)

# Save current analyses
pbpGresOrdered <- resLFC_pbpG[order(resLFC_pbpG$pvalue),]
write.csv(as.data.frame(pbpGresOrdered), file="DES_pbpG_2024-04-04.csv")

lpsBresOrdered <- resLFC_lpsB[order(resLFC_lpsB$pvalue),]
write.csv(as.data.frame(lpsBresOrdered), file="DES_lpsB_2024-04-04.csv")

### Volcano plots
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

pbpG_data_to_plot <- as.data.frame(pbpGresOrdered)

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
pbpG_data_to_plot$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
pbpG_data_to_plot$diffexpressed[pbpG_data_to_plot$log2FoldChange > 0.6 & pbpG_data_to_plot$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
pbpG_data_to_plot$diffexpressed[pbpG_data_to_plot$log2FoldChange < -0.6 & pbpG_data_to_plot$padj < 0.05] <- "DOWN"
# Explore a bit
head(pbpG_data_to_plot[order(pbpG_data_to_plot$padj) & pbpG_data_to_plot$diffexpressed == 'DOWN', ])

# Labels for points
pbpG_data_to_plot$gene_symbol <- row.names(pbpG_data_to_plot)
pbpG_data_to_plot$delabel <- ifelse(pbpG_data_to_plot$gene_symbol %in% head(pbpG_data_to_plot[order(pbpG_data_to_plot$padj), "gene_symbol"], 30), pbpG_data_to_plot$gene_symbol, NA)

ggplot(data = pbpG_data_to_plot, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 1) + 
  scale_color_manual(values = c("#F28169", "grey", "#00AFBB"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +# to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 90), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Condition', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"(adjusted p-value)")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  geom_label_repel(aes(label = delabel),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   max.overlaps = Inf)
