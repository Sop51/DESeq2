library(dplyr)
library(purrr)
library(pasilla)
library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(apeglm)
library(ggplot2)
library(org.Dr.eg.db)

# ------------------ formatting the data for analysis ----------------------#
# Read the file
countData <- read.delim('/Users/sm2949/Desktop/fragmentCounts.txt', header=TRUE, sep = "\t", skip=1)
# subset to only include gene id and counts
countData <- countData[, c(1, 7:12)]
# rename the columns for clarity
colnames(countData) <- c("Geneid", "WF1", "WF2", "WF3", "WM1", "WM2", "WM3")
# set Geneid as the row names
rownames(countData) <- countData$Geneid
# drop the Geneid column
countData <- countData[, -which(names(countData) == "Geneid")]

# create a metadata df
metadata <- data.frame(
    sample = c('WF1', 'WF2', 'WF3', 'WM1', 'WM2','WM3'),
    sex = c('F','F','F','M','M','M')
)
  
# set the row names as sample names 
rownames(metadata) <- metadata$sample
metadata$sample <- NULL

# convert sex to a factor
metadata$sex <- as.factor(metadata$sex)

# ----------------------- DESEQ2 -----------------------------#
# create the deseq object on ~ sex controlling for fish samples
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = metadata,
                              design = ~ sample + sex)

# set M as the reference level
dds$sex <- relevel(dds$sex, ref = "M")

# pre filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# get normalized counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# log normalize the counts to use for heatmap plotting
rld <- rlog(dds, blind = FALSE)

# perform differential expression 
dds <- DESeq(dds)
# perform shrinkage
resLFC <- lfcShrink(dds, coef="sex_F_vs_M", type="apeglm")
# convert to a dataframe
results_df <- as.data.frame(resLFC)
# save the results to a csv
write.csv(results_df,file='/Users/sm2949/Desktop/DESEQ_MvF_results.csv',na='')

# remove rows with NA adjusted pvals
results_df <- results_df[!is.na(results_df$padj), ]

# map ensembl ids to gene symbols
gene_symbols <- mapIds(org.Dr.eg.db, keys = rownames(results_df), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# add the gene symbols to the results dataframe
results_df$gene_symbol <- gene_symbols[match(rownames(results_df), names(gene_symbols))]

#--------------------- performing KEGG and GO --------------------#
# filter results for significant genes with padj < 0.05
top_de_genes <- resLFC[which(resLFC$padj < 0.05), ]

# split into upregulated and downregulated genes
upregulated_genes <- top_de_genes_sorted[top_de_genes_sorted$log2FoldChange > 0, ]
downregulated_genes <- top_de_genes_sorted[top_de_genes_sorted$log2FoldChange < 0, ]

# get the top 500 upregulated and downregulated genes
top_upregulated_genes <- head(upregulated_genes, 1000)
top_downregulated_genes <- head(downregulated_genes,1000)

# get entrez IDs for top 500 upregulated and downregulated genes
upregulated_entrez_ids <- mapIds(org.Dr.eg.db, 
                                 keys = rownames(top_upregulated_genes), 
                                 column = "ENTREZID", 
                                 keytype = "ENSEMBL", 
                                 multiVals = "first")

downregulated_entrez_ids <- mapIds(org.Dr.eg.db, 
                                   keys = rownames(top_downregulated_genes), 
                                   column = "ENTREZID", 
                                   keytype = "ENSEMBL", 
                                   multiVals = "first")

# perform KEGG enrichment for top 500 upregulated and downregulated genes
kegg_up <- enrichKEGG(gene = upregulated_entrez_ids,
                      organism = "dre",  # Zebrafish KEGG code
                      pvalueCutoff = 0.05)

kegg_down <- enrichKEGG(gene = downregulated_entrez_ids,
                        organism = "dre",  # Zebrafish KEGG code
                        pvalueCutoff = 0.05)

# perform GO enrichment for top 500 upregulated and downregulated genes
go_up <- enrichGO(gene = upregulated_entrez_ids,
                  OrgDb = org.Dr.eg.db, 
                  keyType = "ENTREZID", 
                  ont = "BP", 
                  pvalueCutoff = 0.05)

go_down <- enrichGO(gene = downregulated_entrez_ids,
                    OrgDb = org.Dr.eg.db, 
                    keyType = "ENTREZID", 
                    ont = "BP", 
                    pvalueCutoff = 0.05)

# create individual plots
kegg_up_plot <- dotplot(kegg_up, showCategory = 15) + 
  ggtitle("Upregulated KEGG Terms") + 
  theme_minimal()

kegg_down_plot <- dotplot(kegg_down, showCategory = 15) + 
  ggtitle("Downregulated KEGG Terms") + 
  theme_minimal()

go_up_plot <- dotplot(go_up, showCategory = 15) + 
  ggtitle("Upregulated GO Terms") + 
  theme_minimal()

go_down_plot <- dotplot(go_down, showCategory = 15) + 
  ggtitle("Downregulated GO Terms") + 
  theme_minimal()

# arrange the plots side by side
grid.arrange(kegg_up_plot, kegg_down_plot, go_up_plot, go_down_plot, ncol = 2)
