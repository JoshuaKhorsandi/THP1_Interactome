library(edgeR)
library(Glimma)
library(tidyverse)
library(Homo.sapiens)
library(RColorBrewer)
library(pheatmap)
library(biomaRt)

# Data Organization ----------------------------
counts_matrix <- read.csv("../data/raw_counts_GRCh38_GREIS.csv") %>% 
  column_to_rownames(var = "X")
  # dplyr::select(c(1:5, 7:9)) |> # compare uninfected to WT infected cells, triplicate

group <- as.factor(c(rep("UN", 3), rep("WT", 3)))
x <- DGEList(counts_matrix, group = group)
x$samples

# first column contains old gene symbols, updated version for HGNC ids
# geneid <- rownames(x) # Entrez gene ids from Salmon
# genes <- select(Homo.sapiens, keys = geneid, columns = c("SYMBOL", "TXCHROM"), 
                # keytype = "ENSEMBL")
# genes <- genes[!duplicated(genes$ENSEMBL),]
head(genes)
x

# Normalization and Filtering --------------------------------
cpm <- cpm(x)
lcpm <- cpm(x, log=T) # log2(CPM + 2/L) to add "prior" value
L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6

table(rowSums(x$counts==0)==6) # zero counts across all samples

# Figure 1 - The Density of log-CPM values before/after filtering for low genes
samplenames <- colnames(x$counts)
lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")

# A - Raw Data
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

# Filter
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

# B - Filtered Data
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

# calculate normalization with TMM (Oshlack 2010)
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

# Figure 2 MDS Plot
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
plotMDS(lcpm, labels=group)
title(main="Sample groups MDS")

# Negative Binomial Models for Differential Expression -------------
x <- estimateDisp(x)
et <- exactTest(x)
summary(decideTests(et))

results <- topTags(et, n = Inf)
write.csv(results$table, file = "out/UN_vs_WT_DGE.csv", row.names = F, quote = F)
save(x, file = "out/DGEList_UNvsWT.Rdata")
save(et, file = "out/exactTest_UNvsWT.Rdata")

# Interactive glimma figures
glimmaMA(et)
glimmaVolcano(et, dge = x_nb)

# Figure 3, heatmaps
top_genes <- topTags(et, n = 100)
lcpm_norm <- as.data.frame(cpm(x,log=TRUE))

i <- which(et$genes$gene_symbol %in% top_genes$table$gene_symbol)

heat_plot <- pheatmap(lcpm_norm[i,], 
                      col = brewer.pal(10, 'RdYlGn'),
                      cluster_rows = T, cluster_cols = T, trace = "none",
                      angle_col = 45, 
                      fontsize_row = 6, fontsize_col = 8, show_rownames = T, 
                      main = " Uninfected vs. WT THP1", labels_row = et$genes$gene_symbol, 
                      labels_col = group)


# Exports to Pathway Analysis --------------------
# See - Expression Map Protocol

# g:Profiler - gene list
tt_exact_test <- topTags(et,n=nrow(x))
sig_genes = which(tt_exact_test$table$FDR < 0.005)
length(sig_genes)

topgenes_qval <- rownames(tt_exact_test$table)[sig_genes]
head(topgenes_qval)
write.table(topgenes_qval, "out/qval0005_UNvsWT.txt", sep = "\t", row.names = F, col.names = F, quote = F)

# GSEA ranks file 
ranks_RNAseq = sign(tt_exact_test$table$logFC) * -log10(tt_exact_test$table$PValue)

# create the ranks file
ranks_RNAseq <- cbind(rownames(tt_exact_test$table), ranks_RNAseq)
colnames(ranks_RNAseq) <- c("GeneID","rank")
ranks_RNAseq <- ranks_RNAseq[order(as.numeric(ranks_RNAseq[,2]),decreasing = TRUE),]
head(ranks_RNAseq) # the top2 values have pval ~ 0 so produce Inf rank

write.table(ranks_RNAseq, file = "out/UNvsWT.rnk", sep = "\t", quote = F, row.names = F)

# Expression Map file 
options(RCurlOptions=list(followlocation=TRUE, postredir=2L)) # bio mart redirect 
normalized_expression_RNAseq <- cpm(x)

EM_expressionFile_RNAseq <- data.frame(Description = x$genes$gene_symbol, normalized_expression_RNAseq) # set the description to be the symbol
head(EM_expressionFile_RNAseq)

write.table(EM_expressionFile_RNAseq ,file = "out/EM_expressionFile_RNAseq.txt", sep = "\t", row.names = T, quote = F)



