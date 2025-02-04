library(edgeR)
library(Homo.sapiens)
library(gage)
library(pathview)
library(EnhancedVolcano)

# GO and KEGG analysis with edgeR ----------
load("out/DGEList_UNvsWT.Rdata")
load("out/exactTest_UNvsWT.Rdata")

# convert to entrez ids 
entrez_ids <- select(org.Hs.eg.db, keys=rownames(et),
                 columns=c("ENSEMBL","ENTREZID"), keytype="ENSEMBL")

et$table <- merge(et$table, entrez_ids, by.x = 0, by.y = "ENSEMBL")
et$genes <- merge(et$genes, entrez_ids, by.x = 1, by.y = "ENSEMBL")

# GO enrichment
go <- goana(et, geneid = 'ENTREZID')
write.csv(go, file = "out/GO_UNvsWT.csv")
topGO(go, ontology = "BP", n = 20)

# KEGG enrichment
kk <- kegga(et, geneid = "ENTREZID")
topKEGG(kk, sort = "Up")
write.csv(kk, file = "out/KEGG_UNvsWT.csv")

# GAGE and Pathview ------------
data("kegg.gs") # load the kegg databases

edger.fc=et$table$logFC
names(edger.fc)=et$table$ENTREZID

# Analysis for up-regulated pathways
fc.kegg.p <- gage(edger.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
sel <- fc.kegg.p$greater[, "q.val"] < 0.1 &
  !is.na(fc.kegg.p$greater[, "q.val"]) # select for significant pathways
path.ids <- rownames(fc.kegg.p$greater)[sel]

sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 &
   !is.na(fc.kegg.p$less[,"q.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)

# visualize the top 3 pathways
pv.out.list <- sapply(path.ids2[1:3], function(pid) pathview(
  gene.data = edger.fc, pathway.id = pid, 
  low = c("gene" = "red"), high = "green"))


plin_pathway <- pathview(gene.data = edger.fc, pathway.id = "hsa04923", 
                         low = c("gene" = "red"), high = "green")

