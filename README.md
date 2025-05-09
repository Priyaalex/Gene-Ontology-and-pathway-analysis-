# Gene-Ontology-and-pathway-analysis-
This repository contains an R script designed to perform Gene Ontology enrichment and pathway analysis for a given set of genes. 
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
###Cluster Profiler
library(org.Hs.eg.db)
library(clusterProfiler)
x <- read.table("GENELIST.txt", header=TRUE, sep="\t", as.is=TRUE)
m <- x$GENESYMBOL
#n <- select(org.Hs.eg.db, m, c("ENTREZID", "GENENAME"), "ALIAS")
n <- AnnotationDbi::select (org.Hs.eg.db, m, c("ENTREZID", "GENENAME"), "ALIAS")
head(n)
write.csv(n, file="gene_entrezid.csv")
o <- n$ENTREZID
aa <- enrichGO(o, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.05, readable = TRUE)
dotplot(aa) 
write.csv(aa, file="BP.csv")
aa <- enrichGO(o, 'org.Hs.eg.db', ont="MF", pvalueCutoff=0.05, readable = TRUE)
dotplot(aa) 
write.csv(aa, file="MF.csv")
aa <- enrichGO(o, 'org.Hs.eg.db', ont="CC", pvalueCutoff=0.05, readable = TRUE)
dotplot(aa)  
write.csv(aa, file="CC.csv")
ap <- enrichKEGG(o, organism = "hsa", pvalueCutoff = 0.05)
dotplot(ap) 
write.csv(ap, file="KEGG.csv")
savehistory(file = "command.Rhistory")
