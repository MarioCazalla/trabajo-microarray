library(org.Hs.eg.db)
library(AnnotationDbi)
library(data.table)
library(hgu133plus2.db)

gsi <- fread("/home/guille/Desktop/Omics/microarrays/trabajo-microarray/GSEA_files/gsi-notch.txt", header = F)
gsi <- gsi[!duplicated(gsi$V1),]

AnnotationDbi::select(org.Hs.eg.db, gsi$V1, keytype = "SYMBOL", columns = c("ENTREZID", "ENSEMBL"))


#write.table(gsi$V1, file = "/home/guille/Desktop/Omics/microarrays/trabajo-microarray/GSEA_files/gsi-notch.txt", quote = F, row.names = F, col.names = F)

gsi <- t(gsi)
write.table(t(gsi), file = "/home/guille/Desktop/Omics/microarrays/trabajo-microarray/GSEA_files/gsi-notch.gmt", sep = "\t", quote = F, row.names = F, col.names = F)

lqf <- gsi <- fread("/home/guille/Desktop/Omics/microarrays/trabajo-microarray/GSEA_files/los_que_faltaban", header = F)
AnnotationDbi::select(org.Hs.eg.db, c("RASGRP1", "PTPRC", "EVI2A", "METTL7A", "RORC", "CD69", "ELOVL4"), keytype = "SYMBOL", columns = c("ENTREZID", "ENSEMBL"))

AnnotationDbi::select(hgu133plus2.db, "207226_at", keytype = "PROBEID", columns = c("SYMBOL", "ENSEMBL"))


gsh <- fread("/home/guille/Desktop/Omics/microarrays/trabajo-microarray/GSEA_files/gene_set_hortensia.txt", header = F)

AnnotationDbi::select(org.Hs.eg.db, gsh$V1, keytype = "SYMBOL", columns = c("ENTREZID", "ENSEMBL"))
