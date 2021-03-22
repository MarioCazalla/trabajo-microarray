### Establecer directorio de trabajo y cargar librerias
setwd("/home/guille/Desktop/Omics/microarrays/trabajo-microarray/")

library(affy)
library(limma)
library(genefilter)
library(ArrayTools)

library(Biobase)
library(GEOquery)
library(annotate)
library(hgu133plus2.db)

library(topGO)
library(clusterProfiler)
library(ggplot2)
library(stringr)

#2. Import targets.txt file. generates a data frame:
targets <- readTargets("targets.txt", row.names="FileName")
targets.kopt <- readTargets("targets_kopt.txt", row.names="FileName")
targets.hpb <- readTargets("targets_hpb.txt", row.names="FileName")

#3. Import .CEL files
data <- ReadAffy(filenames=targets$FileName)
data.kopt <- ReadAffy(filenames=targets.kopt$FileName)
data.hpb <- ReadAffy(filenames=targets.hpb$FileName)


#4. Normalize with RMA
normalize.data <- function(affybatch){
  expresso(affybatch,
           bg.correct = TRUE, 
           bgcorrect.method="rma",
           normalize = TRUE, 
           normalize.method="quantiles", 
           pmcorrect.method="pmonly", 
           summary.method="medianpolish",
           verbose = TRUE)
}

eset <- normalize.data(data)
eset.kopt <- normalize.data(data.kopt)
eset.hpb <- normalize.data(data.hpb)

ematrix <- exprs(eset)
ematrix.kopt <- exprs(eset.kopt)
ematrix.hpb <- exprs(eset.hpb)

esetpd <- pData(eset)
esetpd.kopt <- pData(eset.kopt)
esetpd.hpb <- pData(eset.hpb)

labels<-c(rep("KOPT_K1_DMSO",3), rep("HPB_ALL_DMSO",3),
          rep("KOPT_K1_SAHM1",3), rep("HPB_ALL_SAHM1",3))
labels.kopt <- c(rep("KOPT_K1_DMSO",3), rep("KOPT_K1_SAHM1",3))
labels.hpb <- c(rep("HPB_ALL_DMSO",3), rep("HPB_ALL_SAHM1",3))

colnames(ematrix) <- labels
colnames(ematrix.kopt) <- labels.kopt
colnames(ematrix.hpb) <- labels.hpb
rm(labels, labels.hpb, labels.kopt)

#Filtering:
iqr.filtering <- function(expressionset){
  filtered <- varFilter(expressionset, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)
  return(filtered)
}
esetIQR <- iqr.filtering(eset)
esetIQR.kopt <- iqr.filtering(eset.kopt)
esetIQR.hpb <- iqr.filtering(eset.hpb)

###### PCA analysis ######

library(PCAtools)

pca <- pca(ematrix, metadata = esetpd)
pca.kopt <- pca(ematrix.kopt, metadata = esetpd.kopt)
pca.hpb <- pca(ematrix.hpb, metadata = esetpd.hpb)
a <- list(pca, pca.kopt, pca.hpb)

for(i in a){
  jpeg(paste("figuras/",deparse(substitute(i)),"_screeplot.jpg", sep = ""), width = 600, height = 600)
  screeplot(i)
  dev.off()
  jpeg(paste("figuras/",deparse(substitute(i)),"_plot.jpg", sep = ""), width = 600, height = 600)
  biplot(pca1, lab = labels)
  dev.off()
}



#5. Generate BOXPLOTS before and after normalization

#boxplot for raw data
jpeg("figuras/boxplot_not_normalised.jpg", width = 600, height = 600)
boxplot(data,
        main="Boxplot Before Normalization",
        col = "lightgrey",
        las=2)
dev.off()

#boxplot for normalized data
exprseset <- as.data.frame(exprs(eset))		

jpeg("figuras/boxplot_normalised.jpg", width = 600, height = 600)
boxplot(data.frame(exprseset),
        main="Boxplot After Normalization (log scale)",
        col = "white",
        las=2)
dev.off()
rm(exprseset)


#######Differential expression analysis.#######

#7. Design matrix.
# Para el análisis conjunto por bloques:
cell_type_blocking <- factor(targets$CellType)
treat_blocking <- factor(targets$Treatment)
design_blocking <- model.matrix(~cell_type_blocking+treat_blocking)
rownames(design_blocking)<-targets$FileName

# Para el análisis de líneas celulares por separado y comparando ambas:
design<-cbind(KOPT_K1_DMSO=c(1,1,1,0,0,0,0,0,0,0,0,0),
              HPB_ALL_DMSO=c(0,0,0,1,1,1,0,0,0,0,0,0),
              KOPT_K1_SAHM1=c(0,0,0,0,0,0,1,1,1,0,0,0),
              HPB_ALL_SAHM1=c(0,0,0,0,0,0,0,0,0,1,1,1)
)

rownames(design)<-targets$FileName

#8. Contrasts matrix.
#Cambios en la expresion diferentes en las lineas celulares
cont.matrix_lines<-makeContrasts((KOPT_K1_SAHM1 - KOPT_K1_DMSO) - (HPB_ALL_SAHM1 - HPB_ALL_DMSO), levels=design)
cont.matrix_kopt<-makeContrasts(KOPT_K1_SAHM1 - KOPT_K1_DMSO, levels=design)
cont.matrix_hpb<-makeContrasts(HPB_ALL_SAHM1 - HPB_ALL_DMSO, levels=design)


#9. Obtaining differentially expressed genes (DEGs)
# Primero para el analisis por bloques:
fit<-lmFit(esetIQR, design_blocking)
fit2<-eBayes(fit)
toptableIQR_blocking<-topTable(fit2, coef="treat_blockingSAHM1", number=dim(exprs(esetIQR))[1], adjust.method="BH", sort.by="p")
toptableIQR_blocking_nsig <- nrow(toptableIQR_blocking[toptableIQR_blocking$adj.P.Val<0.05, ])
rm(fit); rm(fit2)

deg <- function(expression_set, design_matrix, contrasts_matrix){
  fit<-lmFit(expression_set,design_matrix) 
  fit2<-contrasts.fit(fit, contrasts_matrix)
  fit2<-eBayes(fit2)
  toptable <- topTable(fit2, number=dim(exprs(expression_set))[1], adjust.method="BH", sort.by="p")
  significant_results <- nrow(toptable[toptable$adj.P.Val<0.05, ])
  
  results<-list(toptable, significant_results)
  return(results)
}

res_lines <- deg(esetIQR, design, cont.matrix_lines)
toptableIQR_lines <- res_lines[[1]]
toptableIQR_lines_nsig <- res_lines[[2]]
rm(res_lines)

res_kopt <- deg(esetIQR, design, cont.matrix_kopt)
toptableIQR_kopt <- res_kopt[[1]]
toptableIQR_kopt_nsig <- res_kopt[[2]]
rm(res_kopt)

res_hpb <- deg(esetIQR, design, cont.matrix_hpb)
toptableIQR_hpb <- res_hpb[[1]]
toptableIQR_hpb_nsig <- res_hpb[[2]]
rm(res_hpb)

#Ver cuantos cambios significativos son mayores/menores en KOPT-K1 que en HPB-ALL
toptableIQR_lines$sign <- sign(toptableIQR_lines$logFC)
toptableIQR_lines_filtered <- toptableIQR_lines[toptableIQR_lines$adj.P.Val<=0.05,]
sign_count_lines=table(toptableIQR_lines_filtered$sign)
sign_count_lines
rm(toptableIQR_lines_filtered)
toptableIQR_lines$sign <-NULL

# Comparar los analisis independientes kopt y hpb:
toptableIQR_kopt_sig <- toptableIQR_kopt[toptableIQR_kopt$adj.P.Val<=0.05,]
toptableIQR_hpb_sig <- toptableIQR_hpb[toptableIQR_hpb$adj.P.Val<=0.05,]
a <- rownames(toptableIQR_kopt_sig)
b <- rownames(toptableIQR_hpb_sig)

intersect_kopt_hpb <- intersect(a, b)
length(intersect_kopt_hpb) # prints '1110'
rm(a); rm(b)
toptableIQR_hpb_nsig - length(intersect_kopt_hpb) # prints '1717'
toptableIQR_kopt_nsig - length(intersect_kopt_hpb) # prints '6289'

##10. Save results
#save(toptableIQR,file="CellLinesResults.RData")
#load("/home/guille/Desktop/Omics/microarrays/Trabajo_final/CellLinesResults.RData")

######Anotación######

anno <- function(toptable){
  affy_IDs <- as.character(row.names(toptable))
  
  GENE_SYMBOLS <- AnnotationDbi::select(hgu133plus2.db, keys=affy_IDs, columns="SYMBOL", keytype="PROBEID")
  ENSEMBL_IDs <-AnnotationDbi::select(hgu133plus2.db, keys=affy_IDs, columns="ENSEMBL", keytype="PROBEID")
  
  GENE_SYMBOLS = GENE_SYMBOLS[!duplicated(GENE_SYMBOLS["PROBEID"]),]
  ENSEMBL_IDs = ENSEMBL_IDs[!duplicated(ENSEMBL_IDs["PROBEID"]),]
  
  toptable$SYMBOL <- GENE_SYMBOLS$SYMBOL
  toptable$ENSEMBL <- ENSEMBL_IDs$ENSEMBL
  toptable$rank <- c(1:nrow(toptable))
  toptable<-toptable[, c(9,8,7,1,2,3,4,5,6)]
  
  return(toptable)
}

toptableIQR_blocking <- anno(toptableIQR_blocking)
toptableIQR_lines <- anno(toptableIQR_lines)
toptableIQR_kopt <- anno(toptableIQR_kopt)
toptableIQR_hpb <- anno(toptableIQR_hpb)


######ORA: Overrepresentation Analysis: GO terms######

GO_analysis <- function (genes, universe, ontology){
  clusterProfiler::enrichGO(gene          = genes,
           universe      = universe,
           OrgDb         = hgu133plus2.db,
           keyType       = 'PROBEID',
           ont           = ontology, # "BP", "MF", "CC" o "ALL"
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.05)
}

universe <- rownames(ematrix)

FC_sign_subsetting <- function(toptable){
  toptable$logFCsign <- sign(toptable$logFC)
  pos <- toptable[toptable$logFCsign==1,]
  neg <- toptable[toptable$logFCsign==-1,]
  res <- list(pos, neg)
  return(res)
}

# GO tems usando el diseño por bloques:

#~~~~~~Blocking~~~~~~
toptableIQR_blocking_sig <- toptableIQR_blocking[toptableIQR_blocking$adj.P.Val<=0.05,]
toptableIQR_blocking_sig_signs <- FC_sign_subsetting(toptableIQR_blocking_sig)
toptableIQR_blocking_sig_pos <- toptableIQR_blocking_sig_signs[[1]]
toptableIQR_blocking_sig_neg <- toptableIQR_blocking_sig_signs[[2]]
rm(toptableIQR_blocking_sig_signs)

GO_blocking_pos <- GO_analysis(rownames(toptableIQR_blocking_sig_pos), universe, "BP")
GO_blocking_neg <- GO_analysis(rownames(toptableIQR_blocking_sig_neg), universe, "BP")

jpeg("figuras/GOBP_blocking_pos.jpg", width = 900, height = 600)
dotplot(GO_blocking_pos, showCategory = 10, font.size = 10, title = "GO:BP Upregulation. Blocking KOPT-K1, HPB-ALL")
dev.off()

jpeg("figuras/GOBP_blocking_neg.jpg", width = 900, height = 600)
dotplot(GO_blocking_neg, showCategory = 10, font.size = 10, title = "GO:BP Downregulation. Blocking KOPT-K1, HPB-ALL")
dev.off()


# GO terms usando los genes significativos del analisis kopt y hpb:

#~~~~~~KOPT-K1~~~~~~
toptableIQR_kopt_sig_signs <- FC_sign_subsetting(toptableIQR_kopt_sig)
toptableIQR_kopt_sig_pos <- toptableIQR_kopt_sig_signs[[1]]
toptableIQR_kopt_sig_neg <- toptableIQR_kopt_sig_signs[[2]]
rm(toptableIQR_kopt_sig_signs)

GO_kopt_pos <- GO_analysis(rownames(toptableIQR_kopt_sig_pos), universe, "BP")
GO_kopt_neg <- GO_analysis(rownames(toptableIQR_kopt_sig_neg), universe, "BP")

jpeg("figuras/GOBP_kopt_pos.jpg", width = 900, height = 600)
dotplot(GO_kopt_pos, showCategory = 10, font.size = 10, title = "GO:BP Upregulation KOPT-K1")
dev.off()

jpeg("figuras/GOBP_kopt_neg.jpg", width = 900, height = 600)
dotplot(GO_kopt_neg, showCategory = 10, font.size = 10, title = "GO:BP Downregulation KOPT-K1")
dev.off()


#~~~~~~HPB-ALL~~~~~~
toptableIQR_hpb_sig_signs <- FC_sign_subsetting(toptableIQR_hpb_sig)
toptableIQR_hpb_sig_pos <- toptableIQR_hpb_sig_signs[[1]]
toptableIQR_hpb_sig_neg <- toptableIQR_hpb_sig_signs[[2]]
rm(toptableIQR_hpb_sig_signs)

GO_hpb_pos <- GO_analysis(rownames(toptableIQR_hpb_sig_pos), universe, "BP")
GO_hpb_neg <- GO_analysis(rownames(toptableIQR_hpb_sig_neg), universe, "BP")

jpeg("figuras/GOBP_hpb_pos.jpg", width = 900, height = 600)
dotplot(GO_hpb_pos, showCategory = 10, font.size = 10, title = "GO:BP Upregulation HPB-ALL")
dev.off()

jpeg("figuras/GOBP_hpb_neg.jpg", width = 900, height = 600)
dotplot(GO_hpb_neg, showCategory = 10, font.size = 10, title = "GO:BP Downregulation HPB-ALL")
dev.off()


#~~~~~~Intersect~~~~~~
# GO terms usando la interseccion de kopt y hpb: vademás de los genes
# de la intersección, necesitamos que tengan el mismo sentido en logFC

intersect_pos <- intersect(rownames(toptableIQR_kopt_sig_pos), rownames(toptableIQR_hpb_sig_pos)) #87
intersect_neg <- intersect(rownames(toptableIQR_kopt_sig_neg), rownames(toptableIQR_hpb_sig_neg)) #385
# Esto significa que de los 1110 DEGs de la intersección, sólo comparten el sentido de logFC 472.

GO_intersect_pos <- GO_analysis(intersect_pos, universe, "BP")
GO_intersect_neg <- GO_analysis(intersect_neg, universe, "BP")

jpeg("figuras/GOBP_intersect_pos.jpg", width = 900, height = 600)
dotplot(GO_intersect_pos, showCategory = 10, font.size = 10, title = "GO:BP Intersect Upregulation KOPT-K1, HPB-ALL")
dev.off()

jpeg("figuras/GOBP_intersect_neg.jpg", width = 900, height = 600)
dotplot(GO_intersect_neg, showCategory = 10, font.size = 10, title = "GO:BP Intersect Downregulation KOPT-K1, HPB-ALL")
dev.off()







######Preparing GSEA input files######

#Crea el directorio 'GSEA_files' en caso de que no exista
dir.create("./GSEA_files/", showWarnings = FALSE)

#~~~~~~.GCT~~~~~~
# Vamos a hacer un .gct con los datos de expresión normalizados de todas la muestras
# También hacemos un subset al expressionset por cell line hacer .gct de cada línea
eset_kopt <- eset[,c(1:3, 7:9)]
eset_hpb <- eset[,c(4:6, 10:12)]
#View(exprs(eset_hpb))
#View(exprs(eset_hpb))

output.gct(eset, filename = "./GSEA_files/GSE18198_base.gct")
output.gct(eset_kopt, filename = "./GSEA_files/GSE18198_kopt")
output.gct(eset_hpb, filename = "./GSEA_files/GSE18198_hpb")

#~~~~~~.RNK~~~~~~
# Primero vamos a crear el .rnk para el análisis por bloques:
fit3<-lmFit(eset, design_blocking)
fit4<-eBayes(fit3)
toptable_blocking<-topTable(fit4, coef="treat_blockingSAHM1", number=dim(exprs(eset))[1], adjust.method="BH", sort.by="p")
rm(fit3); rm(fit4)

# Vamos a crear el fichero .rnk comparando todos los genes para la comparación
# de líneas, para ello sacamos topTable con todos los genes
res_rnk <- deg(eset, design, cont.matrix_lines)
toptable_rnk <- res_rnk[[1]]
toptable_rnk_nsig <- res_rnk[[2]]
rm(res_rnk)



#Vamos a utilizar como métrica para el ranking el adjpval,
#pero teniendo en cuenta el sentido de la expresión diferencial

output.rnk <- function(toptable, file_path){
  x <- toptable[,c("adj.P.Val","logFC")]
  x$name<-rownames(toptable)
  x <- x[,c(3,1,2)]
  x$logP=-log10(x$adj.P.Val)
  x$sign<-sign(x$logFC)
  x$metric<-(x$logP/x$sign)
  
  y<-x[,c("name", "metric")]
  write.table(y, file = file_path,
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  return(y)
}

rnk_table <- output.rnk(toptable_rnk, "GSEA_files/GSE18198_cells.rnk")


#~~~~~~~~~.CLS~~~~~~~~~
#Ahora vamos a generar los .cls de cada linea celular:
treat_labels <- c(rep("DMSO", 6), rep("SAHM1", 6))
cell_labels <- c(rep("KOPT-K1", 3), rep("HPB-ALL", 3), rep("KOPT-K1", 3), rep("HPB-ALL", 3))
phenoData(eset)$treat <- treat_labels
phenoData(eset)$cell_line <- cell_labels
esetpd <- pData(eset)
esetpd_kopt <- esetpd[esetpd$cell_line=="KOPT-K1",]
esetpd_hpb <- esetpd[esetpd$cell_line=="HPB-ALL",]


output.cls(target = esetpd_kopt,
           variable = "treat",
           filename = "./GSEA_files/GSE18198_kopt")
output.cls(target = esetpd_kopt,
           variable = "treat",
           filename = "./GSEA_files/GSE18198_hpb")

#~~~~~~~~~.CHIP~~~~~~~~~
#Vamos a crear el fichero .chip. Primero, anotamos todos los genes,
#por ejemplo, usando toptable_rnk

toptable_rnk <- anno(toptable_rnk)

output.chip <- function(toptable, file_path){
  chip <- as.data.frame(toptable[,"SYMBOL"])
  colnames(chip)[1]<-"Gene Symbol"
  chip$"Probe Set ID"<-rownames(toptable)
  chip$"Gene Title"<-rep(NA, nrow(chip))
  chip <- chip[,c(2,1,3)]
  
  write.table(chip, file = file_path,
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  return(chip)
}

chip <- output.chip(toptable_rnk, "GSEA_files/GSE18198_base.chip")

