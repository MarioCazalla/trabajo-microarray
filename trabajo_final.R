### Establecer directorio de trabajo y cargar librerias
setwd("/home/guille/Desktop/Omics/microarrays/Trabajo_final/")

library("affy")
library("limma")
library("genefilter")

#2. Import targets.txt file. generates a data frame:
targets <- readTargets("targets.txt", row.names="FileName")

#3. Import .CEL files
data <- ReadAffy(filenames=targets$FileName)

#4. Normalize with RMA 
eset <- expresso(data,
                 bg.correct = TRUE, 
                 bgcorrect.method="rma",
                 normalize = TRUE, 
                 normalize.method="quantiles", 
                 pmcorrect.method="pmonly", 
                 summary.method="medianpolish",
                 verbose = TRUE,
)


#PCA analysis:
library(PCAtools)
ematrix <- exprs(eset)
pca1 <- pca(ematrix, metadata = ematrix_pd)

labels <- c()
for (label in c("KOPT_K1-C", "HPB_ALL-C", "KOPT_K1-T", "HPB_ALL-T" )){
  labels <- c(labels, rep(label,3))
}

dir.create("./figuras/", showWarnings = FALSE)

jpeg("figuras/pca_screeplot.jpg", width = 600, height = 600)
screeplot(pca1)
dev.off()

jpeg("figuras/pca_plot.jpg", width = 600, height = 600)
biplot(pca1, lab = labels)
dev.off()



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


#6.Data filtering using IQR.
esetIQR <- varFilter(eset, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)



#######Differential expression analysis.#######

#7. Design matrix.
cell_type <- factor(targets$CellType)
treat <- factor(targets$Treatment)
design <- model.matrix(~cell_type+treat)
rownames(design)<-targets$FileName


#9. Obtaining differentially expressed genes (DEGs)
#Linear model and eBayes 
fit<-lmFit(esetIQR, design)
fit2<-eBayes(fit)

#Table with DEGs results.
toptableIQR<-topTable(fit2, coef="treatSAHM1", number=dim(exprs(esetIQR))[1], adjust.method="BH", sort.by="p")
nrow(toptableIQR[toptableIQR$adj.P.Val<0.05, ])

##10. Save results
save(toptableIQR,file="MyResults.RData")


library(Biobase)
library(GEOquery)
library(annotate)
library(hgu133plus2.db)

#Cargamos el dataframe que hemos sacado antes (opcional si no lo tenemos ya):
load("/home/guille/Desktop/Omics/microarrays/Trabajo_final/MyResults.RData")

# Sacamos los IDs de affymetrix
affy_IDs <- as.character(row.names(toptableIQR))
# Convertimos de ID affymetrix a GENE SYMBOL
GENE_SYMBOLS <- select(hgu133plus2.db, keys=affy_IDs, columns="SYMBOL", keytype="PROBEID")
ENSEMBL_IDs <- select(hgu133plus2.db, keys=affy_IDs, columns="ENSEMBL", keytype="PROBEID")
# Eliminamos duplicados por multiple mapping generados en la conversion
GENE_SYMBOLS = GENE_SYMBOLS[!duplicated(GENE_SYMBOLS["PROBEID"]),]
ENSEMBL_IDs = ENSEMBL_IDs[!duplicated(ENSEMBL_IDs["PROBEID"]),]
# Añadimos las conversiones al dataframe
toptableIQR$SYMBOL <- GENE_SYMBOLS$SYMBOL
toptableIQR$ENSEMBL <- ENSEMBL_IDs$ENSEMBL
toptableIQR$rank <- c(1:nrow(toptableIQR))
# Reordenamos las columnas
toptableIQR<-toptableIQR[, c(9,7,8,1,2,3,4,5,6)]
View(toptableIQR)




######Preparing GSEA input files######

#Crea el directorio 'GSEA_files' en caso de que no exista
dir.create("./GSEA_files/", showWarnings = FALSE)

#Añadimos a phenoData que muestras son controles y cuales tratamientos
treat_labels <- c(rep("DMSO", 6), rep("SAHM1", 6))
phenoData(eset)$treat <- treat_labels
esetpd <- pData(eset)

#Utilizamos las funciones output.gct() y output.cls() del paquete ArrayTools
#utilizamos eset porque es un ExpressionSet normalizado

library(ArrayTools)
output.gct(eset, filename = "./GSEA_files/GSE18198_base")
output.cls(target = esetpd,
           variable = "treat",
           filename = "./GSEA_files/GSE18198_base")

# Vamos a crear el fichero .rnk, para ello hacemos topTable con todos los genes
fit3<-lmFit(eset, design)
fit4<-eBayes(fit3)
toptable_all<-topTable(fit4, coef="treatSAHM1", number=dim(exprs(eset))[1], adjust.method="BH", sort.by="p")

affy_IDs_all <- as.character(row.names(toptable_all))
GENE_SYMBOLS_all <- select(hgu133plus2.db, keys=affy_IDs_all, columns="SYMBOL", keytype="PROBEID")
ENSEMBL_IDs_all <- select(hgu133plus2.db, keys=affy_IDs_all, columns="ENSEMBL", keytype="PROBEID")

GENE_SYMBOLS_all = GENE_SYMBOLS_all[!duplicated(GENE_SYMBOLS_all["PROBEID"]),]
ENSEMBL_IDs_all = ENSEMBL_IDs_all[!duplicated(ENSEMBL_IDs_all["PROBEID"]),]

toptable_all$SYMBOL <- GENE_SYMBOLS_all$SYMBOL
toptable_all$ENSEMBL <- ENSEMBL_IDs_all$ENSEMBL
toptable_all$rank <- c(1:nrow(toptable_all))

toptable_all<-toptable_all[, c(9,7,8,1,2,3,4,5,6)]
View(toptable_all)

#Vamos a utilizar como métrica para el ranking el adjpval,
#pero teniendo en cuenta el sentido de la expresión diferencial

x<-data.frame(adj.P.Val = numeric(54675), logFC = numeric(54675))
x$name<-rownames(toptable_all)
x$adj.P.Val <- toptable_all$adj.P.Val
x$logFC <- toptable_all$logFC
x$logP=-log10(x$adj.P.Val)
x$sign<-sign(x$logFC)
x$metric<-(x$logP/x$sign)

y<-x[,c("name", "metric")]
write.table(y, file = "GSEA_files/GSE18198_base.rnk",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#Vamos a crear el fichero .chip
chip <- as.data.frame(toptable_all[,"SYMBOL"])
colnames(chip)[1]<-"Gene Symbol"
chip$"Probe Set ID"<-rownames(toptable_all)
chip$"Gene Title"<-rep(NA, nrow(chip))
chip <- chip[,c(2,1,3)]

write.table(chip, file = "GSEA_files/GSE18198_base.chip",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

