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
library(org.Hs.eg.db)

library(topGO)
library(clusterProfiler)
library(ggplot2)
library(stringr)
library(rrvgo)
library(pathview)


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

#colnames(ematrix) <- labels
#colnames(ematrix.kopt) <- labels.kopt
#colnames(ematrix.hpb) <- labels.hpb
rm(labels.hpb, labels.kopt)

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

get.pca.plots <- function(pcaobj){
  
  jpeg(paste("figuras/",as.character(deparse(substitute(pcaobj))),"_screeplot.jpg", sep = ""), width = 600, height = 600)
  print(screeplot(pcaobj))
  dev.off()

  jpeg(paste("figuras/",as.character(deparse(substitute(pcaobj))),"_plot.jpg", sep = ""), width = 600, height = 600)
  print(biplot(pcaobj, lab = labels))
  dev.off()
}

get.pca.plots(pca)
rm(pca)

#5. Generate BOXPLOTS before and after normalization

get.boxplots <- function(eset, data){
  
  jpeg(paste("figuras/boxplot_", as.character(deparse(substitute(data))), "_not_normalised.jpg", sep = ""), width = 600, height = 600)
  boxplot(data,
          main="Boxplot Before Normalization",
          col = "lightgrey",
          las=2)
  dev.off()
  
  exprseset <- as.data.frame(exprs(eset))
  
  jpeg(paste("figuras/boxplot_", as.character(deparse(substitute(eset))), "_normalised.jpg", sep = ""), width = 600, height = 600)
  boxplot(data.frame(exprseset),
          main="Boxplot After Normalization (log scale)",
          col = "white",
          las=2)
  dev.off()
}

get.boxplots(eset, data)
get.boxplots(eset.kopt, data.kopt)
get.boxplots(eset.hpb, data.hpb)
rm(data, data.hpb, data.kopt)



################ Differential expression analysis ##################
####################################################################

#7. Design matrix.
# Para el análisis con tipo celular como covariable:
cell.type.cov <- factor(targets$CellType)
treat.cov <- factor(targets$Treatment)
design.cov <- model.matrix(~cell.type.cov+treat.cov)
rownames(design.cov)<-targets$FileName

rm(cell.type.cov, treat.cov)

# Para el análisis de líneas celulares por separado y comparando ambas:
design.kopt <- cbind(KOPT_K1_DMSO=c(1,1,1,0,0,0),
                     KOPT_K1_SAHM1=c(0,0,0,1,1,1))
rownames(design.kopt) <- targets.kopt$FileName

design.hpb <- cbind(HPB_ALL_DMSO=c(1,1,1,0,0,0),
                    HPB_ALL_SAHM1=c(0,0,0,1,1,1))
rownames(design.hpb) <- targets.hpb$FileName


#8. Contrasts matrix.
#Sólo necesarias para las líneas celulares analizadas por separado
cont.matrix.kopt <- makeContrasts(KOPT_K1_SAHM1 - KOPT_K1_DMSO, levels=design.kopt)
cont.matrix.hpb <- makeContrasts(HPB_ALL_SAHM1 - HPB_ALL_DMSO, levels=design.hpb)


#9. Obtaining differentially expressed genes (DEGs)
# Primero para el analisis con covariable:
fit<-lmFit(esetIQR, design.cov)
fit<-eBayes(fit)
toptableIQR.cov<-topTable(fit, coef="treat.covSAHM1", number=dim(exprs(esetIQR))[1], adjust.method="BH", sort.by="p")
toptableIQR.cov.nsig <- nrow(toptableIQR.cov[toptableIQR.cov$adj.P.Val<0.05, ])
rm(fit)

#Ahora hacemos una función para las líneas por separado
generateToptable <- function(expression_set, design_matrix, contrasts_matrix){
  fit<-lmFit(expression_set,design_matrix) 
  fit<-contrasts.fit(fit, contrasts_matrix)
  fit<-eBayes(fit)
  toptable <- topTable(fit, number=dim(exprs(expression_set))[1],
                       adjust.method="BH", sort.by="p")
  significant_results <- nrow(toptable[toptable$adj.P.Val<0.05, ])
  
  results <- list(toptable, significant_results)
  return(results)
}

res.kopt <- generateToptable(esetIQR.kopt, design.kopt, cont.matrix.kopt)
toptableIQR.kopt <- res.kopt[[1]]
toptableIQR.kopt.nsig <- res.kopt[[2]]
rm(res.kopt)

res.hpb <- generateToptable(esetIQR.hpb, design.hpb, cont.matrix.hpb)
toptableIQR.hpb <- res.hpb[[1]]
toptableIQR.hpb.nsig <- res.hpb[[2]]
rm(res.hpb)

######Anotación######

anno <- function(toptable){
  affy_IDs <- as.character(row.names(toptable))
  
  GENE_SYMBOLS <- AnnotationDbi::select(hgu133plus2.db, keys=affy_IDs, columns="SYMBOL", keytype="PROBEID")
  ENSEMBL_IDs <-AnnotationDbi::select(hgu133plus2.db, keys=affy_IDs, columns="ENSEMBL", keytype="PROBEID")
  ENTREZ_IDs <-AnnotationDbi::select(hgu133plus2.db, keys=affy_IDs, columns="ENTREZID", keytype="PROBEID")
  
  
  GENE_SYMBOLS = GENE_SYMBOLS[!duplicated(GENE_SYMBOLS["PROBEID"]),]
  ENSEMBL_IDs = ENSEMBL_IDs[!duplicated(ENSEMBL_IDs["PROBEID"]),]
  ENTREZ_IDs = ENTREZ_IDs[!duplicated(ENTREZ_IDs["PROBEID"]),]
  
  
  toptable$SYMBOL <- GENE_SYMBOLS$SYMBOL
  toptable$ENSEMBL <- ENSEMBL_IDs$ENSEMBL
  toptable$ENTREZ <- ENTREZ_IDs$ENTREZID
  toptable$rank <- c(1:nrow(toptable))
  toptable<-toptable[, c(10,9,8,7,1,2,3,4,5,6)]
  
  return(toptable)
}

toptableIQR.cov <- anno(toptableIQR.cov)
toptableIQR.kopt <- anno(toptableIQR.kopt)
toptableIQR.hpb <- anno(toptableIQR.hpb)


# Comparar los analisis independientes kopt y hpb:
toptableIQR.kopt.sig <- toptableIQR.kopt[toptableIQR.kopt$adj.P.Val<=0.05,]
toptableIQR.hpb.sig <- toptableIQR.hpb[toptableIQR.hpb$adj.P.Val<=0.05,]
a <- rownames(toptableIQR.kopt.sig)
b <- rownames(toptableIQR.hpb.sig)

intersect_kopt_hpb <- intersect(a, b)
length(intersect_kopt_hpb) # prints '714'
rm(a); rm(b)
toptableIQR.hpb.nsig - length(intersect_kopt_hpb) # prints '1391'
toptableIQR.kopt.nsig - length(intersect_kopt_hpb) # prints '7051'



############# GeneOntology: BiologicalProcess (GO:BP) ##############
####################################################################

GO_analysis <- function (genes, universe, ontology){
  clusterProfiler::enrichGO(gene = genes,
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


#~~~~~~Covariable~~~~~~

toptableIQR.cov.sig <- toptableIQR.cov[toptableIQR.cov$adj.P.Val<=0.05,]
toptableIQR.cov.sig.signs <- FC_sign_subsetting(toptableIQR.cov.sig)
toptableIQR.cov.sig.pos <- toptableIQR.cov.sig.signs[[1]]
toptableIQR.cov.sig.neg <- toptableIQR.cov.sig.signs[[2]]
rm(toptableIQR.cov.sig.signs)

GO.cov.pos <- GO_analysis(rownames(toptableIQR.cov.sig.pos), universe, "BP")
GO.cov.pos.df <- as.data.frame(GO.cov.pos)
GO.cov.neg <- GO_analysis(rownames(toptableIQR.cov.sig.neg), universe, "BP")
GO.cov.neg.df <- as.data.frame(GO.cov.neg)

#~~~~~~KOPT-K1~~~~~~
toptableIQR.kopt.sig.signs <- FC_sign_subsetting(toptableIQR.kopt.sig)
toptableIQR.kopt.sig.pos <- toptableIQR.kopt.sig.signs[[1]]
toptableIQR.kopt.sig.neg <- toptableIQR.kopt.sig.signs[[2]]
rm(toptableIQR.kopt.sig.signs)

GO.kopt.pos <- GO_analysis(rownames(toptableIQR.kopt.sig.pos), universe, "BP")
GO.kopt.pos.df <- as.data.frame(GO.kopt.pos)
GO.kopt.neg <- GO_analysis(rownames(toptableIQR.kopt.sig.neg), universe, "BP")
GO.kopt.neg.df <- as.data.frame(GO.kopt.neg)

#~~~~~~HPB-ALL~~~~~~
toptableIQR.hpb.sig.signs <- FC_sign_subsetting(toptableIQR.hpb.sig)
toptableIQR.hpb.sig.pos <- toptableIQR.hpb.sig.signs[[1]]
toptableIQR.hpb.sig.neg <- toptableIQR.hpb.sig.signs[[2]]
rm(toptableIQR.hpb.sig.signs)

GO.hpb.pos <- GO_analysis(rownames(toptableIQR.hpb.sig.pos), universe, "BP")
GO.hpb.pos.df <- as.data.frame(GO.hpb.pos)
GO.hpb.neg <- GO_analysis(rownames(toptableIQR.hpb.sig.neg), universe, "BP")
GO.hpb.neg.df <- as.data.frame(GO.hpb.neg)


########################## GO clustering ############################
# calculateSimMatrix() calcula la similaridad semántica de los GO terms:

go_clustering <- function(GO_BP_df){
  # Matriz con la similitud entre todas las combinaciones de GO:BP terms
  simMatrix <- calculateSimMatrix(GO_BP_df$ID,
                                  orgdb="org.Hs.eg.db",
                                  ont="BP",
                                  method="Rel")
  
  # reduceSimMatrix() agrupa los GO terms en función de la similaridad
  # Primero vamos a asignar scores a los GO terms en función de su significancia estadística:
  scores <- setNames(-log10(GO_BP_df$qvalue), GO_BP_df$ID)
  # Ahora hacemos el clustering:
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Hs.eg.db")
  
  #Vamos a verlo representado:
  jpeg(paste("/home/guille/Desktop/Omics/microarrays/trabajo-microarray/figuras/GO_clusters/", as.character(deparse(substitute(GO_BP_df))), ".jpg", sep = ""),
       res = 300, width = 4000, height = 4000)
  print(treemapPlot(reducedTerms))
  dev.off()
}

go_clustering(GO.cov.pos.df)
go_clustering(GO.cov.neg.df)
go_clustering(GO.kopt.pos.df)
go_clustering(GO.kopt.neg.df)
go_clustering(GO.hpb.pos.df)
go_clustering(GO.hpb.neg.df)


########################### GSEA files #############################
####################################################################

#Crea el directorio 'GSEA_files' en caso de que no exista
dir.create("./GSEA_files/", showWarnings = FALSE)

#~~~~~~.GCT~~~~~~
# Vamos a hacer un .gct con los datos de expresión normalizados de todas las
# muestras y de cada línea por separado


output.gct(eset, filename = "./GSEA_files/GSE18198_base")
output.gct(eset.kopt, filename = "./GSEA_files/GSE18198_kopt")
output.gct(eset.hpb, filename = "./GSEA_files/GSE18198_hpb")

#~~~~~~.RNK~~~~~~
# Primero vamos a crear el .rnk para el análisis con covariable.
# Para gsea, usamos los datos sin filtrar
fit2 <- lmFit(eset, design.cov)
fit2 <- eBayes(fit2)
toptable.cov.gsea <-topTable(fit2, coef="treat.covSAHM1", number=dim(exprs(eset))[1], adjust.method="BH", sort.by="p")
rm(fit2)

# Vamos a utilizar como métrica para el ranking el adjpval,
# pero teniendo en cuenta el sentido de la expresión diferencial

output.rnk <- function(toptable, file_path){
  x <- toptable[,c("adj.P.Val","logFC")]
  x$name <- rownames(toptable)
  x <- x[,c(3,1,2)]
  x$logP = -log10(x$adj.P.Val)
  x$sign <- sign(x$logFC)
  x$metric <- (x$logP/x$sign)
  
  y<-x[,c("name", "metric")]
  write.table(y, file = file_path,
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  return(y)
}

output.rnk(toptable.cov.gsea, "GSEA_files/GSE18198_cov.rnk")


#~~~~~~~~~.CLS~~~~~~~~~

# Ahora vamos a generar los .cls de cada linea celular:
# En ambos casos tenemos 3 controles seguidos de 3 ttos, nos vale el mismo .cls
treat.labels <- c(rep("DMSO", 3), rep("SAHM1", 3))
phenoData(eset.kopt)$treat <- treat.labels
esetpd.gsea <- pData(eset.kopt)

output.cls(target = esetpd.gsea,
           variable = "treat",
           filename = "/home/guille/Desktop/Omics/microarrays/trabajo-microarray/GSEA_files/GSE18198_single_cell_line")

#~~~~~~~~~.CHIP~~~~~~~~~
#Vamos a crear el fichero .chip. Primero, anotamos todos los genes,
#por ejemplo, usando toptable_rnk

toptable.cov.gsea <- anno(toptable.cov.gsea)

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

output.chip(toptable.cov.gsea, "GSEA_files/GSE18198.chip")



######################### KEGG pathways ############################
####################################################################

# Buscamos pathways enriquecidas significativas

kegg.cov <- enrichKEGG(toptableIQR.cov$ENTREZ,
                       organism = "hsa",
                       keyType = "ncbi-geneid",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       universe = unique(toptable.cov.gsea$ENTREZ),
                       minGSSize = 10,
                       maxGSSize = 500,
                       qvalueCutoff = 0.2,
                       use_internal_data = FALSE) # usa la db online más reciente

kegg.cov.df <- as.data.frame(kegg.cov) # df con los resultados significativos


#~~~~~~~~~~~~~~~~~~~~~~~~~~~Pathview~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("/home/guille/Desktop/Omics/microarrays/trabajo-microarray/kegg/")

# Pathview necesita una matriz numerica con los entrezids como rownames
# y las muestras como colnames.

# Tomamos los datos de la tabla de gsea, que es la que no ha sido filtrada
pathview_df.cov <- toptable.cov.gsea[, c("ENTREZ", "logFC", "adj.P.Val")]
pathview_df.cov <- pathview_df.cov[!is.na(pathview_df.cov$ENTREZ),]
pathview_df.cov <- pathview_df.cov[!duplicated(pathview_df.cov$ENTREZ),]
pathview_df.cov$sig <- NA # creamos una columna vacía que luego rellenaremos

# Quitamos esta pathway porque da error con pathview:
kegg.cov.df <- kegg.cov.df[-which(kegg.cov.df$ID=="hsa05206"),]

# Esta función crea una columna adicional que pathview leerá como una segunda
# muestra, pero en realidad nos va a indicar la significancia (valor rojo, >0)
# o la no significancia (valor blanco, ==0).
# Además, el valor positivo para los significantes lo toma del mayor logFC en
# valor absoluto del conjunto de genes que forman cada una de las rutas, para
# evitar modificar el rango de expresión de dichos genes y mantener la escala
# de colores original.

kegg_exprs_sig <- function(kegg_df){
  
for(i in 1:length(kegg_df[,"geneID"])){
  print("#########GENERATING PATHWAY##########")
  a <- unlist(strsplit(kegg_df[i,"geneID"], split = "/"))
  
  # Sacar el máximo logFC de pathview_df.cov de cada grupo de genes:
  maxi <- -Inf
  for(j in a){
    if(abs(pathview_df.cov[pathview_df.cov["ENTREZ"]==j,"logFC"]) > maxi){
      maxi <- abs(pathview_df.cov[pathview_df.cov["ENTREZ"]==j,"logFC"])}
  }
  
  # Si el padj es significativo, entonces le mete maxi, si no, un 0
  for(k in 1:length(pathview_df.cov$adj.P.Val)){
    if(pathview_df.cov$adj.P.Val[k] < 0.05){
      pathview_df.cov$sig[k] <- maxi
    }
    else{
      pathview_df.cov$sig[k] <- 0
    }
  }
  
  # Convertir a matrix y darle el formato que requiere pathview
  pathview_matrix.cov <- as.matrix(pathview_df.cov)
  pathview_matrix.cov <- pathview_matrix.cov[,-c(1,3)]
  pathview_matrix.cov <- apply(pathview_matrix.cov, 2, as.numeric)
  row.names(pathview_matrix.cov) <- pathview_df.cov$ENTREZ
  
  # Sacamos la pathway:
  pathview(gene.data = pathview_matrix.cov,
           pathway.id = kegg_df[i,"ID"],
           species = "hsa",
           gene.idtype = "entrez",
           kegg.dir = "/home/guille/Desktop/Omics/microarrays/trabajo-microarray/kegg/kegg_data/",
           low = "green", mid = "white", high = "red", na.col = "grey",
           multi.state = T, same.layer = T)
}
rm(a, i, j, k, maxi)
}

# Corremos la función con las pathways significativas de enrichKEGG()
kegg_exprs_sig(kegg.cov.df)


# Manualmente miramos rutas reguladas por NOTCH1 o importantes en cáncer.
setwd("/home/guille/Desktop/Omics/microarrays/trabajo-microarray/kegg/others/")

manpaths <- c("hsa04151", #PI3K-AKT
              "hsa04010", #MAPK
              "hsa05230", #Central carbon metabolism in cancer
              "hsa04612", #Antigen processing and presentation
              "hsa05235", #PD-L1 expression and PD-1 checkpoint pathway in cancer
              "hsa04110", #Cell cycle
              "hsa04210") #Apoptosis

kegg.cov.others.df <- kegg.cov@result[kegg.cov@result[,"ID"] %in% manpaths, ]
kegg_exprs_sig(kegg.cov.others.df)

# Estas rutas no aparecen ni siquiera como no significativas en los resultados,
# pero visto los resultados de GO:BP, decidimos mirarlas
setwd("/home/guille/Desktop/Omics/microarrays/trabajo-microarray/kegg/others2/all_genes/")

manpaths2 <- c("hsa00230", #Purine metabolism
               "hsa00240", #Pirimidine metabolism
               "hsa00190", #Oxidative phosphorylation
               "hsa03030") #DNA replication 
               
# a) Con todos los genes, significantivos y no significantes
for(i in manpaths2){
  pathview(gene.data = pathview_matrix.cov[,1],
           pathway.id = i,
           species = "hsa",
           gene.idtype = "entrez",
           kegg.dir = "/home/guille/Desktop/Omics/microarrays/trabajo-microarray/kegg/kegg_data/",
           low = "green", mid = "white", high = "red", na.col = "grey")
}

# b) Sólo con los genes significativos
setwd("/home/guille/Desktop/Omics/microarrays/trabajo-microarray/kegg/others2/significant_genes/")

manpaths2.sig <- toptable.cov.gsea[toptable.cov.gsea$adj.P.Val<0.05,"logFC"]
names(manpaths2.sig) <- toptable.cov.gsea[toptable.cov.gsea$adj.P.Val<0.05,"ENTREZ"]

for(i in manpaths2){
  pathview(gene.data = pathview_matrix.cov[,1],
           pathway.id = i,
           species = "hsa",
           gene.idtype = "entrez",
           kegg.dir = "/home/guille/Desktop/Omics/microarrays/trabajo-microarray/kegg/kegg_data/",
           low = "green", mid = "white", high = "red", na.col = "grey")
}
