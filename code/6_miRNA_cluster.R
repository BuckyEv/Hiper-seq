### The primary function of this code is to perform unsupervised clustering of miRNA expression patterns. 
### The input for this code is a compiled collection of miRNA expression data from the two datasets (in CSV format). 
### The output will be the clustering results for these two datasets.

clean_data <- function(sample){
  samplename <- dir(sample)
  samplefile <- samplename[grep('expre',samplename)]
  mtx <- lapply(samplefile, function(x){
    x <- read.csv(paste0(sample, '/', x),sep='\t')
    x <- split(x,x$X.miRNA)
    x <- lapply(x, function(df){df[1,]})
    x <- do.call(rbind,x)
  })
  mtx <- do.call(cbind, mtx)
  mtx <- mtx[,grep("read_count",colnames(mtx))]
  mtx <- as.data.frame(mtx)
  samplefile <- gsub('miRNAs_expressed_all_samples_','',samplefile)
  samplefile <- gsub('_1_5M_trimmed.csv','',samplefile)
  colnames(mtx) <- samplefile
  return(mtx)
}

BF <- clean_data('BF')
MQ <- clean_data('MQ')



mtx <- as.data.frame(cbind(BF,MQ))

#colnames(mtx) <- c(paste('H2228',1:10,sep='_'),paste('H1975',1:10,sep='_'),paste('A549',1:10,sep='_'))
library(Seurat)

srat <- CreateSeuratObject(mtx,min.cells = 5)
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 20000)

plot1 <- VariableFeaturePlot(srat)
top10 <- head(VariableFeatures(srat), 10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes <- rownames(x = srat)
srat <- ScaleData(object = srat, features = all.genes)
srat <- RunPCA(srat, features = VariableFeatures(srat),npcs = 10)

srat <- RunTSNE(object = srat, dims.use = 1:10, perplexity =1)
DimPlot(object = srat, reduction = "tsne")

srat <- RunUMAP(srat, reduction = "pca", dims = 1:5, n.neighbors = 14)
DimPlot(srat, repel = TRUE, label.size = 3, reduction = "umap")

library(ggplot2)
UMAP <- srat@reductions[["umap"]]@cell.embeddings
UMAP <- as.data.frame(UMAP)
UMAP$Species <- gsub('-.*','',rownames(UMAP))
p <- ggplot(UMAP, aes(x = UMAP_1, y = UMAP_2 ,color=Species)) +
  geom_point() +
  stat_ellipse(level = 0.9)+
  theme(axis.line.x=element_line(linetype=1,color="black",size=3))
p + labs(color = 'Color') + theme_bw() + 
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line=element_line(color="grey50")) +
  geom_text(aes(x = UMAP_1, y = UMAP_2, label = rownames(UMAP)))


TSNE <- srat@reductions[["tsne"]]@cell.embeddings
TSNE <- as.data.frame(TSNE)
TSNE$Species <- gsub('-.*','',rownames(TSNE))
p1 <- ggplot(TSNE, aes(x = tSNE_1, y = tSNE_2 ,color=Species)) +
  geom_point() +
  stat_ellipse(level = 0.9)
p1 + labs(color = 'Color') + theme_bw() + 
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line=element_line(color="grey50")) +
  geom_text(aes(x = tSNE_1, y = tSNE_2, label = rownames(TSNE)))

