### The functionality of this code is to generate a volcano plot depicting the miRNA expression patterns in the BF and MQ datasets. 
### The input for this code is a compiled collection of miRNA expression data from the two datasets (in CSV format). 
### The output will be volcano plots illustrating the differential expression of miRNAs in these two datasets.

rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(EnhancedVolcano)

data_sample <- read.csv(file="BFfinal.csv", header= T)
data_sample <- data_sample[!duplicated(data_sample$miRNA),]
rownames(data_sample) <- data_sample$miRNA
data_sample <- subset(data_sample, select = -miRNA)
data <- (data_sample)

### Create Seurat object for each dataset
BF <- CreateSeuratObject(counts = data,
                          project = "BF", 
                          min.cells = 8,
                          min.features = 0) 


data_sample1 <- read.csv(file="MQfinal.csv", header= T)

data_sample1 <- data_sample1[!duplicated(data_sample1$miRNA),]
rownames(data_sample1) <- data_sample1$miRNA
data_sample1 <- subset(data_sample1, select = -miRNA)
data1 <- (data_sample1)

MQ <- CreateSeuratObject(counts = data1,
                          project = "MQ", 
                          min.cells = 8, 
                          min.features = 0) 

mgall.combined <- merge(BF, y = MQ, add.cell.ids = c("BF", "MQ"), project = "MGALL")

table(mgall.combined$orig.ident)

### normalization
mgall.combined <- NormalizeData(mgall.combined, normalization.method = "LogNormalize", scale.factor = 1e4)
mgall.combined <- FindVariableFeatures(mgall.combined, selection.method = "vst", nfeatures = 2000) #返回两千个高变基

markers_2 <- FindMarkers(mgall.combined, ident.1="BF", ident.2="MQ")
head(x = markers_2)

### create figure
EnhancedVolcano(markers_2 ,
                lab = rownames(markers_2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'CART Vol',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                subtitle = NULL,
                caption = NULL,
                legendPosition = 'bottom')+ coord_flip()

pdf('vol.pdf')



rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(EnhancedVolcano)


### Alternatively, the differential expression results can be directly used for further analysis. 
### The input data source would be the differential expression results obtained through the DESeq2 program (in a CSV file format).

data_sample <- read.table(file="/home/disk/naxing/project5/2.18/DESeq2.txt", header= T)

data <- as.data.frame(lapply(as.data.frame(data_sample),as.numeric))

pdf("bench_query_sort.pdf", height = 9, width = 10, onefile = F)  
EnhancedVolcano(data ,
                lab = rownames(data_sample),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'CART Vol',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 4.0,
                subtitle = NULL,
                caption = NULL,
                widthConnectors = 1.0,
                col = c('#6b6b6b', '#558fe6', 'purple', '#dd2e2e'),
                colAlpha = 4/5,
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 4.0,
                colConnectors = 'black')+ coord_flip()
dev.off()
