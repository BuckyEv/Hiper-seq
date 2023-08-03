### The main functionality of this code is to perform differential gene expression analysis on two datasets, BF and MQ, and visualize the expression patterns of differentially expressed genes in a heatmap. 
### The input for this code is the miRNA alignment results obtained after processing with miRDeep2 in the form of a CSV file. 
### The output includes the differential expression results of genes and a heatmap depicting the expression patterns of differentially expressed genes.

library(DESeq2)

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
  samplefile <- gsub('_L._1_5M_trimmed.csv','',samplefile)
  colnames(mtx) <- samplefile
  return(mtx)
}

BF <- clean_data('BF')
MQ <- clean_data('MQ')

mtx <- as.data.frame(cbind(BF, MQ))
mtx$sum <- rowSums(mtx)
mtx <- mtx[which(!mtx$sum == 0),]
mtx$sum <- NULL

coldata <- data.frame(condition = factor(rep(c('BF', 'MQ'), c(20,20)), levels = c('BF', 'MQ')))
dds <- DESeqDataSetFromMatrix(countData = mtx, colData = coldata, design= ~condition)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 3, parallel = FALSE)
res <- results(dds1, contrast = c('condition', 'BF', 'MQ'))
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)

### Filtering differentially expressed genes
res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.05),'sig'] <- 'up'
res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.05),'sig'] <- 'down'
res1[which(abs(res1$log2FoldChange) <= 1 | res1$padj >= 0.05),'sig'] <- 'none'

### Output the final list of selected differentially expressed genes
res1_select <- subset(res1, sig %in% c('up', 'down'))
res_up <- subset(res1, sig == 'up')
res_down <- subset(res1, sig == 'down')
write.table('SSC_up_DE_gene.txt',x=rownames(res_up), quote = F, col.names = F, row.names = F)
write.table('SSC_down_DE_gene.txt',x=rownames(res_down), quote = F, col.names = F, row.names = F)
write.table(res1_select, file = 'DESeq2.select.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(res1, file = 'DESeq2.txt', sep = '\t', col.names = NA, quote = FALSE)

### create figure
library(pheatmap)

de_mirna <- rownames(res1_select)

mtx <- log2(t(t(mtx)/colSums(mtx)) * 1000000 + 0.25)
mtx <- mtx[which(rownames(mtx) %in% de_mirna),]
mtx <- t(scale(t(mtx)))

### breaks
bk <- c(seq(-4,-0.1,by=0.01),seq(0,4,by=0.01))

group_sample <- as.data.frame(do.call(rbind,strsplit(colnames(mtx),'-'))[,1])
rownames(group_sample) <- colnames(mtx)
colnames(group_sample) <- 'Cell line'
colors=list(Group=c(BF="#7882A4", MQ="#C0A080"))

pheatmap(mtx,
         scale="none",
         color = c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)), 
         legend_breaks=seq(-8,8,2),
         border=FALSE,
         cellwidth = 15, 
         cellheight = 12,
         fontsize = 8,
         filename = "test-0.01.pdf",
         cluster_cols = FALSE,
         annotation_col=group_sample, 
         annotation_colors=colors,
         breaks=bk)

