### This code primarily aims to demonstrate the inter-group differential expression of miRNAs. 
### It takes the miRNA alignment results in the form of a CSV file, obtained after processing with miRDeep2, as input. 
### It outputs a correlation heatmap and trend lines.

library(ggplot2)
library(Seurat)
unloadNamespace('Seurat')
library(ggrepel)
mir_raw <- dir()[grep('expre',dir())]
mir_raw <- lapply(mir_raw, function(x){
  x <- read.csv(x,sep='\t')
  x <- split(x,x$X.miRNA)
  x <- lapply(x, function(df){df[1,]})
  x <- do.call(rbind,x)
})
mtx <- do.call(cbind,mir_raw)
mtx <- mtx[,grep("read_count",colnames(mtx))]
mtx <- as.data.frame(mtx)
mtx$sum <- rowSums(mtx)
mtx <- mtx[-which(mtx$sum==0),]
frac <- log2(t(t(mtx)/colSums(mtx)) * 1000000 + 0.25)


std <- apply(frac, 1, sd)
mn <- apply(frac, 1, mean)

dat <- as.data.frame(cbind(std,mn))
dat$miRNA <- rownames(dat)
#dat$miRNA <- gsub('\\.','-',dat$miRNA)
p <- ggplot(data = dat,aes(x=mn,y=std))+
  geom_point(color='#393e46',size=2)+
  geom_smooth(method = "loess", se=FALSE, linewidth=2,
              color="#e84545", formula = y ~ x)+
  xlab('Log2(miRNA level)')+
  ylab('Std Log2(miRNA level)')

p + geom_label_repel(data = subset(dat, dat$std > 5.5 | dat$mn > 16 |dat$mn > 14 &dat$std >4),
                   aes(label = miRNA), fill = NA, label.size = NA,size = 6,
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.8, "lines"), 
                   segment.color = "black", show.legend = FALSE,
                   segment.alpha=1)+
  theme_classic(base_size = 16)
                 
