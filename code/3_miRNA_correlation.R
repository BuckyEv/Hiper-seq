### The main functionality of this code is to demonstrate the correlation between bulk data and single-cell data, aiming to validate the similarity between the results from single-cell analysis and multi-cell analysis. 
### The input for this code is also the miRNA alignment results in the form of a CSV file obtained after processing with miRDeep2. 
### The output includes a correlation heatmap, trend lines, and the correlation coefficient (R) representing the degree of correlation.

library(ggplot2)

bulk <- dir('bulk/')
sc <- dir('sc')
bulk <- read.csv(paste0('bulk/',bulk),sep='\t')[,1:3]
bulk <- split(bulk,bulk$X.miRNA)
bulk <- lapply(bulk, function(df){df[1,]})
bulk <- do.call(rbind,bulk)
 
sc <- lapply(sc, function(x){
  x <- paste0('sc/',x)
  x <- read.csv(x,sep='\t')
  x <- split(x,x$X.miRNA)
  x <- lapply(x, function(df){df[1,]})
  x <- do.call(rbind,x)
})

sc <- do.call(cbind,sc)

sc_mtx <- sc[,grep('read_count',colnames(sc))]
bulk_mtx <- bulk[,grep('read_count',colnames(bulk))]

sc_b <- rowSums(sc_mtx) 
sc_bulk <- cbind(bulk_mtx,sc_b)

mtx <- sc_bulk[-which(sc_bulk[,1] == 0|sc_bulk[,2] == 0),]

rpm <- t(t(mtx)/colSums(mtx) * 1000000)

rpm[,1] <- log(rpm[,1]+1)
rpm[,2] <- log(rpm[,2]+1)


rpm <- as.data.frame(rpm)
colnames(rpm) <- c('bulk','sc')

cor.test(rpm[,1],rpm[,2])

p <- ggplot(rpm, aes(x=bulk,y=sc))+geom_point(color='#393e46',size=2)+geom_smooth(linewidth=2,color="#e84545",method="loess",se=F)
p <- p+theme_bw()+theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(color="grey50"))
p + annotate(geom="text", x=11, y=8, label="cor = 0.8875139") + labs(x = 'miRNA level in bulk (RPM)', y = 'miRNA level in chip (RPM)')+theme(axis.line.y=element_line(linetype=1,color="black",linewidth=0.8))+ theme(axis.line.x=element_line(linetype=1,color="black",linewidth=0.8))+ theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))
