library(ggplot2)
library(reshape2)
library(seriation)
library(RColorBrewer)

ca = commandArgs(trailing=T)
df.file = ca[1]
output.pdf = ca[2]

df = read.table(df.file, header=T, quote="\"")

fpkm.matrix = acast(df, Gene ~ Sample, value.var="FPKM")

gene.dist = dist(fpkm.matrix)
#gene.clust = hclust(gene.dist)
#gene.order = rownames(fpkm.matrix)[gene.clust$order]
gene.ser = seriate(gene.dist, method="OLO")
gene.order = rownames(fpkm.matrix)[get_order(gene.ser)]

sample.dist = dist(t(fpkm.matrix))
#sample.clust = hclust(sample.dist)
#sample.order = colnames(fpkm.matrix)[sample.clust$order]
sample.ser = seriate(sample.dist, method="OLO")
sample.order = colnames(fpkm.matrix)[get_order(sample.ser)]

ggplot(df, aes(x=Sample, y=Gene, fill=FPKM)) +
    geom_tile() +
    scale_x_discrete("", limits=sample.order) + 
    scale_y_discrete(limits=gene.order) +
    scale_fill_gradientn("FPKM", colours=c("white", brewer.pal(8, "YlOrRd"))) +
    theme_bw() +
    theme(text=element_text(size=20)) +
    theme(axis.text.x=element_text(angle=315, hjust=0, vjust=1), axis.ticks.y=element_blank(), axis.text.y=element_blank())

ggsave(output.pdf)

# scale_fill_gradientn("FPKM", colours=c("white", brewer.pal(8, "YOrRd"))) +
