library(ggplot2)
library(reshape2)

ca = commandArgs(trailing=T)
df.file = ca[1]
output.pdf = ca[2]

df = read.table(df.file, header=T, quote="\"")

fpkm.matrix = acast(df, Gene ~ Sample, value.var="FPKM")

gene.dist = dist(fpkm.matrix)
gene.clust = hclust(gene.dist)
gene.order = rownames(fpkm.matrix)[gene.clust$order]

sample.dist = dist(t(fpkm.matrix))
sample.clust = hclust(sample.dist)
sample.order = colnames(fpkm.matrix)[sample.clust$order]

ggplot(df, aes(x=Sample, y=Gene, fill=FPKM)) +
    geom_tile() +
    scale_x_discrete("", limits=sample.order) + 
    scale_y_discrete(limits=gene.order) +
    scale_fill_gradient("FPKM", low="white", high="tomato3") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=315, hjust=0, vjust=1), axis.ticks.y=element_blank(), axis.text.y=element_blank())

ggsave(output.pdf)
