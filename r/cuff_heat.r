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

gene_breaks = (1:nrow(fpkm.matrix))-0.5

ggplot(df, aes(x=Sample, y=Gene, fill=FPKM)) +
 geom_tile() +
 scale_x_discrete("", limits=sample.order) + 
 scale_y_discrete(breaks=gene_breaks, limits=gene.order, labels=rep("",nrow(fpkm.matrix))) +
 scale_fill_gradient("FPKM", low="white", high="darkred") +
 theme(text = element_text(size=30), axis.text.x=element_text(size=30)) +
 theme_bw()

ggsave(output.pdf)
