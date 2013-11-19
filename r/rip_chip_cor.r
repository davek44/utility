library(ggplot2)
library(reshape2)

ca = commandArgs(trailing=T)
df.file = ca[1]
output.pdf = ca[2]

df = read.table(df.file, header=T, quote="\"")

cor.matrix = acast(df, ChIP ~ RIP, value.var="Cor")

chip.dist = dist(cor.matrix)
chip.clust = hclust(chip.dist)
chip.order = rownames(cor.matrix)[chip.clust$order]

rip.dist = dist(t(cor.matrix))
rip.clust = hclust(rip.dist)
rip.order = colnames(cor.matrix)[rip.clust$order]

ggplot(df, aes(x=RIP, y=ChIP, fill=Cor)) +
 geom_tile() +
 scale_x_discrete(limits=rip.order) +
 scale_y_discrete(limits=chip.order) +
 scale_fill_gradient2(low = "blue", high = "red") +
 theme_bw()

ggsave(output.pdf)
