library(ggplot2)
library(reshape2)

ca = commandArgs(trailing=T)
df.file = ca[1]
output.pdf = ca[2]

df = read.table(df.file, header=T, quote="\"")

# this is broken

df$dist = 1 - df$Correlation
sample12.matrix = acast(df, Sample1 ~ Sample2, value.var="dist", fill=0)
sample21.matrix = acast(df, Sample2 ~ Sample1, value.var="dist", fill=0)
sample.dist = as.dist(sample12.matrix+sample21.matrix)
sample.clust = hclust(sample.dist)
sample.order = rownames(sample12.matrix)[sample.clust$order]
sample.order = rownames(sample12.matrix)

ggplot(df, aes(x=Sample1, y=Sample2, fill=Correlation)) +
 geom_tile() +
 scale_x_discrete("", limits=sample.order) + 
 scale_y_discrete("", limits=sample.order) +
 scale_fill_gradient(low="white", high="darkred") +
 theme_bw()

ggsave(output.pdf)
