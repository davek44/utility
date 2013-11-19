library(ggplot2)

ca = commandArgs(trailing=T)
df.file = ca[1]
plot.title = ca[2]
output.pdf = ca[3]

df = read.table(df.file, header=T, quote="\"")

annotation.order.all = c('Intergenic','Introns','3\'UTR','5\'UTR','CDS','lncRNA','Pseudogene','rRNA','smallRNA')
annotation.order = annotation.order.all[annotation.order.all %in% df$annotation]
df$annotation = factor(df$annotation, levels=annotation.order)

ggplot(df, aes(x=annotation, y=ratio)) +
 geom_bar(stat="identity") +
 scale_x_discrete("Annotation") +
 scale_y_continuous("log2 feature% / length%") +
 ggtitle(plot.title) +
 theme_bw() +
 theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())

ggsave(output.pdf)
