library(ggplot2)

ca = commandArgs(trailing=T)
df.file = ca[1]
output.pdf = ca[2]

df = read.table(df.file, header=T, quote="\"")

annotation.order.all = c('Intergenic','Introns','Exons','3\'UTR','5\'UTR','CDS','lncRNA','Pseudogene','rRNA','smallRNA')
annotation.order = annotation.order.all[annotation.order.all %in% df$annotation]
df$annotation = factor(df$annotation, levels=annotation.order)

ggplot(df, aes(x=dummy, y=count, fill=annotation)) +
 geom_bar(stat="identity", width=1) +
 coord_polar(theta="y") +
 scale_x_discrete("") +
 scale_y_continuous("") +
 scale_fill_discrete("Annotation") +
 theme_bw() +  
 theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank())

ggsave(output.pdf)
