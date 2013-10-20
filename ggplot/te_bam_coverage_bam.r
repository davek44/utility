library(ggplot2)

ca = commandArgs(trailing=T)
df.file = ca[1]
plot.title = ca[2]
output.pdf = ca[3]

df = read.table(df.file, header=T)

ggplot(df, aes(x=indexes, y=coverage)) +
 geom_histogram(stat="identity") +
 scale_x_continuous("TE index") +
 scale_y_continuous("") +
 ggtitle(plot.title) +
 theme_bw()

ggsave(output.pdf)
