library(ggplot2)

ca = commandArgs(trailing=T)
df.file = ca[1]
output.pdf = ca[2]

df = read.table(df.file, header=T, quote="\"")

ggplot(df, aes(x=Index, y=Feature, fill=Coverage)) +
 geom_tile() +
 scale_fill_gradient(low="white", high="darkred") +
 theme_bw()

ggsave(output.pdf)
