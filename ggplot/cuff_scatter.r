library(ggplot2)

ca = commandArgs(trailing=T)
df.file = ca[1]
output.pdf = ca[2]
cond1 = ca[3]
cond2 = ca[4]

df = read.table(df.file, header=T, quote="\"")

ggplot(df, aes(x=fpkm1, y=fpkm2, colour=qval)) +
 geom_point(size=1.5, alpha=.3) +
 scale_x_continuous(paste(cond1, "log2 FPKM")) +
 scale_y_continuous(paste(cond2, "log2 FPKM")) +
 geom_abline(intercept=0, slope=1, linetype=2) +
 ggtitle(paste(cond1, "vs", cond2)) +
 theme_bw()

ggsave(output.pdf)
