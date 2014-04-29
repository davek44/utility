library(ggplot2)

ca = commandArgs(trailing=T)
df.file = ca[1]
output.pdf = ca[2]
cond1 = ca[3]
cond2 = ca[4]

df = read.table(df.file, header=T, quote="\"")

x.min = min(df$fpkm1)
x.max = quantile(df$fpkm2,.997)

y.min = min(df$fpkm2)
y.max = quantile(df$fpkm2, .997)

qval.unique = unique(df$qval)
if (length(qval.unique) == 1) {
    gp = ggplot(df, aes(x=fpkm1, y=fpkm2))
} else {
    gp = ggplot(df, aes(x=fpkm1, y=fpkm2, colour=qval))
}

gp +
    geom_point(size=1.5, alpha=.3) +
    scale_x_continuous(paste(cond1, "log2 FPKM"), lim=c(x.min,x.max)) +
    scale_y_continuous(paste(cond2, "log2 FPKM"), lim=c(y.min,y.max)) +
    geom_abline(intercept=0, slope=1, linetype=2) +
    ggtitle(paste(cond1, "vs", cond2)) +
    theme_bw() +
    theme(text=element_text(size=16))

ggsave(output.pdf)
