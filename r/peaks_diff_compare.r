library(ggplot2)

ca = commandArgs(trailing=T)
df.file = ca[1]
output.pdf = ca[2]

df = read.table(df.file, header=T, quote="\"")

x.min = quantile(df$Test_stat, .005, na.rm=T)
x.max = quantile(df$Test_stat, .995, na.rm=T)
x.min = min(x.min, -x.max)
x.max = max(-x.min, x.max)

ggplot(df, aes(x=Test_stat, colour=Peak)) +
 geom_density() +
 scale_x_continuous(limits=c(x.min,x.max)) +
 theme_bw()

ggsave(output.pdf)
