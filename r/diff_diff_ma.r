library(ggplot2)

ca = commandArgs(trailing=T)
df.file = ca[1]
output.pdf = ca[2]

df = read.table(df.file, header=T, quote="\"")

x.min = quantile(df$avg, .002)
x.max = quantile(df$avg, .998)

y.min = quantile(df$minus, .002)
y.max = quantile(df$minus, .998)

ggplot(df, aes(x=avg, y=minus)) +
    geom_point(size=1.5, alpha=.3) +
    scale_x_continuous("Avg test stat", lim=c(x.min,x.max)) +
    scale_y_continuous("Test stat 1 - 2", lim=c(y.min,y.max)) +
    geom_hline(y=0) +
    theme_bw() +
    theme(text=element_text(size=16))

ggsave(output.pdf)
