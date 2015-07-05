library(ggplot2)

ca = commandArgs(trailing=T)
df.file = ca[1]
output.pdf = ca[2]

df = read.table(df.file, header=T, quote="\"")

x.min = quantile(df$diff1, .002)
x.max = quantile(df$diff1, .998)

y.min = quantile(df$diff2, .002)
y.max = quantile(df$diff2, .998)

ggplot(df, aes(x=diff1, y=diff2)) +
    geom_point(size=1.5, alpha=.3) +
    stat_smooth(method="lm") +
    scale_x_continuous("Test stat 1", lim=c(x.min,x.max)) +
    scale_y_continuous("Test stat 2", lim=c(y.min,y.max)) +
    geom_abline(intercept=0, slope=1, linetype=2) +
    theme_bw() +
    theme(text=element_text(size=25))

# stat_smooth(method="lm") +

ggsave(output.pdf)
