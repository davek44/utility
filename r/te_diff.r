library(ggplot2)

ca = commandArgs(trailing=T)
df.file = ca[1]
output.pdf = ca[2]

df = read.table(df.file, header=T, quote="\"")

xmin = quantile(df$fold, .007, na.rm=T)
xmax = quantile(df$fold, .993, na.rm=T)

ggplot(df, aes(x=fold, color=class)) +
    stat_ecdf(size=1.5, alpha=0.8) +
    scale_x_continuous("log2 RIP/Input", limits=c(xmin,xmax)) +
    scale_y_continuous("") +
    scale_color_manual("", values=c("#F46D43", "#66BD63")) +
    theme_bw() +
    theme(text=element_text(size=20)) +
    theme(legend.justification=c(1,0), legend.position=c(1,0))

#     scale_x_continuous("Differential expression test statistic", limits=c(xmin,xmax)) +
#     scale_color_brewer("", palette="Set1") +

ggsave(output.pdf)
