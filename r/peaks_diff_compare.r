library(ggplot2)

ca = commandArgs(trailing=T)
df.file = ca[1]
out.pre = ca[2]

df = read.table(df.file, header=T, quote="\"")

x.min = quantile(df$Test_stat, .003, na.rm=T)
x.max = quantile(df$Test_stat, .997, na.rm=T)

ggplot(df, aes(x=Test_stat, colour=Peak)) +
    geom_line(stat="density") +
    scale_x_continuous(limits=c(x.min,x.max)) +
    scale_color_brewer(palette="Set1") +
    theme_bw() +
    theme(text=element_text(size=16)) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))

out.pdf = paste(out.pre, "_dens.pdf", sep="")
ggsave(out.pdf)


ggplot(df, aes(x=Test_stat, colour=Peak)) +
    stat_ecdf() +
    scale_x_continuous(limits=c(x.min,x.max)) +
    scale_y_continuous("") +
    scale_color_brewer(palette="Set1") +
    theme_bw() +
    theme(text=element_text(size=16)) +
    theme(legend.justification=c(1,0), legend.position=c(1,0))

out.pdf = paste(out.pre, "_cdf.pdf", sep="")
ggsave(out.pdf)
