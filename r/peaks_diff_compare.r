library(ggplot2)

ca = commandArgs(trailing=T)
df.file = ca[1]
out.pre = ca[2]
rbp = ca[3]
tstat = ca[4]

df = read.table(df.file, header=T, quote="\"")

x.min = quantile(df$RIP, .005, na.rm=T)
x.max = quantile(df$RIP, .995, na.rm=T)

if(tstat == "True") {
    x.lab = paste(rbp, "fRIP/input diff stat")
} else {
    x.lab = paste(rbp, "log2 fRIP/input")
}

ggplot(df, aes(x=RIP, color=CLIP)) +
    geom_line(stat="density", size=1.5, alpha=0.8) +
    scale_x_continuous(x.lab, limits=c(x.min,x.max)) +
    scale_color_brewer(palette="Set1") +
    theme_bw() +
    theme(text=element_text(size=25)) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))

out.pdf = paste(out.pre, "_dens.pdf", sep="")
ggsave(out.pdf)


ggplot(df, aes(x=RIP, color=CLIP)) +
    stat_ecdf(size=1.5, alpha=0.8) +
    scale_x_continuous(x.lab, limits=c(x.min,x.max)) +
    scale_y_continuous("") +
    scale_color_brewer(palette="Set1") +
    theme_bw() +
    theme(text=element_text(size=25)) +
    theme(legend.justification=c(1,0), legend.position=c(1,0))

out.pdf = paste(out.pre, "_cdf.pdf", sep="")
ggsave(out.pdf)
