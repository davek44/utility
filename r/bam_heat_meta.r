library(ggplot2)

ca = commandArgs(trailing=T)
df.file = ca[1]
output.pdf = ca[2]

df = read.table(df.file, header=T, quote="\"")

if (ncol(df) == 2) {
    gp = ggplot(df, aes(x=Index, y=Coverage))
} else {
    gp = ggplot(df, aes(x=Index, y=Coverage, color=Type)) +
        scale_color_brewer(palette="Set1")
}

gp +
    geom_point() +
    geom_smooth() +
    theme_bw() +
    theme(text=element_text(size=16)) +
    theme(legend.justification=c(1,0), legend.position=c(1,0))

ggsave(output.pdf)
