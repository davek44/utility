library(ggplot2)

ca = commandArgs(trailing=T)
df.file = ca[1]
output.pdf = ca[2]

df = read.table(df.file, header=T, quote="\"")

if (ncol(df) == 2) {
    gp = ggplot(df, aes(x=Index, y=Coverage))
} else {
    gp = ggplot(df, aes(x=Index, y=Coverage, colour=Type))
}

gp +
    geom_point() +
    geom_smooth() +
    theme_bw()

ggsave(output.pdf)
