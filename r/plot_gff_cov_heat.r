library(ggplot2)

ca = commandArgs(trailing=T)
df.file = ca[1]
output.pdf = ca[2]
control = ca[3]

df = read.table(df.file, header=T, quote="\"")

gp = ggplot(df, aes(x=Index, y=Anchor, fill=Coverage)) +
    geom_tile()

if(control == "True") {
    gp = gp + scale_fill_gradient2(low="#377eb8", high="#e41a1c")
} else {
    gp = gp + scale_fill_gradient(low="white", high="#e41a1c")
}

gp +
    scale_y_discrete("") +
    theme_bw() +
    theme(text=element_text(size=16)) +
    theme(axis.ticks.y=element_blank(), axis.text.y=element_blank())

ggsave(output.pdf)
