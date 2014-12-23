library(ggplot2)

ca = commandArgs(trailing=T)
df.file = ca[1]
output.pdf = ca[2]

df = read.table(df.file, header=T, quote="\"")

ggplot(df, aes(x=D1, y=D2, label=Label, color=Sample)) +
	geom_point(size=3, alpha=0.8) +
    theme_bw() +
    theme(text=element_text(size=22)) +
    coord_fixed()
    

ggsave(output.pdf)

# theme(legend.justification=c(1,0), legend.position=c(1,0))
# geom_text(size=5) +
