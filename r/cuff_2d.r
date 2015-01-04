library(ggplot2)

ca = commandArgs(trailing=T)
df.file = ca[1]
output.pdf = ca[2]
square = ca[3]

df = read.table(df.file, header=T, quote="\"")

if (square == "True") {
	d1.span = max(df$D1) - min(df$D1)
	d2.span = max(df$D2) - min(df$D2)
	plot.ratio = d1.span/d2.span
} else {
	plot.ratio = 1
}


ggplot(df, aes(x=D1, y=D2, label=Label, color=Sample)) +
	geom_point(size=3, alpha=0.8) +
	scale_x_continuous("") +
	scale_y_continuous("") +
	scale_color_discrete("") +
	theme_bw() +
    theme(text=element_text(size=22)) +
    coord_fixed(ratio=plot.ratio) +
    theme(legend.justification=c(1,0), legend.position=c(1,0))
    
ggsave(output.pdf)

# theme(legend.justification=c(1,0), legend.position=c(1,0))
# geom_text(size=5) +
