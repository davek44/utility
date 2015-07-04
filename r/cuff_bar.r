library(ggplot2)
library(RColorBrewer)

ca = commandArgs(trailing=T)
df.file = ca[1]
y.min = as.numeric(ca[2])
y.max = ca[3]
output.pdf = ca[4]

df = read.table(df.file, header=T, quote="\"")

if (y.max == "None") {
	y.max = max(df$conf_hi)
} else {
	y.max = as.numeric(y.max)
}

color.count = length(unique(df$Sample))
get.pal = colorRampPalette(brewer.pal(11, "Spectral"))

ggplot(df, aes(x=Sample, y=FPKM, ymin=conf_lo, ymax=conf_hi, fill=Sample)) +
	geom_bar(stat="identity", fill=get.pal(color.count)) +
	geom_errorbar(width=0.5) +
	scale_x_discrete("", limits=df$Sample) +
	guides(fill=FALSE) +
	theme_bw() +
	theme(text=element_text(size=22)) +
	theme(axis.text.x=element_text(angle=315, hjust=0, vjust=1))

# 	scale_fill_brewer(palette="Spectral") +

ggsave(output.pdf)
