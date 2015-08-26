library(ggplot2)
library(plyr)

ca = commandArgs(trailing=T)
df.file = ca[1]
out.pdf = ca[2]
scale = as.numeric(ca[3])

df = read.table(df.file, header=T)

ggplot(df, aes(x=TEs, y=stat_mid, ymin=stat_low, ymax=stat_hi)) +
    geom_pointrange() +
    stat_smooth(se=FALSE, color="black", lty=2) +
    scale_y_continuous("log2 fRIP/input") +
    theme_bw() +
    theme(text=element_text(size=(sqrt(scale)*28)))

ggsave(out.pdf, scale=scale)
