library(ggplot2)

ca = commandArgs(trailing=T)
df.file = ca[1]
output.pdf = ca[2]

df = read.table(df.file, header=T, quote="\"")

ggplot(df, aes(x=peak_i, y=cov)) +
    geom_point() +
    scale_x_continuous("Peak index") +
    scale_y_continuous("Coverage") +
    theme_bw() +
    theme(text=element_text(size=20))

ggsave(output.pdf)
