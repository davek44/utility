library(ggplot2)
library(plyr)

ca = commandArgs(trailing=T)
df.file = ca[1]
out.pre = ca[2]

df = read.table(df.file, header=T, quote="\"")

# unnormalized
if (ncol(df) == 2) {
    gp = ggplot(df, aes(x=Index, y=Coverage))
} else {
    gp = ggplot(df, aes(x=Index, y=Coverage, color=Type)) +
        scale_color_manual("", values=c("#F46D43", "#66BD63"))

        #scale_color_brewer(palette="Set1")
}

gp +
    geom_point() +
    geom_smooth() +
    theme_bw() +
    theme(text=element_text(size=20)) +
    theme(legend.justification=c(1,0), legend.position=c(1,0))

ggsave(paste(out.pre,"_raw.pdf",sep=""))


# normalized
if (ncol(df) == 2) {
    cov.sum = sum(df$Coverage)
    df$Coverage.Norm = df$Coverage / cov.sum

    gp = ggplot(df, aes(x=Index, y=Coverage.Norm))
} else {
    control.sum = sum(df[df$Type=="Control",]$Coverage)
    primary.sum = sum(df[df$Type!="Control",]$Coverage)

    df$Coverage.Norm = df$Coverage
    for (i in 1:nrow(df)) {
        if (df[i,"Type"] == "Control") {
            df[i,"Coverage.Norm"] = df[i,"Coverage"] / control.sum
        } else {
            df[i,"Coverage.Norm"] = df[i,"Coverage"] / primary.sum
        }
    }

    gp = ggplot(df, aes(x=Index, y=Coverage.Norm, color=Type)) +
        scale_color_manual("", values=c("#F46D43", "#66BD63"))
        #scale_color_brewer(palette="Set1")
}

gp +
    geom_point() +
    geom_smooth() +
    scale_y_continuous("Normalized coverage") +
    theme_bw() +
    theme(text=element_text(size=20)) +
    theme(legend.justification=c(1,0), legend.position=c(1,0))

ggsave(paste(out.pre,"_norm.pdf",sep=""))
