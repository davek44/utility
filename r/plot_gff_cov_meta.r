library(ggplot2)
library(plyr)

ca = commandArgs(trailing=T)
df.file = ca[1]
out.pre = ca[2]
smooth.span = as.numeric(ca[3])
label.primary = ca[4]
label.control = ca[5]

df = read.table(df.file, header=T, quote="\"")

# unnormalized
if (ncol(df) == 2) {
    gp = ggplot(df, aes(x=Index, y=Coverage))
} else {
    gp = ggplot(df, aes(x=Index, y=Coverage, color=Type)) +
        scale_x_continuous("% in transcript") +
        scale_color_manual("", values=c("#F46D43", "#66BD63"), breaks=c("Primary","Control"), labels=c(label.primary, label.control))

        #scale_color_brewer(palette="Set1")
}

gp +
    geom_point() +
    stat_smooth(method="loess", span=smooth.span) +
    theme_bw() +
    theme(text=element_text(size=25)) +
    theme(legend.justification=c(1,0), legend.position=c(1,0))

ggsave(paste(out.pre,"_raw.pdf",sep=""))


# normalized
if (ncol(df) > 2) {
    # the values are so low, I want to boost them up
    fudge=10

    control.sum = sum(df[df$Type=="Control",]$Coverage)
    primary.sum = sum(df[df$Type!="Control",]$Coverage)

    df$Coverage.Norm = df$Coverage
    for (i in 1:nrow(df)) {
        if (df[i,"Type"] == "Control") {
            df[i,"Coverage.Norm"] = fudge * df[i,"Coverage"] / control.sum
        } else {
            df[i,"Coverage.Norm"] = fudge * df[i,"Coverage"] / primary.sum
        }
    }

    ggplot(df, aes(x=Index, y=Coverage.Norm, color=Type)) +
        scale_x_continuous("% in transcript") +
        scale_color_manual("", values=c("#F46D43", "#66BD63"), breaks=c("Primary","Control"), labels=c(label.primary, label.control)) +
        geom_point() +
        stat_smooth(method="loess", span=smooth.span) +
        scale_y_continuous("Normalized coverage") +
        theme_bw() +
        theme(text=element_text(size=25)) +
        theme(legend.justification=c(1,0), legend.position=c(1,0))

# scale_color_brewer(palette="Set1")

    ggsave(paste(out.pre,"_norm.pdf",sep=""))
}
