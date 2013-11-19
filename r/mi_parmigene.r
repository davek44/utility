library(parmigene)

ca = commandArgs(trailing=T)
df.file = ca[1]

df = read.table(df.file, header=T, quote="\"")

mi = knmi(df$A, df$B)

cat(mi)
