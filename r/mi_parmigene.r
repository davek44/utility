library(parmigene)

ca = commandArgs(trailing=T)
df.file = ca[1]

df = read.table(df.file, header=T, quote="\"")

mi = knnmi(df$A, df$B, k=5)

cat(mi)
