library(cummeRbund)

cuff = readCufflinks()

############################################
# density plot
############################################
csDensity(genes(cuff), rep=T, pseudocount=0.1)
ggsave("density.pdf")

############################################
# dendrogram
############################################
pdf("dendro.pdf")
csDendro(genes(cuff), rep=T, pseudocount=1)
dev.off()

############################################
# MDS
############################################
MDSplot(genes(cuff), rep=T, pseudocount=1) +
	coord_fixed() +
	theme_bw() +
    theme(text=element_text(size=15))

ggsave("mds.pdf")

############################################
# PCA
############################################
PCAplot(genes(cuff), x="PC1", y="PC2", rep=T, pseudocount=1) +
	coord_fixed() +
	theme_bw() +
    theme(text=element_text(size=15))

ggsave("pca12.pdf")

PCAplot(genes(cuff), x="PC2", y="PC3", rep=T, pseudocount=1) +
	coord_fixed() +
	theme_bw() +
    theme(text=element_text(size=15))

ggsave("pca23.pdf")
