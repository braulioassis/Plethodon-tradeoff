##########################################
# Correlations for WGCNA modules and PCA #
##########################################

df <- read.csv("gene_counts_model_eigenvalues.csv")
delta <- read.csv("delta_ri_and_vo2_with_pc.csv")
df <- merge(df, delta, by.x = "ID", by.y = "id", all.x = TRUE)

wgcna <- read.csv("WGCNA.csv")
wgcna <- wgcna[, c(2,3)]
vars <- colnames(df[10:36])

pc1.correlations <- sapply(vars, function(x) cor(df$pc1, df[[x]]))
getP <- sapply(vars, function(x) cor.test(df$pc1, df[[x]]))
getP <- t(getP)
pc1.cor <- as.data.frame(pc1.correlations)
pc1.cor$Module <- row.names(pc1.cor)
getP <- as.data.frame(getP)
getP$Module <- row.names(getP)

for (i in pc1.cor$Module) {
  pc1.cor$P_pc1[pc1.cor$Module == i] <- getP$p.value[getP$Module == i]  
}

pc2.correlations <- sapply(vars, function(x) cor(df$pc2, df[[x]]))
getP <- sapply(vars, function(x) cor.test(df$pc2, df[[x]]))
getP <- t(getP)
pc2.cor <- as.data.frame(pc2.correlations)
pc2.cor$Module <- row.names(pc2.cor)
getP <- as.data.frame(getP)
getP$Module <- row.names(getP)

for (i in pc2.cor$Module) {
  pc2.cor$P_pc2[pc2.cor$Module == i] <- getP$p.value[getP$Module == i]  
}

pca.cor <- cbind(pc1.cor, pc2.cor)
pca.cor <- pca.cor[, c(2, 1, 3, 4, 6)]
pca.cor$P_pc1 <- as.numeric(pca.cor$P_pc1)
pca.cor$P_pc2 <- as.numeric(pca.cor$P_pc2)

write.csv(pca.cor, "WGCNA module correlations with PCA.csv", row.names = F, quote = F)

################################
# Single gene PCA correlations #
################################

df <- read.csv("gene_counts_model_eigenvalues.csv", header = T)
delta <- read.csv("delta_ri_and_vo2_with_pc.csv")
df <- merge(df, delta, by.x = "ID", by.y = "id", all.x = TRUE)

wgcna <- read.csv("WGCNA.csv", header = T)
wgcna <- wgcna[, c(2,3)]
vars <- colnames(df[37:9382])

pc1.correlations <- sapply(vars, function(x) cor(df$pc1, df[[x]]))
pc1.cor <- as.data.frame(pc1.correlations)
pc1.cor$Transcript <- row.names(pc1.cor)
pc1.cor <- merge(pc1.cor, wgcna, by.x = "Transcript", by.y = "SeqName", all.x = TRUE)
pc1.cor$abs.pc1.correlations <- abs(pc1.cor$pc1.correlations)
pc1.cor <- pc1.cor[, c(1, 3, 2, 4)]

pc1topone <- pc1.cor[pc1.cor$abs.pc1.correlations >= quantile(pc1.cor$abs.pc1.correlations, 0.99), ]

pc1.cor$pc1topone <- ifelse(pc1.cor$Transcript %in% pc1topone$Transcript, 1, 0)

pc1.cor <- pc1.cor[order(pc1.cor$abs.pc1.correlations, decreasing = T), ]

write.csv(full, "Transcript correlations with PC1.csv", row.names = F)
