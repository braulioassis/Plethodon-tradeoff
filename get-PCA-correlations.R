##########################################
# Correlations for WGCNA Modules and PCA #
##########################################

# --- Load data ---
df <- read.csv("gene_counts_model_eigenvalues.csv")         # Gene expression eigenvalues per sample
delta <- read.csv("delta_ri_and_vo2_with_pc.csv")           # Phenotypic/physiological data with PCA values
df <- merge(df, delta, by.x = "ID", by.y = "id", all.x = TRUE)  # Merge by sample ID

# --- Load WGCNA module information ---
wgcna <- read.csv("WGCNA.csv")
wgcna <- wgcna[, c(2, 3)]  # Keep only columns 2 (SeqName) and 3 (ModuleColor)

# --- Define which variables correspond to WGCNA modules in df ---
vars <- colnames(df[10:36])  # Columns containing module eigengenes

# -------------------------------
# Correlation: PC1 vs Module Eigengenes
# -------------------------------

# Compute correlation between PC1 and each module eigengene
pc1.correlations <- sapply(vars, function(x) cor(df$pc1, df[[x]]))

# Compute associated p-values
getP <- sapply(vars, function(x) cor.test(df$pc1, df[[x]]))
getP <- t(getP)  # Transpose result for easier handling

# Combine results into a data frame
pc1.cor <- as.data.frame(pc1.correlations)
pc1.cor$Module <- row.names(pc1.cor)
getP <- as.data.frame(getP)
getP$Module <- row.names(getP)

# Extract p-values from cor.test() results
for (i in pc1.cor$Module) {
  pc1.cor$P_pc1[pc1.cor$Module == i] <- getP$p.value[getP$Module == i]  
}

# -------------------------------
# Correlation: PC2 vs Module Eigengenes
# -------------------------------

pc2.correlations <- sapply(vars, function(x) cor(df$pc2, df[[x]]))
getP <- sapply(vars, function(x) cor.test(df$pc2, df[[x]]))
getP <- t(getP)
pc2.cor <- as.data.frame(pc2.correlations)
pc2.cor$Module <- row.names(pc2.cor)
getP <- as.data.frame(getP)
getP$Module <- row.names(getP)

# Extract p-values for PC2 correlations
for (i in pc2.cor$Module) {
  pc2.cor$P_pc2[pc2.cor$Module == i] <- getP$p.value[getP$Module == i]  
}

# -------------------------------
# Combine and clean correlation results
# -------------------------------

# Merge PC1 and PC2 correlations side-by-side
pca.cor <- cbind(pc1.cor, pc2.cor)
pca.cor <- pca.cor[, c(2, 1, 3, 4, 6)]  # Reorder columns for readability
pca.cor$P_pc1 <- as.numeric(pca.cor$P_pc1)
pca.cor$P_pc2 <- as.numeric(pca.cor$P_pc2)

# Save WGCNA module vs PCA correlation results
write.csv(pca.cor, "WGCNA module correlations with PCA.csv", row.names = FALSE, quote = FALSE)


################################
# Single Gene PCA Correlations #
################################

# --- Reload data for gene-level analysis ---
df <- read.csv("gene_counts_model_eigenvalues.csv", header = TRUE)
delta <- read.csv("delta_ri_and_vo2_with_pc.csv")
df <- merge(df, delta, by.x = "ID", by.y = "id", all.x = TRUE)

wgcna <- read.csv("WGCNA.csv", header = TRUE)
wgcna <- wgcna[, c(2, 3)]  # Keep transcript and module color

# Define which columns correspond to individual gene expression values
vars <- colnames(df[37:9382])

# -------------------------------
# Correlation: PC1 vs Individual Genes
# -------------------------------

# Compute correlation for each gene
pc1.correlations <- sapply(vars, function(x) cor(df$pc1, df[[x]]))
pc1.cor <- as.data.frame(pc1.correlations)
pc1.cor$Transcript <- row.names(pc1.cor)

# Add module annotation from WGCNA
pc1.cor <- merge(pc1.cor, wgcna, by.x = "Transcript", by.y = "SeqName", all.x = TRUE)

# Compute absolute correlation values (for magnitude-based filtering)
pc1.cor$abs.pc1.correlations <- abs(pc1.cor$pc1.correlations)
pc1.cor <- pc1.cor[, c(1, 3, 2, 4)]  # Reorder for readability

# Identify top 1% most strongly correlated transcripts
pc1topone <- pc1.cor[pc1.cor$abs.pc1.correlations >= quantile(pc1.cor$abs.pc1.correlations, 0.99), ]

# Mark top transcripts in full correlation table
pc1.cor$pc1topone <- ifelse(pc1.cor$Transcript %in% pc1topone$Transcript, 1, 0)

# Sort by correlation strength (descending)
pc1.cor <- pc1.cor[order(pc1.cor$abs.pc1.correlations, decreasing = TRUE), ]

# -------------------------------
# Save gene-level correlation results
# -------------------------------
write.csv(pc1.cor, "Transcript correlations with PC1.csv", row.names = FALSE)
