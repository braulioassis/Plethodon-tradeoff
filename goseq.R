##############################################
# GO Enrichment Analysis for Top Correlated Genes
#   This script performs GO enrichment analysis using the goseq package.
#   It adjusts for gene length bias, filters categories by size,
#   and tests enrichment for the top correlated genes in WGCNA modules.
##############################################

library(goseq)
library(plyr)

#-----------------------------
# 1. Load Input Data
#-----------------------------

# Top correlated genes (e.g., top 10% based on correlation with traits)
df <- read.csv("Top 10percent correlations for GO enrichment.csv")

# WGCNA module assignments for all genes
skin_all_genes <- read.csv("WGCNA.csv")

# Gene lengths for bias correction
gene_lengths <- read.csv("gene_lengths.csv")
colnames(gene_lengths) <- c("gene_id", "length")

# Merge WGCNA gene list with gene lengths
skin_all_genes <- merge(skin_all_genes, gene_lengths, by.y = "gene_id")

# Extract gene length vector for goseq bias correction
length_bias <- as.vector(skin_all_genes[, length(skin_all_genes)])

#-----------------------------
# 2. Prepare GO Category Data
#-----------------------------

# Load gene-to-GO term mapping
skin_categories <- read.csv("new_categories.csv")

# Keep only GO terms for genes present in the dataset
skin_categories <- skin_categories[skin_categories[['id']] %in% skin_all_genes$gene_id, , drop = FALSE]

# Keep only unique gene–GO associations
skin_categories <- unique(skin_categories[, 1:2])

# Convert gene list to a vector
skin_all_genes <- as.vector(skin_all_genes[, 1])

#-----------------------------
# 3. Define Gene Set of Interest
#-----------------------------

# Choose which correlation set to analyze:
#   Options: "driver", "vo2", or "ri"
#   Change the string below to select the target trait.
gene_vector <- df$Transcript[df$decoupler == "driver"]

# Create binary vector (1 = target gene, 0 = background)
gene_vector <- as.integer(skin_all_genes %in% gene_vector)
names(gene_vector) <- skin_all_genes

#-----------------------------
# 4. Filter GO Categories by Size
#-----------------------------

# Rename columns for goseq compatibility
colnames(skin_categories) <- c("gene_id", "go_term")

# Count how many genes belong to each GO term
go_counts <- ddply(skin_categories, c("go_term"), summarise, N = length(go_term))

# Remove overly small or large categories to avoid bias
go_counts <- go_counts[-which(go_counts$N < 9), ]     # Remove terms with <9 genes
go_counts <- go_counts[-which(go_counts$N > 500), ]   # Remove terms with >500 genes

# Keep only GO terms within the filtered range
skin_categories <- skin_categories[skin_categories[['go_term']] %in% go_counts$go_term, , drop = FALSE]

#-----------------------------
# 5. Run GOseq Enrichment Analysis
#-----------------------------

# Estimate probability weighting function (PWF) to correct for gene length bias
pwf <- nullp(gene_vector, bias.data = length_bias)
rownames(pwf) <- skin_all_genes

# Run GO enrichment for Biological Process (BP), Molecular Function (MF), and Cellular Component (CC)
go <- goseq(
  pwf,
  gene2cat = skin_categories,
  test.cats = c("GO:BP", "GO:MF", "GO:CC"),
  use_genes_without_cat = TRUE
)

#-----------------------------
# 6. Save Output
#-----------------------------

# Export GO enrichment results to CSV
write.csv(go, "GO drivers top 10 percent 267 of 6953 genes.csv")
