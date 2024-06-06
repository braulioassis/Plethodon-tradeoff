library(goseq)
library(plyr)
setwd("C:/Users/Braulio/Desktop/Tradeoff/Data/")

df <- read.csv("Top 10percent decouplers for GO enrichment.csv")
skin_all_genes <- read.csv("WGCNA.csv")

gene_lengths <- read.csv("gene_lengths.csv")
colnames(gene_lengths) <- c("gene_id","length")
skin_all_genes <- merge(skin_all_genes,gene_lengths, by.y = "gene_id")
length_bias <- as.vector(skin_all_genes[,length(skin_all_genes)])

skin_categories <- read.csv("new_categories.csv")
skin_categories <- skin_categories[skin_categories[['id']] %in% skin_all_genes$gene_id, ,drop = FALSE]
skin_categories <- unique(skin_categories[ , 1:2 ] )

skin_all_genes <- as.vector(skin_all_genes[,1])
gene_vector <- df$Transcript[df$decoupler == "ri"]
gene_vector <- as.integer(skin_all_genes%in%gene_vector)
names(gene_vector) <- skin_all_genes

colnames(skin_categories) <- c("gene_id","go_term")
go_counts <- ddply(skin_categories, c("go_term"), summarise,N = length(go_term))
go_counts <- go_counts[-which(go_counts$N < 9),]
go_counts <- go_counts[-which(go_counts$N > 500),]
skin_categories <- skin_categories[skin_categories[['go_term']] %in% go_counts$go_term, ,drop = FALSE]

pwf <- nullp(gene_vector,bias.data = length_bias)
rownames(pwf) <- skin_all_genes
go <- goseq(pwf,gene2cat=skin_categories,test.cats=c("GO:MF","GO:CC","GO:BP"),use_genes_without_cat=TRUE)

write.csv(go,"! GO RI decouplers top 10 percent 429 of only known genes.csv")