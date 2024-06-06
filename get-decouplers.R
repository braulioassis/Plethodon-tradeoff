setwd("C:/Users/Braulio/Desktop/Tradeoff/Data")

# Analysis keeping in uncharacterized transcripts for single gene evaluation

df <- read.csv("gene_counts_model_eigenvalues.csv", header = T)
wgcna <- read.csv("WGCNA.csv", header = T)
wgcna <- wgcna[, c(2,3)]
vars <- colnames(df[37:9382])

ri.correlations <- sapply(vars, function(x) cor(df$delta_ri, df[[x]]))
ricor <- as.data.frame(ri.correlations)
ricor$Transcript <- row.names(ricor)
ricor <- merge(ricor, wgcna, by.x = "Transcript", by.y = "SeqName", all.x = TRUE)
ricor$ri.correlations <- abs(ricor$ri.correlations)
ricor <- ricor[, c(1,3,2)]

ritopone <- ricor[ricor$ri.correlations >= quantile(ricor$ri.correlations, 0.99), ]
ribottomninety <- ricor[ricor$ri.correlations <= quantile(ricor$ri.correlations, 0.9), ]

ricor$ritopone <- ifelse(ricor$Transcript %in% ritopone$Transcript, 1, 0)
ricor$ribottomninety <- ifelse(ricor$Transcript %in% ribottomninety$Transcript, 1, 0)

vo2.correlations <- sapply(vars, function(x) cor(df$delta_vo2, df[[x]]))
vo2cor <- as.data.frame(vo2.correlations)
vo2cor$Transcript <- row.names(vo2cor)
vo2cor$vo2.correlations <- abs(vo2cor$vo2.correlations)

vo2topone <- vo2cor[vo2cor$vo2.correlations >= quantile(vo2cor$vo2.correlations, 0.99), ]
vo2bottomninety <- vo2cor[vo2cor$vo2.correlations <= quantile(vo2cor$vo2.correlations, 0.9), ]

vo2cor$vo2topone <- ifelse(vo2cor$Transcript %in% vo2topone$Transcript, 1, 0)
vo2cor$vo2bottomninety <- ifelse(vo2cor$Transcript %in% vo2bottomninety$Transcript, 1, 0)

full <- merge(ricor,vo2cor, by = "Transcript")
full$decoupler <- "none"
full$decoupler[full$ritopone == 1 & full$vo2bottomninety == 1] <- "ri"
full$decoupler[full$vo2topone == 1 & full$ribottomninety == 1] <- "vo2"
full$decoupler[full$vo2topone == 1 & full$ritopone == 1] <- "driver"

write.csv(full, "Top 1percent decouplers for all transcript correlations.csv", row.names = F)

# Analysis removing uncharacterized transcripts for GO enrichment

df <- read.csv("gene_counts_model_eigenvalues.csv", header = T)
wgcna <- read.csv("WGCNA.csv", header = T)
wgcna <- wgcna[, c(2,3)]
vars <- colnames(df[37:9382])
# To get only "known" genes
wgcna <- subset(wgcna, !grepl("-NA-", geneSymbol))
wgcna <- subset(wgcna, !grepl("hypoth", geneSymbol))
knowngenes <- wgcna$SeqName
vars <- vars[vars %in% knowngenes]

ri.correlations <- sapply(vars, function(x) cor(df$delta_ri, df[[x]]))
ricor <- as.data.frame(ri.correlations)
ricor$Transcript <- row.names(ricor)
ricor <- merge(ricor, wgcna, by.x = "Transcript", by.y = "SeqName", all.x = TRUE)
ricor$ri.correlations <- abs(ricor$ri.correlations)
ricor <- ricor[, c(1,3,2)]

ritopten <- ricor[ricor$ri.correlations >= quantile(ricor$ri.correlations, 0.9), ]
ribottomninety <- ricor[ricor$ri.correlations <= quantile(ricor$ri.correlations, 0.9), ]

ricor$ritopten <- ifelse(ricor$Transcript %in% ritopten$Transcript, 1, 0)
ricor$ribottomninety <- ifelse(ricor$Transcript %in% ribottomninety$Transcript, 1, 0)

vo2.correlations <- sapply(vars, function(x) cor(df$delta_vo2, df[[x]]))
vo2cor <- as.data.frame(vo2.correlations)
vo2cor$Transcript <- row.names(vo2cor)
vo2cor$vo2.correlations <- abs(vo2cor$vo2.correlations)

vo2topten <- vo2cor[vo2cor$vo2.correlations >= quantile(vo2cor$vo2.correlations, 0.9), ]
vo2bottomninety <- vo2cor[vo2cor$vo2.correlations <= quantile(vo2cor$vo2.correlations, 0.9), ]

vo2cor$vo2topten <- ifelse(vo2cor$Transcript %in% vo2topten$Transcript, 1, 0)
vo2cor$vo2bottomninety <- ifelse(vo2cor$Transcript %in% vo2bottomninety$Transcript, 1, 0)

full <- merge(ricor,vo2cor, by = "Transcript")
full$decoupler <- "none"
full$decoupler[full$ritopten == 1 & full$vo2bottomninety == 1] <- "ri"
full$decoupler[full$vo2topten == 1 & full$ribottomninety == 1] <- "vo2"
full$decoupler[full$vo2topten == 1 & full$ritopten == 1] <- "driver"

write.csv(full, "Top 10percent decouplers for GO enrichment.csv", row.names = F)
