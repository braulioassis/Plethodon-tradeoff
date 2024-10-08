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
ricor$abs.ri.correlations <- abs(ricor$ri.correlations)
ricor <- ricor[, c(1,3,2,4)]

ritopone <- ricor[ricor$abs.ri.correlations >= quantile(ricor$abs.ri.correlations, 0.99), ]
ribottomninety <- ricor[ricor$abs.ri.correlations <= quantile(ricor$abs.ri.correlations, 0.9), ]

ricor$ritopone <- ifelse(ricor$Transcript %in% ritopone$Transcript, 1, 0)
ricor$ribottomninety <- ifelse(ricor$Transcript %in% ribottomninety$Transcript, 1, 0)

vo2.correlations <- sapply(vars, function(x) cor(df$delta_vo2, df[[x]]))
vo2cor <- as.data.frame(vo2.correlations)
vo2cor$Transcript <- row.names(vo2cor)
vo2cor$abs.vo2.correlations <- abs(vo2cor$vo2.correlations)

vo2topone <- vo2cor[vo2cor$abs.vo2.correlations >= quantile(vo2cor$abs.vo2.correlations, 0.99), ]
vo2bottomninety <- vo2cor[vo2cor$abs.vo2.correlations <= quantile(vo2cor$abs.vo2.correlations, 0.9), ]

vo2cor$vo2topone <- ifelse(vo2cor$Transcript %in% vo2topone$Transcript, 1, 0)
vo2cor$vo2bottomninety <- ifelse(vo2cor$Transcript %in% vo2bottomninety$Transcript, 1, 0)

full <- merge(ricor,vo2cor, by = "Transcript")
full$decoupler <- "none"
full$decoupler[full$ritopone == 1 & full$vo2bottomninety == 1] <- "ri"
full$decoupler[full$vo2topone == 1 & full$ribottomninety == 1] <- "vo2"
full$decoupler[full$vo2topone == 1 & full$ritopone == 1] <- "driver"

write.csv(full, "Top 1percent of all transcript correlations.csv", row.names = F)

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

write.csv(full, "Top 10percent correlations for GO enrichment.csv", row.names = F)

#################################
# Get correlations and p-values #
#################################

# Delta ri and delta vo2

setwd("C:/Users/Braulio/Desktop/Tradeoff/Data")

df <- read.csv("gene_counts_model_eigenvalues.csv", header = T)
wgcna <- read.csv("WGCNA.csv")
wgcna$geneSymbol <- gsub(",","",wgcna$geneSymbol)
vars <- colnames(df[37:9382])

correlations <- sapply(vars, function(x) cor(df$delta_ri, df[[x]]))
ricor <- as.data.frame(correlations)
ricor$Transcript <- row.names(ricor)
ricor$correlations <- abs(ricor$correlations)

correlations <- sapply(vars, function(x) cor(df$delta_vo2, df[[x]]))
vo2cor <- as.data.frame(correlations)
vo2cor$Transcript <- row.names(vo2cor)
vo2cor$correlations <- abs(vo2cor$correlations)

for (i in ricor$Transcript) {
  ricor$Gene[ricor$Transcript == i] <- wgcna$geneSymbol[wgcna$gene_id == i]  
}

for (i in vo2cor$Transcript) {
  vo2cor$Gene[vo2cor$Transcript == i] <- wgcna$geneSymbol[wgcna$gene_id == i]  
}

ct <- sapply(vars, function(x) cor.test(df$delta_ri, df[[x]]))
ct <- t(ct)
ct <- as.data.frame(ct)
ct$Transcript <- row.names(ct)

for (i in ricor$Transcript) {
  ricor$P[ricor$Transcript == i] <- ct$p.value[ct$Transcript == i]  
}

ct <- sapply(vars, function(x) cor.test(df$delta_vo2, df[[x]]))
ct <- t(ct)
ct <- as.data.frame(ct)
ct$Transcript <- row.names(ct)

for (i in vo2cor$Transcript) {
  vo2cor$P[vo2cor$Transcript == i] <- ct$p.value[ct$Transcript == i]  
}

rifirstpercent <- ricor[ricor$correlations >= quantile(ricor$correlations, 0.99), ]
rilowpercent <- ricor[ricor$correlations <= quantile(ricor$correlations, 0.9), ]
vo2firstpercent <- vo2cor[vo2cor$correlations >= quantile(vo2cor$correlations, 0.99), ]
vo2lowpercent <- vo2cor[vo2cor$correlations <= quantile(vo2cor$correlations, 0.9), ]

ridecouplers <- rifirstpercent[rifirstpercent$Transcript %in% vo2lowpercent$Transcript,]
ridecouplers$P <- as.numeric(ridecouplers$P)
ridecouplers$abs.Pearson <- ridecouplers$correlations
ridecouplers <- ridecouplers[, c(2, 3, 5, 4)]
ridecouplers <- ridecouplers[order(ridecouplers$P), ]
write.csv(ridecouplers, "Delta ri top 1 percent and delta VO2 bottom 90 percent.csv",
          row.names = F, quote = F)

vo2decouplers <- vo2firstpercent[vo2firstpercent$Transcript %in% rilowpercent$Transcript,]
vo2decouplers$P <- as.numeric(vo2decouplers$P)
vo2decouplers$abs.Pearson <- vo2decouplers$correlations
vo2decouplers <- vo2decouplers[, c(2, 3, 5, 4)]
vo2decouplers <- vo2decouplers[order(vo2decouplers$P), ]
write.csv(vo2decouplers, "Delta vo2 top 1 percent and delta ri bottom 90 percent.csv",
          row.names = F, quote = F)

# Drivers

setwd("C:/Users/Braulio/Desktop/Tradeoff/Data")

df <- read.csv("gene_counts_model_eigenvalues.csv", header = T)
wgcna <- read.csv("WGCNA.csv")
wgcna$geneSymbol <- gsub(",","",wgcna$geneSymbol)
vars <- colnames(df[37:9382])

ri.correlations <- sapply(vars, function(x) cor(df$delta_ri, df[[x]]))
ricor <- as.data.frame(ri.correlations)
ricor$Transcript <- row.names(ricor)
ricor$ri.abs.correlations <- abs(ricor$ri.correlations)

vo2.correlations <- sapply(vars, function(x) cor(df$delta_vo2, df[[x]]))
vo2cor <- as.data.frame(vo2.correlations)
vo2cor$Transcript <- row.names(vo2cor)
vo2cor$vo2.abs.correlations <- abs(vo2cor$vo2.correlations)

for (i in ricor$Transcript) {
  ricor$Symbol[ricor$Transcript == i] <- wgcna$geneSymbol[wgcna$gene_id == i]  
}

for (i in vo2cor$Transcript) {
  vo2cor$Symbol[vo2cor$Transcript == i] <- wgcna$geneSymbol[wgcna$gene_id == i]  
}

ct <- sapply(vars, function(x) cor.test(df$delta_ri, df[[x]]))
ct <- t(ct)
ct <- as.data.frame(ct)
ct$Transcript <- row.names(ct)

for (i in ricor$Transcript) {
  ricor$P_ri[ricor$Transcript == i] <- ct$p.value[ct$Transcript == i]  
}

ct <- sapply(vars, function(x) cor.test(df$delta_vo2, df[[x]]))
ct <- t(ct)
ct <- as.data.frame(ct)
ct$Transcript <- row.names(ct)

for (i in vo2cor$Transcript) {
  vo2cor$P_vo2[vo2cor$Transcript == i] <- ct$p.value[ct$Transcript == i]  
}

rifirstpercent <- ricor[ricor$ri.abs.correlations >= quantile(ricor$ri.abs.correlations, 0.99), ]
vo2firstpercent <- vo2cor[vo2cor$vo2.abs.correlations >= quantile(vo2cor$vo2.abs.correlations, 0.99), ]

driversri <- rifirstpercent[rifirstpercent$Transcript %in% vo2firstpercent$Transcript,]
driversvo2 <- vo2firstpercent[vo2firstpercent$Transcript %in% rifirstpercent$Transcript,]
drivers <- cbind(driversri, driversvo2)

drivers$Gene <- driversri$Symbol
drivers <- drivers[, c(2, 11, 1, 3, 5, 6, 8, 10)]
drivers$Corr.sum <- drivers$ri.abs.correlations + drivers$vo2.abs.correlations
drivers <- drivers[order(drivers$Corr.sum, decreasing = T), ]
drivers$P_vo2 <- as.numeric(drivers$P_vo2)
drivers$P_ri <- as.numeric(drivers$P_ri)

write.csv(drivers, "Drivers top 1 percent delta ri and delta VO2.csv",
          row.names = F, quote = F)
