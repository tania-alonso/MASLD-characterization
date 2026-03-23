getwd()

# PCA
# Paquets necessaris
library(readxl)
library(dplyr)

# llegir dades AT+EVs
data <- readxl::read_excel("ATEVsn_resultados_proteómica_final.xlsx", sheet = "Proteins")
mostres <- data[, grep("^Sample", colnames(data))]
healthy_samples <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6", "Sample7", "Sample8", "Sample10")
F0_samples <- c("Sample11", "Sample12", "Sample13", "Sample14", "Sample15", "Sample16", "Sample17",
                 "Sample18", "Sample19", "Sample22", "Sample23",  "Sample24", "Sample25")
F4_samples <- c("Sample26",  "Sample27", "Sample28", "Sample29", "Sample30", "Sample31")

 data <- as.data.frame(data)  
 data
 rownames(data) <- data$`Entry Name`
  data
# llegir dades AT
# data <- readxl::read_excel("AT_REANALYSED.xlsx", sheet = "SEL INFO")
# head(data)
# 
# data <- as.data.frame(data)  
# rownames(data) <- data$`Entry Name`
# 
# # 2. Noms de les columnes que contenen mostres
# sample_cols <- c("P393","P426","P497","P499","P538",
#                   "P376","P556 AT","P366","1468 AT","P360",
#                   "P495","P396","P409","P397","P401")
# 
# healthy_samples <- c("P393","P426","P497","P499","P538")
# F0_samples <- c("P376","P556 AT","P366","1468 AT","P360")
# F4_samples <- c("P495","P396","P409","P397","P401")

# 3. Converteix a numeric (substitueix comes per punts)
# data2 <- data %>%
#    mutate(across(all_of(sample_cols),
#                 ~as.numeric(gsub(",", ".", .))))
# 
# # 4. Extraiem només les columnes de mostres
# mostres <- data2[, sample_cols]

mostres
mostres <- as.data.frame(lapply(mostres, function(x) as.numeric(gsub(",", ".", x))), check.names = FALSE)
rownames(mostres) <- data$`Entry Name`
mostres_log <- log2(mostres + 1)
# Comptar valors no NA per proteïna
mostres_log[is.na(mostres_log)] <- 0
sum(is.na(mostres_log))
mostres
rownames(mostres_log) <- rownames(mostres)   # proteïnes
colnames(mostres_log) <- colnames(mostres)   # samples

mostres_log

# var_prot <- apply(mostres_log, 1, var)       # variància per proteïna
# mostres_top <- mostres_log[var_prot > quantile(var_prot, 0.75), ]  # top 25% més variable
# pca_top <- prcomp(t(mostres_top), center=TRUE, scale.=FALSE)

# 2. PCA sense escalar (scale=FALSE)
mostres_t <- t(as.matrix(mostres_log))
mostres_t
pca <- prcomp(mostres_t, center = TRUE, scale. = FALSE)
summary(pca)

# 2.2 Trobar el nombre òptim de clusters amb Elbow i Silhouette
library(factoextra)
fviz_nbclust(pca$x[,1:2], kmeans, method = "wss")       # Elbow
fviz_nbclust(pca$x[,1:2], kmeans, method = "silhouette") # Silhouette


# 3. Definim clusters manualment (pots fer clustering amb k-means si vols)
# Exemple: fem 3 clusters amb k-means
set.seed(123)
km <- kmeans(pca$x[,1:2], centers = 3)
clusters <- km$cluster

# 4. Plot PCA amb colors per cluster

# Assignem els grups correctament segons els noms de les columnes
grups <- rep(NA, ncol(mostres_log))

grups[colnames(mostres) %in% healthy_samples] <- "SA"
grups[colnames(mostres) %in% F0_samples] <- "F0"
grups[colnames(mostres) %in% F4_samples] <- "F4"
# plot per grups 
grups_factor <- factor(grups, levels = c("SA", "F0", "F4"))

colors <- c("red","green","blue")[grups_factor]
# Comprova
colnames(mostres)
colnames(mostres) %in% healthy_samples
table(grups)
grups

# 4. Plot PCA amb colors per cluster
plot(pca$x[,1], pca$x[,2], 
     col = clusters, pch = 19,
     xlab = paste0("PC1 (", round(summary(pca)$importance[2,1]*100,1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca)$importance[2,2]*100,1), "%)"),
     main = "PCA Proteòmica amb 3 clusters")
text(pca$x[,1], pca$x[,2], labels = colnames(mostres), pos = 3, cex = 0.7)

# PCA ja calculat (suposem que pca <- prcomp(t(mostres_log), ...))
plot(pca$x[,1], pca$x[,2],
     col = colors,
     pch = 19,
     xlab = "PC1", ylab = "PC2",
     main = "PCA Proteòmica amb grups biològics reals")
legend("topright", legend=c("SA","F0","F4"), col=c("red","green","blue"), pch=19)
text(pca$x[,1], pca$x[,2], labels = colnames(mostres), pos = 3, cex = 0.7)


# identificar quines proteines son les mes importants a pc1 i pc2
pca$rotation
loadings <- pca$rotation
rownames(loadings)
top_PC1 <- sort(abs(loadings[,1]), decreasing = TRUE)[1:20]
top_PC1
top_PC2 <- sort(abs(loadings[,2]), decreasing = TRUE)[1:20]
top_PC2
names(top_PC1)
names(top_PC2)

proteines_clau <- names(top_PC1)[1:20]
proteines_clau2 <- names(top_PC2)[1:20]

mostres_log[proteines_clau, ]

heatmap(as.matrix(mostres_log[proteines_clau, ]),
        scale = "row")

library(pheatmap)
# Data frame per pheatmap
ncol(mostres_log)
annotation_col <- data.frame(Group = factor(grups))
rownames(annotation_col) <- colnames(mostres_log)

# Heatmap amb anotacions
pheatmap(mostres_log[proteines_clau, ],
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         color = colorRampPalette(c("blue","yellow","red"))(50),
         main = "Top 20 proteïnes de PC1 amb grups biològics")


# Heatmap amb anotacions
pheatmap(mostres_log[proteines_clau2, ],
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         color = colorRampPalette(c("blue","yellow","red"))(50),
         main = "Top 20 proteïnes de PC2 amb grups biològics")
