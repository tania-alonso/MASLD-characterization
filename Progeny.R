library(readxl)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(DESeq2)
library(progeny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(limma)

# 1. Càrrega de dades
counts_df <- read_excel("C:/Users/Usuario/Documents/TÀNIA/UOC Master/TFM/Obtenció DEGs/test-counts.xlsx")

# 2. Assignació de Symbols (Traducció d'Ensembl a Gene Symbol)
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = counts_df$gene_id, 
                       column = "SYMBOL", 
                       keytype = "ENSEMBL", 
                       multiVals = "first")
counts_df$symbol <- gene_symbols

# 3. Eliminem NAs i duplicats de symbols (imprescindible per PROGENy)
counts_df_clean <- counts_df %>% 
  filter(!is.na(symbol)) %>% 
  distinct(symbol, .keep_all = TRUE)

# 4. Preparació de la matriu i metadades
sample_cols <- c("CTL_1","CTL_2","NASH_1","NASH_2",
                 "P548_F0","P554_F0","P599_FO","P511_HY","P579_HY","P585_HY")

clean_counts <- as.matrix(counts_df_clean[, sample_cols])
rownames(clean_counts) <- counts_df_clean$symbol

colData <- data.frame(
  sample = sample_cols,
  condition = factor(c("CTL", "CTL", "NASH", "NASH", "F0", "F0", "F0", "HY", "HY", "HY"))
)
rownames(colData) <- colData$sample

# 5. Normalització amb DESeq2 (VST)
dds <- DESeqDataSetFromMatrix(countData = clean_counts, 
                              colData = colData, 
                              design = ~ condition)
vsd <- vst(dds, blind = FALSE)
expression_matrix <- assay(vsd) 

# 6. Execució de PROGENy
pathway_activity <- progeny(expression_matrix, 
                            scale = TRUE, 
                            organism = "Human", 
                            top = 500)

# Anàlisi diferencial d'activitats ---

# 1. Creem una matriu d'activitats transposada (gens a les files)
act <- t(pathway_activity)

# 2. Definim el disseny experimental per limma
design <- model.matrix(~ condition, data = colData)
design
fit <- lmFit(act, design)

contr.matrix <- makeContrasts(
  NASHvsCTL = conditionF0-conditionHY,
  levels = design
)

cfit <- contrasts.fit(fit, contr.matrix)
efit <- eBayes(cfit)

pathway_results <- topTable(efit, number = Inf)

pathway_results$pathway <- rownames(pathway_results)
pathway_results$neg_log10_padj <- -log10(pathway_results$adj.P.Val)

pathway_results <- pathway_results %>%
  filter(!is.na(logFC)) %>%
  mutate(pathway = reorder(pathway, logFC))

ggplot(pathway_results, aes(x = logFC, y = pathway)) +
  geom_segment(aes(x = 0, xend = logFC, yend = pathway), color = "grey85") +
  geom_point(aes(size = neg_log10_padj, color = logFC)) +
  scale_color_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0) +
  theme_minimal() +
  labs(x = "PROGENy activity (logFC)", title = "F0 vs HY") 



library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Transformem i calculem mitjanes assegurant-nos de treure els NA
df_plot <- as.data.frame(pathway_activity) %>%
  mutate(condition = colData$condition) %>%
  pivot_longer(cols = -condition, names_to = "pathway", values_to = "score") %>%
  group_by(condition, pathway) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = 'drop') %>%
  filter(!is.na(mean_score)) #  eliminem files sense valor

# 2. Comprovació
print(head(df_plot)) 
# Hauries de veure una taula amb números a la columna mean_score.
# Si la taula surt buida (0 files), el problema ve de 'pathway_activity'.

# 3. El Gràfic
ggplot(df_plot, aes(x = mean_score, y = pathway, color = mean_score)) +
  geom_point(aes(size = abs(mean_score))) +
  scale_color_gradient2(low = "blue", mid = "grey", high = "red") +
  facet_wrap(~condition) +
  theme_minimal()
