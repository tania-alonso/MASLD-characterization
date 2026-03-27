getwd()
setwd("C:/Users/Usuario/Documents/TÀNIA/UOC Master/TFM/Obtenció DEGs")
list.files()

library(readxl)
library(readxl)
library(DESeq2)

# -------------------------------
# 1️⃣ Llegir dades
# -------------------------------
counts <- read_excel("test-counts.xlsx")
counts <- as.data.frame(counts)

# Posar gene_id com a rownames
rownames(counts) <- counts$gene_id

# Crear lookup correcte
gene_name_lookup <- counts$gene_name
names(gene_name_lookup) <- rownames(counts)

# Treure la columna gene_id
counts <- counts[ , !(names(counts) %in% c("gene_id"))]

# Assegurar que totes les columnes de comptatge siguin numèriques
numeric_cols <- setdiff(names(counts), "gene_name")
counts[ , numeric_cols] <- lapply(counts[ , numeric_cols], as.numeric)

# -------------------------------
# 2️⃣ Crear matriu només amb les mostres que necessitem
# -------------------------------
# , 
counts_small <- counts[, c( "CTL_1",	"CTL_2",	"NASH_1", "NASH_2")]

coldata_small <- data.frame(
  condition = c( 	 "CTL", "CTL" , "NASH", "NASH")
)
rownames(coldata_small) <- colnames(counts_small)


# IMPORTANT: fixar nivell de referència
coldata_small$condition <- factor(coldata_small$condition)
coldata_small$condition <- relevel(coldata_small$condition, ref = "CTL")
# DESeq2 calcula log2FoldChange = ALTRE GRUP / GRUP DE REF

# -------------------------------
# 3️⃣ Crear objecte DESeq2
# -------------------------------
dds_small <- DESeqDataSetFromMatrix(
  countData = counts_small,
  colData = coldata_small,
  design = ~ condition
)

# -------------------------------
# 4️⃣ Filtrar gens amb counts baixos (opcional)
# -------------------------------
# dds_small <- dds_small[rowSums(counts(dds_small)) > 10, ]

# -------------------------------
# 5️⃣ Correr DESeq2
# -------------------------------
dds_small <- DESeq(dds_small)

res_small <- results(dds_small, independentFiltering = FALSE, alpha=0.05)

# Eliminar NAs
res_small <- res_small[!is.na(res_small$pvalue) & !is.na(res_small$log2FoldChange), ]

# backgorund genes for enrich
background_genes <- rownames(res_small)
length(background_genes)

# Convertir Ensembl IDs a gene symbols
background_symbols <- gene_name_lookup[background_genes]

# Treure NAs
background_symbols <- background_symbols[!is.na(background_symbols)]

# Treure duplicats
background_symbols <- unique(background_symbols)

length(background_symbols)
       
write.table(background_symbols,
             "background_genes_symbols.txt",
             quote = FALSE,
             row.names = FALSE,
             col.names = FALSE)

# -------------------------------
# 6️⃣ Filtrar DEGs(p-value & fold-change)
# -------------------------------
# UP DEGS 

degs_up_small <- res_small[which(res_small$pvalue <= 0.05 & res_small$log2FoldChange >= 0.585), ]
nrow(degs_up_small)

degs_down_small <- res_small[which(res_small$pvalue <= 0.05 & res_small$log2FoldChange <= -0.585), ]


# Convertir a data frame
degs_up_small <- as.data.frame(degs_up_small)
degs_down_small <- as.data.frame(degs_down_small)

# Afegir columna gene_name
degs_up_small$gene_name <- gene_name_lookup[rownames(degs_up_small)]
degs_down_small$gene_name <- gene_name_lookup[rownames(degs_down_small)]

head(degs_up_small)
nrow(degs_up_small)

# -------------------------------
# 7️⃣ Exportar a CSV amb gene_name
# -------------------------------
 write.csv2(degs_up_small, "DEGs-UP_NASH_vs_CTL.csv", row.names = TRUE)
 write.csv2(degs_down_small, "DEGs-DOWN_NASH_vs_CTL.csv", row.names = TRUE)


# vOlcano plots
library(ggplot2)
library(ggrepel)

res_df <- as.data.frame(res_small)
res_df
nrow(res_df)

# eliminar NA si cal (ja está fet al res_small)
# crear categoria
res_df$group <- "NO"

res_df$group[res_df$pvalue <= 0.05 & res_df$log2FoldChange >= 1] <- "UP"
res_df$group[res_df$pvalue <= 0.05 & res_df$log2FoldChange <= -1] <- "DOWN"
res_df$gene_name <- gene_name_lookup[rownames(res_df)]

# comptar
counts_table <- table(res_df$group)

n_up <- counts_table["UP"]
n_down <- counts_table["DOWN"]
n_no <- counts_table["NO"]

library(dplyr)

# seleccionar top significatius
top_up <- res_df %>%
  filter(group == "UP") %>%
  arrange(pvalue) %>%
  slice_head(n = 10)

top_down <- res_df %>%
  filter(group == "DOWN") %>%
  arrange(pvalue) %>%
  slice_head(n = 10)

top_genes <- rbind(top_up, top_down)


# volcano plot amb números
ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = group)) +
  geom_point(alpha = 0.6, size = 1.5) +
  
  geom_text_repel(
    data = top_genes,
    aes(label = gene_name),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "black",
    max.overlaps = Inf
  ) +
  
  scale_color_manual(
    values = c("DOWN"="blue", "NO"="grey", "UP"="red"),
    labels = c(
      paste0("DOWN (", n_down, ")"),
      paste0("NO (", n_no, ")"),
      paste0("UP (", n_up, ")")
    )
  ) +
  
  geom_vline(xintercept = c(-1, 1), linetype="dashed") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  
  theme_minimal() +
  labs(title = "F0 vs NASH",
       x = "log2 Fold Change",
       y = "-log10(pvalue)",
       color = "")

