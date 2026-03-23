# ---------------------------------
# Anàlisi Proteòmica Completa
# ---------------------------------

# Llibreries
library(tidyverse)
library(readxl)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# -----------------------------
# 1. Llegir les 3 comparacions
# -----------------------------
hy_f0 <- read_excel("ATEVsn_resultados_proteómica_final.xlsx", sheet = "Healthy vs FO")
hy_f4 <- read_excel("ATEVsn_resultados_proteómica_final.xlsx", sheet = "Healthy vs F4")
f0_f4 <- read_excel("ATEVsn_resultados_proteómica_final.xlsx", sheet = "F0 vs F4")

# -----------------------------
# 2. Seleccionar columnes importants
# -----------------------------
hy_f0_sel <- hy_f0 %>%
  dplyr::select(
    Protein_ID = `T: Protein.Group`,
    Entry_Name = `T: Protein.Names`,
    Description = `T: First.Protein.Description`,
    Diff = `N: Student's T-test Difference Healthy_F0`,
    p_value = `p value`,
    starts_with("P")
  ) %>%
  dplyr::mutate(Diff = -Diff)


hy_f4_sel <- hy_f4 %>%
  dplyr::select(
    Protein_ID = `T: Protein.Group`,
    Entry_Name = `T: Protein.Names`,
    Description = `T: First.Protein.Description`,
    Diff = `N: Student's T-test Difference Healthy_F4`,
    p_value = `p value`,
    starts_with("P")
  ) %>%
  dplyr::mutate(Diff = -Diff)


f0_f4_sel <- f0_f4 %>%
  dplyr::select(
    Protein_ID = `T: Protein.Group`,
    Entry_Name = `T: Protein.Names`,
    Description = `T: First.Protein.Description`,
    Diff = `N: Student's T-test Difference F0_F4`,
    p_value = `p value`,
    starts_with("P")
  ) %>%
  dplyr::mutate(Diff = -Diff)


# -----------------------------
# 3. Proteïnes diferencials (filter: p<0.05 i |Diff|>0.585)
# -----------------------------
sig_hy_f0 <- hy_f0_sel %>%
  filter(!is.na(Diff), !is.na(p_value)) %>%
  filter(p_value < 0.05 & abs(Diff) > 0.585) %>%
  mutate(Regulation = ifelse(Diff > 0.585, "UP", "DOWN"))


sig_hy_f4 <- hy_f4_sel %>%
  filter(!is.na(Diff), !is.na(p_value)) %>%
  filter(p_value < 0.05 & abs(Diff) > 0.585) %>%
  mutate(Regulation = ifelse(Diff > 0.585, "UP", "DOWN"))

sig_f0_f4 <- f0_f4_sel %>%
  filter(!is.na(Diff), !is.na(p_value)) %>%
  filter(p_value < 0.05 & abs(Diff) > 0.585) %>%
  mutate(Regulation = ifelse(Diff > 0.585, "UP", "DOWN"))


data.frame(
  Comparison = c("HY vs F0-F1","HY vs F3-F4","F0-F1 vs F3-F4"),
  UP = c(sum(sig_hy_f0$Regulation=="UP"),
         sum(sig_hy_f4$Regulation=="UP"),
         sum(sig_f0_f4$Regulation=="UP")),
  DOWN = c(sum(sig_hy_f0$Regulation=="DOWN"),
           sum(sig_hy_f4$Regulation=="DOWN"),
           sum(sig_f0_f4$Regulation=="DOWN"))
)

# -----------------------------
# 4. Volcano plot
# -----------------------------
plot_volcano <- function(df, title="Volcano Plot") {
  res_df <- df %>%
    mutate(
      group = "NO",
      group = ifelse(p_value <= 0.05 & Diff >= 0.585, "UP", group),
      group = ifelse(p_value <= 0.05 & Diff <= -0.585, "DOWN", group),
      protein_name = Entry_Name
    )
  
  counts_table <- table(res_df$group)
  n_up <- counts_table["UP"]
  n_down <- counts_table["DOWN"]
  n_no <- counts_table["NO"]
  
  top_up <- res_df %>% filter(group=="UP") %>% arrange(p_value) %>% slice_head(n=15)
  top_down <- res_df %>% filter(group=="DOWN") %>% arrange(p_value) %>% slice_head(n=15)
  top_proteins <- bind_rows(top_up, top_down)
  
  ggplot(res_df, aes(x=Diff, y=-log10(p_value), color=group)) +
    geom_point(alpha=0.6, size=1.5) +
    geom_text_repel(
      data=top_proteins,
      aes(label=protein_name),
      size=3, box.padding=0.4, point.padding=0.3,
      segment.color="black", max.overlaps=Inf
    ) +
    scale_color_manual(
      values=c("DOWN"="blue","NO"="grey","UP"="red"),
      labels=c(
        paste0("DOWN (", n_down, ")"),
        paste0("NO (", n_no, ")"),
        paste0("UP (", n_up, ")")
      )
    ) +
    geom_vline(xintercept=c(-0.585,0.585), linetype="dashed") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    theme_minimal() +
    labs(title=title, x="Log2 Fold Change", y="-log10(p-value)", color="")
}

# Exemple de volcano plot Hy vs F0
plot_volcano(hy_f0_sel, title="F0 vs HY")
plot_volcano(hy_f4_sel, title="F4 vs HY")
plot_volcano(f0_f4_sel, title="F4 vs F0")


# -----------------------------
# 5. Heatmap combinant les 3 comparacions
# -----------------------------

# 1. Crear llista única de proteïnes
all_proteins <- unique(c(sig_hy_f0$Entry_Name,
                         sig_hy_f4$Entry_Name,
                         sig_f0_f4$Entry_Name))

# 2. Crear matriu inicial amb zeros
expr_matrix_prot <- matrix(0, nrow=length(all_proteins), ncol=3)
rownames(expr_matrix_prot) <- all_proteins
colnames(expr_matrix_prot) <- c("F0 vs HY","F4 vs HY","F4 vs F0")

# 3. Omplir matriu amb valors Diff
expr_matrix_prot[match(sig_hy_f0$Entry_Name, rownames(expr_matrix_prot)), "F0 vs HY"] <- sig_hy_f0$Diff
expr_matrix_prot[match(sig_hy_f4$Entry_Name, rownames(expr_matrix_prot)), "F4 vs HY"] <- sig_hy_f4$Diff
expr_matrix_prot[match(sig_f0_f4$Entry_Name, rownames(expr_matrix_prot)), "F4 vs F0"] <- sig_f0_f4$Diff

# 4. Substituir NA per 0
expr_matrix_prot[is.na(expr_matrix_prot)] <- 0

# 5. Seleccionar top 100 proteïnes amb major canvi absolut
top_proteins <- names(sort(apply(abs(expr_matrix_prot), 1, max), decreasing=TRUE))[1:100]
expr_top_prot <- expr_matrix_prot[top_proteins, ]

# 6. Calcular Z-score per fila
expr_z_prot <- t(scale(t(expr_top_prot)))

# 7. Generar heatmap
pheatmap(expr_z_prot,
         color=colorRampPalette(c("blue","white","red"))(50),
         border_color=NA,
         fontsize_row=8,
         fontsize_col=12,
         angle_col=0,
         main="Top 100 Proteïnes")



# -----------------------------
# 6. Preparació GO enrichment (opcional)
# -----------------------------
library(clusterProfiler)
library(org.Hs.eg.db)

# Exemple amb proteïnes diferencials UP de HY-F0
proteins_up <- sig_hy_f0$Protein_ID  # utilitzem Protein ID (UniProt)
proteins_up

# Convertir UniProt a Entrez IDs
entrez_ids <- bitr(proteins_up,
                   fromType = "UNIPROT",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
entrez_ids

go_bp <- enrichGO(gene         = entrez_ids$ENTREZID,
                  OrgDb        = org.Hs.eg.db,
                  keyType      = "ENTREZID",
                  ont          = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable     = TRUE)  # converteix Entrez a noms llegibles

head(go_bp)

kegg_enrich <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                          organism     = "hsa",
                          pvalueCutoff = 0.05)

kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db)

head(kegg_enrich)

library(ReactomePA)
library(enrichplot)

reactome_enrich <- enrichPathway(gene = entrez_ids$ENTREZID,
                                 organism = "human",
                                 pvalueCutoff = 0.05,
                                 readable = TRUE)

reactome_enrich
# Plots
barplot(reactome_enrich, showCategory=20, title="Reactome Pathway Enrichment")
dotplot(reactome_enrich, showCategory=20, title="Reactome Pathway Enrichment")


# Barplot GO
barplot(go_bp, showCategory=20, title="GO BP Enrichment - HY vs F0")


# Map Protein ID a Entrez
library(clusterProfiler)
library(org.Hs.eg.db)

sig_hy_f0
up_ids <- sig_hy_f0 %>%
  filter(Regulation == "UP") %>%
  pull(Protein_ID)
length(up_ids)
up_ids

down_ids <- sig_hy_f0 %>%
  filter(Regulation == "DOWN") %>%
  pull(Protein_ID)


entrez_up <- bitr(up_ids,
                  fromType="UNIPROT",
                  toType="ENTREZID",
                  OrgDb=org.Hs.eg.db)
nrow(entrez_up)

entrez_down <- bitr(down_ids,
                    fromType="UNIPROT",
                    toType="ENTREZID",
                    OrgDb=org.Hs.eg.db)

go_up <- enrichGO(
  gene = entrez_up$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "MF",
  readable = TRUE
)
head(go_up)
nrow(go_up)

go_down <- enrichGO(
  gene = entrez_down$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "MF",
  readable = TRUE
)

df_up <- as.data.frame(go_up)
df_down <- as.data.frame(go_down)
head(df_up)
head(df_down)
length(df_up)
length(df_down)

df_up$group <- "Up"
df_down$group <- "Down"
length(df_down)
length(df_up)

df_up
df_up <- na.omit(df_up)
df_down <- na.omit(df_down)
length(df_up)
df_up <- head(df_up, 10)
df_down <- df_down[2:11, ]

df <- rbind(df_up, df_down)

df$value <- -log10(df$p.adjust)
df

# posar els DOWN a la dreta
df$value[df$group == "Down"] <- df$value[df$group == "Down"] * -1

df$term <- paste0(df$Description, " (", df$ID, ")")
library(stringr)

# Convertir la primera lletra a majúscula
df$term <- paste0(str_to_upper(substr(df$Description, 1, 1)), substr(df$Description, 2, nchar(df$Description)), " (", df$ID, ")")


df
colnames(df)
head(df)
# Eliminar files amb NA o NANA
df <- df[!is.na(df$term) & df$term != "NANA (NA)", ]
# Ordena els termes segons value
df$term <- factor(df$term, levels = df$term[order(df$value)])

offset <- 0.15

library(ggplot2)

ggplot(df, aes(x=value, y=term, fill=group)) +
  geom_col(width=0.7) +
  geom_vline(xintercept=0, size=0.6) +
  
  geom_text(
    aes(
      label = term,
      x = ifelse(group == "Up", 0 - offset, 0 + offset)
    ),
    hjust = ifelse(df$group == "Up", 1, 0),
    size = 3
  ) +
  
  scale_fill_manual(values=c("Up"="red", "Down"="blue")) +
  labs(x="-log10(adj p-value)", y=NULL, title="F0-F1 vs HY") +
  
  theme_classic() +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.line.y = element_blank(),
    legend.position="none"
  )

