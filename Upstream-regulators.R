# trobar upstream regulators dels DEGS del bulk RNA-seq
# nomes te sentit buscar dels DEGS UP

# Carregar llibreries
library(readxl)
library(DESeq2)
library(dorothea)
library(viper)
library(tibble)
library(dplyr)
library(openxlsx)

# Llegir el fitxer de counts
counts_df <- read_excel("C:/Users/Usuario/Documents/TÀNIA/UOC Master/TFM/Obtenció DEGs/test-counts.xlsx")  # gens a les files
head(counts_df)
# Crear matriu d'expressió

# Seleccionar només les columnes de mostres
sample_cols <- c("CTL_1","CTL_2","NASH_1","NASH_2",
                 "P548_F0","P554_F0","P599_FO","P511_HY","P579_HY","P585_HY")

# ordenar mostres 
expression_matrix <- as.matrix(counts_df[, sample_cols])
rownames(expression_matrix) <- counts_df$gene_name
expression_matrix <- as.matrix(expression_matrix)

#  Crear colData (informació de mostres)
colData <- data.frame(
  sample = sample_cols,
  condition = c("CTL", "CTL", "NASH", "NASH", "F0", "F0", "F0", "HY", "HY", "HY")
)
rownames(colData) <- colData$sample
# ordenar mostres 
expression_matrix <- expression_matrix[, rownames(colData)]

#  Normalització amb DESeq2
dds <- DESeqDataSetFromMatrix(countData = expression_matrix,
                              colData = colData,
                              design = ~ condition)
vsd <- vst(dds)
expression_matrix <- assay(vsd)
expression_matrix


#  Carregar regulons de TFs (DoRothEA)
#  per cada TF mirar quins targets té

data(dorothea_hs, package = "dorothea")

# Filtrem per confiança i convertim directament al format VIPER
regulons_df <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B", "C"))

regulons <- regulons_df %>%
  dorothea::df2regulon()

regulons_df


#  Calcul de l'activitat dels TFs amb VIPER 
# Mirar els targets de cada TF i comprova si estan globalment activats 

tf_activity <- viper(expression_matrix, regulons, 
                     nes = TRUE, 
                     method = "none", # "none" si la matriu ja és normalitzada/vst
                     minsize = 5, 
                     verbose = FALSE)

# 1. Preparem el dataframe d'activitat (passem rownames a columna)
tf_activity_to_save <- as.data.frame(tf_activity) %>%
  rownames_to_column(var = "TF")

# 2. Creem una llista amb els objectes que volem guardar
# Cada element de la llista serà una pestanya de l'Excel
llista_exportar <- list(
  "Activitat_TF" = tf_activity_to_save,
  "Regulons_DoRothEA" = regulons_df
)

# 3. Guardem a la carpeta que vulguis
write.xlsx(llista_exportar, file = "Resultats_TF_Inference.xlsx")



#  Explorar resultats
# Els TFs més actius o més diferencialment reguladors
tf_activity_df <- as.data.frame(tf_activity)
head(tf_activity_df)

# ordenar pel valor d'activitat per veure reguladors més actius
# calcula la mitjana de l'activitat dels TFs de totes les mostres (estadis).
# es a dir de hy, f0, f4 i nash, llavors la mitjana pot ser 0. si puja en un i baixa en un altre. 
tf_activity_df <- tf_activity_df %>% 
  rownames_to_column(var="TF") %>%
  arrange(desc(rowMeans(tf_activity_df)))

# Mostrar els 10 TFs més actius
head(tf_activity_df, 10)
nrow(tf_activity_df)

# heatmap
library(pheatmap)

# Seleccionem els 20 TFs amb més variabilitat (els que més canvien entre mostres)
variança_tf <- apply(tf_activity, 1, var)
top_variat <- names(sort(variança_tf, decreasing = TRUE))[1:50]

# Fem el heatmap
pheatmap(tf_activity[top_variat, ], 
         main = "Activitat dels TFs (Top 20 més variables)",
         clustering_distance_rows = "correlation",
         scale ="row")

# mirar el canvi d'activitat del TFs entre estadis
library(limma)

# 1. Preparem la matriu de disseny 
design <- model.matrix(~ 0 + condition, data = colData)
colnames(design) <- levels(as.factor(colData$condition))

# 2. Ajustem el model lineal sobre les ACTIVITATS (no sobre els counts)
fit <- lmFit(tf_activity, design)

# 3. Definim la comparació (NASH vs CTL)
# hy vs ctl
#  contrast.matrix <- makeContrasts(HY_vs_CTL = HY - CTL, levels = design)
# f0 vs ctl
#  contrast.matrix <- makeContrasts(F0_vs_CTL = F0 - CTL, levels = design)
# # f0 vs hy
# contrast.matrix <- makeContrasts(F0_vs_HY = F0 - HY, levels = design)
# # nash vs ctl
 contrast.matrix <- makeContrasts(NASH_vs_CTL = NASH - CTL, levels = design)
# # f0 vs nash 
# contrast.matrix <- makeContrasts(F0_vs_NASH = F0 - NASH, levels = design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 4. Obtenim la taula de TFs "upstream" més importants
top_regulators <- topTable(fit2, coef = "NASH_vs_CTL", number = Inf)
top_regulators
nrow(top_regulators)

# 5. Filtrem pels que tenen activitat positiva (Upstream regulators dels teus gens UP)
# agafem pvalue pq padj no hi ha cap significatiu 
upstream_up <- top_regulators %>% 
  filter(logFC > 0 & P.Value < 0.05) %>%
  arrange(desc(logFC))

upstream_up
head(upstream_up)
nrow(upstream_up)

library(ggplot2)

top_regulators$TF <- rownames(top_regulators)
upstream_up$TF <- rownames(upstream_up)

ggplot(top_regulators, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = P.Value < 0.05)) +
  geom_text(data = head(upstream_up, 10), aes(label = TF), vjust = 1, hjust = 1) +
  theme_minimal() +
  labs(title = "Reguladors Upstream (F0 vs NASH)",
       x = "Activitat Diferencial (log2FC)",
       y = "-log10 p-valor")


# mirar targets dels upstream regulatos que vulguis 
# agafo els UR amb logFc mes gran que 0.585
# hy vs ctl: RUNX2, ARNT, SIX2, OTX2, POU5F1, HOXB13 (cap esta a la proteomica)
# f0 vs ctl: no hi ha cap amb log2FC mes gran que 0.585 
# f0 vs hy:  no hi ha cap amb log2FC positiu (mes gran que 0)
# nash vs ctl: ATF4, FOX2L, NFATC1, MEF2A, ARNT, ARNTL, KLF6, EPAS1 (cap esta a la proteomica)
# f0 vs nash:  # KLF5, THAP1 , MAZ (cap esta a la proteomica)


#overpap Tfs amb proteomica TA
proteomics_df <- read_excel("C:/Users/Usuario/Documents/TÀNIA/UOC Master/TFM/proteomica TA/AT_REANALYSED.xlsx", sheet = "SEL INFO")
#overpap Tfs amb proteomica TA+Evs
proteomics_df <- read_excel("C:/Users/Usuario/Documents/TÀNIA/UOC Master/TFM/proteomica TA+EVs/ATEVsn_resultados_proteómica_final.xlsx", sheet = "Proteins")

proteins <- proteomics_df$"Entry Name"

tf_activity_df
tf_list <- tf_activity_df$TF
tf_list

common_tf <- intersect(tf_list, proteins)
common_tf

# selected_tfs <- c("NFIC","STAT2","STAT6",
#                  "STAT3","STAT1","NFKB1")

pheatmap(tf_activity[common_tf, ],
         annotation_col = colData["condition", drop=F],
         scale = "row",
         main = "Activitat dels TFs validats (RNA + proteòmica)")

# analitzar targets
regulons_df
targets_tf <- regulons_df %>%
  filter(tf == "CREB1", mor == 1) %>%
  pull(target)
targets_tf

targets_tf <- intersect(targets_tf, rownames(expression_matrix))

tf_targets_expr <- expression_matrix[targets_tf, ]

# eliminar gens amb NA
tf_targets_expr <- tf_targets_expr[complete.cases(tf_targets_expr), ]

# eliminar gens amb variància 0
tf_targets_expr <- tf_targets_expr[apply(tf_targets_expr, 1, var) != 0, ]

pheatmap(tf_targets_expr, scale = "row")



# overlap amb tfs de transcriptomica degs
degs_trans <- read.csv2("C:/Users/Usuario/Documents/TÀNIA/UOC Master/TFM/Obtenció DEGs/DEGs-UP_NASH_vs_HY.csv")
degs_trans
# overlap amb tfs de transcriptomica 
transcriptomics_df <- read_excel("C:/Users/Usuario/Documents/TÀNIA/UOC Master/TFM/Obtenció DEGs/test-counts.xlsx")
transcriptomics_df

genes <- transcriptomics_df$"gene_name"

common_tf <- intersect(tf_list, genes)
common_tf

counts_tfs <- expression_matrix[rownames(expression_matrix) %in% common_tf, ]

# Calculem l'expressió mitjana de cada gen
abundancia_mitjana <- rowMeans(expression_matrix)

# Mirem on cauen els teus TFs en el rànquing general de tots els gens
ranking_abundancia <- data.frame(
  Gene = names(abundancia_mitjana),
  Abundance = abundancia_mitjana
) %>% arrange(desc(Abundance))
ranking_abundancia

# Afegim una columna per saber si és un dels teus TFs d'interès
ranking_abundancia$Is_TF <- ranking_abundancia$Gene %in% common_tf

library(ggplot2)

# Filtrem només els TFs confirmats per fer el gràfic
top_tfs_plot <- ranking_abundancia %>% filter(Is_TF == TRUE)
top_tfs_plot
ggplot(top_tfs_plot, aes(x = reorder(Gene, Abundance), y = Abundance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Abundància absoluta dels TFs seleccionats",
       subtitle = "Basat en counts normalitzats (VST)",
       x = "Factor de Transcripció",
       y = "Expressió mitjana (log2)")


library(ggplot2)
library(dplyr)
library(ggrepel) # Perquè els noms no s'encavalquin

# 1. Unim l'activitat (VIPER) amb l'abundància (VST)
plot_data <- top_regulators %>%
  rownames_to_column(var = "Gene") %>%
  inner_join(ranking_abundancia, by = "Gene")

# 2. Afegim la informació de si és un DEG (aquí pots usar la llista que ja tenies)
# Suposem que 'common_tf' és la llista de TFs que surten a l'Excel
plot_data <- plot_data %>%
  mutate(Is_DEG = ifelse(Gene %in% common_tf, "DEG Confirmat", "Només Activitat"))

ggplot(plot_data, aes(x = Abundance, y = logFC)) +
  # Punts de fons (tots els TFs)
  geom_point(aes(color = Is_DEG, size = -log10(P.Value)))+  
  # Línies de referència (mitjana i activitat zero)
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  
  # Etiquetem els 15 TFs més importants (per activitat absoluta)
  geom_text_repel(data = plot_data %>% arrange(desc(abs(logFC))) %>% head(15),
                  aes(label = Gene), 
                  box.padding = 0.5, 
                  max.overlaps = Inf,
                  size = 3.5,
                  fontface = "bold") +
  
  # Personalització de colors i estètica
  scale_color_manual(values = c("DEG Confirmat" = "#E41A1C", "Només Activitat" = "#377EB8")) +
  theme_minimal() +
  labs(title = "Paisatge Regulador: Abundància vs. Activitat",
       subtitle = "NASH vs CTL (Mostrant TFs Upstream)",
       x = "Abundància (Expressió mitjana VST)",
       y = "Activitat Diferencial (logFC VIPER)",
       color = "Validació Transcriptòmica") +
  theme(legend.position = "bottom")


