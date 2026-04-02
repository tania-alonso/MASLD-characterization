library(dplyr)
library(pheatmap)
library(tidyr)
library(biomaRt)
library(GO.db)

getwd()
setwd("C:/Users/Usuario/Documents/TÀNIA/UOC Master/TFM/proteomica TA+EVs")

# llegir dades AT
# data <- readxl::read_excel("AT_REANALYSED.xlsx", sheet = "SEL INFO")
# str(data)
# head(data)
# 
# # convertir comes decimals a punt
# sample_cols <- c("P393","P426","P497","P499","P538",
#                  "P376","P556 AT","P366","1468 AT","P360",
#                  "P495","P396","P409","P397","P401")
# 
# data <- data %>%
#   mutate(across(all_of(sample_cols),
#                 ~as.numeric(gsub(",", ".", .))))
# 
# str(data)
# # definir mostres per estadi
# healthy_samples <- c("P393","P426","P497","P499","P538")
# F0_samples <- c("P376","P556 AT","P366","1468 AT","P360")
# F4_samples <- c("P495","P396","P409","P397","P401")

# llegir dades AT+EVs
data <- readxl::read_excel("ATEVsn_resultados_proteómica_final.xlsx", sheet = "Proteins")
str(data)
head(data)

# convertir comes decimals a punt
sample_cols <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6", "Sample7", "Sample8", 
                 "Sample10", "Sample11", "Sample12", "Sample13", "Sample14", "Sample15", "Sample16",
                 "Sample17", "Sample18", "Sample19", "Sample22",
                 "Sample23", "Sample24", "Sample25", "Sample26",  "Sample27", "Sample28", "Sample29", "Sample30", "Sample31")

data <- data %>%
  mutate(across(all_of(sample_cols),
                ~as.numeric(gsub(",", ".", .))))

str(data)
# definir mostres per estadi
healthy_samples <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6", "Sample7", "Sample8", "Sample10")
F0_samples <- c("Sample11", "Sample12", "Sample13", "Sample14", "Sample15", "Sample16", "Sample17",
                "Sample18", "Sample19", "Sample22", "Sample23",  "Sample24", "Sample25")
F4_samples <- c("Sample26",  "Sample27", "Sample28", "Sample29", "Sample30", "Sample31")


# calcular mitjana per estadi
data$Healthy_mean <- rowMeans(data[,healthy_samples], na.rm=TRUE)
data$F0_mean <- rowMeans(data[,F0_samples], na.rm=TRUE)
data$F4_mean <- rowMeans(data[,F4_samples], na.rm=TRUE)

data
colnames(data)
# crear matriu d'expressió
expr_matrix <- data %>%
  dplyr::select(`Entry Name`, Healthy_mean, F0_mean, F4_mean)

expr_matrix <- as.data.frame(expr_matrix)

rownames(expr_matrix) <- expr_matrix$`Entry Name`

expr_matrix <- expr_matrix %>%
  dplyr::select(Healthy_mean, F0_mean, F4_mean)

expr_matrix <- as.matrix(expr_matrix)
expr_matrix

# seleccionar proteïnes amb més abundacia
top_proteins <- names(sort(apply(expr_matrix,1,sd), decreasing=TRUE))[1:100]
top_proteins
expr_top <- expr_matrix[top_proteins,]

# z-score
expr_z <- t(scale(t(expr_top)))
expr_z <- expr_z[complete.cases(expr_z),]

# heatmap
pheatmap(expr_z,
         color=colorRampPalette(c("blue","white","red"))(50),
         border_color=NA,
         fontsize_row=7,
         fontsize_col=12,
         angle_col=0,
         main="Proteïnes característiques per estadi")

# -------------------------
# TOP proteïnes per estadi
# -------------------------

expr_df <- as.data.frame(expr_matrix)
expr_df
expr_matrix

expr_df$Healthy_z <- scale(expr_df$Healthy_mean)
expr_df$F0_z <- scale(expr_df$F0_mean)
expr_df$F4_z <- scale(expr_df$F4_mean)

top_healthy <- rownames(expr_df)[order(-expr_df$Healthy_z)][1:100]
top_F0 <- rownames(expr_df)[order(-expr_df$F0_z)][1:100]
top_F4 <- rownames(expr_df)[order(-expr_df$F4_z)][1:100]

top_table <- data.frame(
  Healthy = top_healthy,
  F0 = top_F0,
  F4 = top_F4
)

top_table

# write.csv(top_table, "top_proteins_per_stage.csv", row.names = FALSE)

# -------------------------
# Preparar proteïnes
# -------------------------
protein_map <- data %>%
  dplyr::select(`Entry Name`, `Protein ID`) %>%
  distinct()

top_long <- top_table %>%
  pivot_longer(cols = everything(),
               names_to = "Stage",
               values_to = "Entry Name")

top_long <- top_long %>%
  left_join(protein_map, by="Entry Name")

top_long


# -------------------------
# obtenir funcions GO
# -------------------------

ensembl <- useMart("ensembl",
                   dataset="hsapiens_gene_ensembl",
                   host="https://www.ensembl.org")

annotations <- getBM(
  attributes = c("uniprotswissprot","hgnc_symbol","description","go_id","name_1006"),
  filters = "uniprotswissprot",
  values = unique(top_long$`Protein ID`),
  mart = ensembl
)

annotations
# -------------------------
# unir proteïnes amb funció
# -------------------------

top_func <- top_long %>%
  left_join(annotations,
            by=c("Protein ID"="uniprotswissprot"))

top_func
# -------------------------
# agrupar per funció
# -------------------------
func_table <- top_func %>%
  filter(!is.na(name_1006), name_1006 != "") %>%
  group_by(Stage, name_1006) %>%
  summarise(
    GO_IDs   = paste(unique(go_id[!is.na(go_id) & go_id != ""]), collapse = ", "), # al costat del name_1006
    Proteins = paste(unique(`Protein ID`), collapse = ", "),
    .groups = "drop"
  ) %>%
  dplyr::select(Stage, name_1006, GO_IDs, Proteins)  # ordre de columnes

func_table

# Crear taula final amb funcions com a files i estats com a columnes
final_table <- func_table %>%
  pivot_wider(
    names_from = Stage,
    values_from = Proteins
  ) %>%
  arrange(name_1006)

final_table
write.csv(final_table, "top_functions_per_stage.csv", row.names = FALSE)

# Substituir Protein IDs per Entry Names a func_table
func_table_named <- func_table %>%
  # Separar les proteïnes per fila per fer el mapatge
  tidyr::separate_rows(Proteins, sep = ",\\s*") %>%
  left_join(protein_map, by = c("Proteins" = "Protein ID")) %>%
  # Si no troba entry name, deixar el Protein ID original
  mutate(Proteins = ifelse(is.na(`Entry Name`), Proteins, `Entry Name`)) %>%
  dplyr::select(-`Entry Name`) %>%
  # Tornar a unir per funció i estadi amb guions
  group_by(Stage, name_1006, GO_IDs) %>%
  summarise(Proteins = paste(unique(Proteins), collapse = " - "), .groups = "drop")

func_table_named


# -------------------------
# Filtrar només l'estadi Healthy
# -------------------------
final_table_healthy <- func_table_named %>%
  filter(Stage == "Healthy") %>%
  dplyr::select(name_1006, GO_IDs, Proteins) %>%
  # comptar nombre de proteïnes per funció
  mutate(
    NumProteins_Healthy = sapply(Proteins, function(x) {
      if(is.na(x) | x=="") 0 else length(unlist(strsplit(x, " - ")))
    })
  ) %>%
  arrange(desc(NumProteins_Healthy))

final_table_healthy
# Exportar
write.csv(final_table_healthy, "top_functions+prot_HY.csv", row.names = FALSE)

# -------------------------
# Filtrar només l'estadi F0
# -------------------------
final_table_F0 <- func_table_named %>%
  filter(Stage == "F0") %>%
  dplyr::select(name_1006, GO_IDs, Proteins) %>%
  # comptar nombre de proteïnes per funció
  mutate(
    NumProteins_F0 = sapply(Proteins, function(x) {
      if(is.na(x) | x=="") 0 else length(unlist(strsplit(x, " - ")))
    })
  ) %>%
  arrange(desc(NumProteins_F0))

final_table_F0
write.csv(final_table_F0, "top_functions+prot_F0.csv", row.names = FALSE)

# -------------------------
# Filtrar només l'estadi F4
# -------------------------
final_table_F4 <- func_table_named %>%
  filter(Stage == "F4") %>%
  dplyr::select(name_1006, GO_IDs, Proteins) %>%
  # comptar nombre de proteïnes per funció
  mutate(
    NumProteins_F4 = sapply(Proteins, function(x) {
      if(is.na(x) | x=="") 0 else length(unlist(strsplit(x, " - ")))
    })
  ) %>%
  arrange(desc(NumProteins_F4))

final_table_F4
write.csv(final_table_F4, "top_functions+prot_F4.csv", row.names = FALSE)



# afegir categoria a les taules anteriors 
library(dplyr)
library(tidyr)
library(GO.db)
library(AnnotationDbi)

# healthy

# 1. Separar GO IDs
HY_expanded <- final_table_healthy %>%
  separate_rows(GO_IDs, sep = " - ") %>%
  mutate(GO_IDs = trimws(GO_IDs))

# 2. Obtenir categories GO
go_categories <- AnnotationDbi::select(
  GO.db,
  keys = unique(HY_expanded$GO_IDs),
  columns = c("GOID","ONTOLOGY"),
  keytype = "GOID"
)

# 3. Fer join amb la teva taula
HY_with_cat <- HY_expanded %>%
  left_join(go_categories, by = c("GO_IDs" = "GOID"))
HY_with_cat
nrow(HY_with_cat)

# # 4. Filtrar només Biological Process i Molecular Function
# HY_filtered <- HY_with_cat %>%
#   filter(ONTOLOGY %in% c("MF"))
# 
# HY_filtered
# nrow(HY_filtered)

# Exportar
write.csv(HY_with_cat, "top_functions+prot_HY.csv", row.names = FALSE)


# F0
# 1. Separar GO IDs
f0_expanded <- final_table_F0 %>%
  separate_rows(GO_IDs, sep = " - ") %>%
  mutate(GO_IDs = trimws(GO_IDs))

# 2. Obtenir categories GO
go_categories <- AnnotationDbi::select(
  GO.db,
  keys = unique(f0_expanded$GO_IDs),
  columns = c("GOID","ONTOLOGY"),
  keytype = "GOID"
)

# 3. Fer join amb la teva taula
f0_with_cat <- f0_expanded %>%
  left_join(go_categories, by = c("GO_IDs" = "GOID"))


# # 4. Filtrar només Biological Process i Molecular Function
# f0_filtered <- f0_with_cat %>%
#   filter(ONTOLOGY %in% c("MF"))
# 
# f0_filtered
# nrow(f0_filtered)

# Exportar
write.csv(f0_with_cat, "top_functions+prot_F0.csv", row.names = FALSE)

# F4
# 1. Separar GO IDs
f4_expanded <- final_table_F4 %>%
  separate_rows(GO_IDs, sep = " - ") %>%
  mutate(GO_IDs = trimws(GO_IDs))

# 2. Obtenir categories GO
go_categories <- AnnotationDbi::select(
  GO.db,
  keys = unique(f4_expanded$GO_IDs),
  columns = c("GOID","ONTOLOGY"),
  keytype = "GOID"
)

# 3. Fer join amb la teva taula
f4_with_cat <- f4_expanded %>%
  left_join(go_categories, by = c("GO_IDs" = "GOID"))
f4_with_cat
nrow(f4_with_cat)

# # 4. Filtrar només Biological Process i Molecular Function
# f4_filtered <- f4_with_cat %>%
#   filter(ONTOLOGY %in% c("MF"))
# 
# f4_filtered
# nrow(f4_filtered)


# Exportar
write.csv(f4_with_cat, "top_functions+prot_F4.csv", row.names = FALSE)





# fer analisis per CLUSTERS

data$Stage <- c("Healthy","F0","F4")[
  max.col(data[,c("Healthy_mean","F0_mean","F4_mean")], ties.method = "first")
]

data
colnames(data)

healthy_prot <- data[data$Stage=="Healthy",]
F0_prot <- data[data$Stage=="F0",]
F4_prot <- data[data$Stage=="F4",]

data$Stage <- unlist(data$Stage)
table(data$Stage)

healthy_top <- healthy_prot[order(-healthy_prot$Healthy_mean),]
F0_top <- F0_prot[order(-F0_prot$F0_mean),]
F4_top <- F4_prot[order(-F4_prot$F4_mean),]


healthy_names <- healthy_top$`Protein ID` # `Entry Name` #
length(healthy_names)
F0_names <- F0_top$`Protein ID` # `Entry Name` # 
F4_names <- F4_top$`Protein ID` # `Entry Name` # 
F4_names

max_len <- max(length(healthy_names), length(F0_names), length(F4_names))

healthy_names <- c(healthy_names, rep(NA, max_len - length(healthy_names)))
F0_names <- c(F0_names, rep(NA, max_len - length(F0_names)))
F4_names <- c(F4_names, rep(NA, max_len - length(F4_names)))
F4_names

healthy_prot
stage_table <- data.frame(
  Healthy = healthy_names,
  F0 = F0_names,
  F4 = F4_names
)

stage_table

# write.csv(stage_table, "top_proteins_stage_cluster_ids.csv", row.names = FALSE)

healthy_df <- data.frame(
  ProteinID = healthy_names,
  Stage = "Healthy"
)

healthy_df

F0_df <- data.frame(
  ProteinID = F0_names,
  Stage = "F0"
)
F0_df

F4_df <- data.frame(
  ProteinID = F4_names,
  Stage = "F4"
)
F4_df

expand_proteins <- function(protein_vector, stage_name) {
  proteins <- unlist(strsplit(protein_vector, ";"))
  proteins <- proteins[proteins != ""]
  data.frame(ProteinID = proteins, Stage = stage_name)
}

healthy_df <- expand_proteins(healthy_names, "Healthy")
F0_df <- expand_proteins(F0_names, "F0")
F4_df <- expand_proteins(F4_names, "F4")


# -------------------------
# obtenir funcions GO
# -------------------------
ensembl <- useMart("ensembl",
                   dataset="hsapiens_gene_ensembl",
                   host="https://www.ensembl.org")


annotations <- getBM(
  attributes = c("uniprotswissprot","hgnc_symbol","description","go_id","name_1006"),
  filters = "uniprotswissprot",
  values = unique(healthy_df$ProteinID),
  mart = ensembl
)

annotations
# -------------------------
# unir proteïnes amb funció
# -------------------------

top_func <- healthy_df %>%
  left_join(annotations,
            by=c("ProteinID"="uniprotswissprot"))

top_func
colnames(top_func)

# -------------------------
# Agrupar per funció
# -------------------------
func_tables <- top_func %>%
  filter(!is.na(name_1006) & name_1006 != "") %>%
  
  # crear columna amb nom llegible de la proteïna
  mutate(ProteinName = ifelse(is.na(hgnc_symbol) | hgnc_symbol == "",
                              ProteinID,
                              hgnc_symbol)) %>%
  
  group_by(Stage, name_1006) %>%
  summarise(
    GO_IDs = paste(unique(go_id[!is.na(go_id) & go_id != ""]), collapse = " - "),
    Proteins = paste(unique(ProteinName), collapse = " - "),
    NumProteins = length(unique(ProteinName)),
    .groups = "drop"
  )

# Opcional: ordenar per nombre de proteïnes
cluster_table <- func_tables %>% arrange(desc(NumProteins))
cluster_table


# afegir categoria a les taules anteriors 

# 1. Separar GO IDs
df_expanded <- cluster_table %>%
  separate_rows(GO_IDs, sep = " - ") %>%
  mutate(GO_IDs = trimws(GO_IDs))

# 2. Obtenir categories GO
go_categories <- AnnotationDbi::select(
  GO.db,
  keys = unique(df_expanded$GO_IDs),
  columns = c("GOID","ONTOLOGY"),
  keytype = "GOID"
)

# 3. Fer join amb la teva taula
HY_with_cat <- df_expanded %>%
  left_join(go_categories, by = c("GO_IDs" = "GOID"))
HY_with_cat
nrow(HY_with_cat)

# Exportar
write.csv(HY_with_cat, "cluster_top_functions+prot_HY.csv", row.names = FALSE)

