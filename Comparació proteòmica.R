# Llibreries
library(dplyr)
library(tidyr)
library(readr)
library(readxl)

# 1. LLEGIR FITXERS
vesicles <- read_excel("C:/Users/Usuario/Documents/TÀNIA/UOC Master/TFM/proteomica TA+EVs/cluster_top_proteins_stage.xlsx")
vesicles_DEPs <- read.csv("C:/Users/Usuario/Documents/TÀNIA/UOC Master/TFM/proteomica TA+EVs/DEP_protein_sets.csv", sep=',')
adipose <- read_csv("cluster_top_proteins_stage-TA.csv")
adipose_DEPs <- read.csv("DEP_protein_sets.csv", sep = ';')

head(vesicles)
head(adipose)
head(adipose_DEPs)
head(vesicles_DEPs)
colnames(adipose)

# 2. FUNCIÓ PER NETEJAR DADES
clean_proteins <- function(df, condition) {
  df %>%
    dplyr::select(all_of(condition)) %>%
    mutate(Protein = .[[condition]]) %>%
    filter(!is.na(Protein)) %>%
    separate_rows(Protein, sep = ",") %>%   # separa "A, B, C"
    mutate(Protein = trimws(Protein)) %>%   # elimina espais
    distinct()
}

# 3. CREAR LLISTES PER ESTADI

# ADIPÓS
adipose_healthy <- clean_proteins(adipose, "Healthy")
adipose_F0      <- clean_proteins(adipose, "F0")
adipose_F4      <- clean_proteins(adipose, "F4")

# VESÍCULES
vesicles_healthy <- clean_proteins(vesicles, "Healthy")
vesicles_F0      <- clean_proteins(vesicles, "F0")
vesicles_F4      <- clean_proteins(vesicles, "F4")


# 3.1 Opcio agafar les llistes ja separades
adipose_healthy <- scan("HY_proteins.txt", what = "character")
adipose_F0      <- scan("F0_proteins.txt", what = "character")
adipose_F4      <- scan("F4_proteins.txt", what = "character")

vesicles_healthy <- scan("HY_proteins.txt",  what = "character")
vesicles_F0      <- scan("F0_proteins.txt",  what = "character")
vesicles_F4      <- scan("F4_proteins.txt",  what = "character")

nrow(vesicles_F0)

#  4. TROBAR COMUNES Per estadi entre TA i TA+EVs

common_healthy <- intersect(adipose_healthy, vesicles_healthy)
common_healthy
common_F0      <- intersect(adipose_F0, vesicles_F0)
common_F4      <- intersect(adipose_F4, vesicles_F4)

# diagrama de venn per estadi 
library(VennDiagram)
library(grid)

# F0
venn_F4 <- venn.diagram(
  x = list(
    Adipose = adipose_F4,
    Vesicles = vesicles_F4
  ),
  filename = NULL,
  fill = c("blue", "red"),
  alpha = 0.5,
  main = "F4: Adipose vs Vesicles"
)

venn_F4


# Mirar si les proteines common, estan en els DEPs TA
HY_F0_UP <- na.omit(adipose_DEPs$HY_vs_F0_UP)
HY_F0_DOWN <- na.omit(adipose_DEPs$HY_vs_F0_DOWN)

HY_F4_UP <- na.omit(adipose_DEPs$HY_vs_F4_UP)
HY_F4_DOWN <- na.omit(adipose_DEPs$HY_vs_F4_DOWN)

F0_F4_UP <- na.omit(adipose_DEPs$F0_vs_F4_UP)
F0_F4_DOWN <- na.omit(adipose_DEPs$F0_vs_F4_DOWN)

# Mirar si les proteines common, estan en els DEPs TA+Evs
HY_F0_UP <- na.omit(vesicles_DEPs$HY_vs_F0_UP)
HY_F0_DOWN <- na.omit(vesicles_DEPs$HY_vs_F0_DOWN)

HY_F4_UP <- na.omit(vesicles_DEPs$HY_vs_F4_UP)
HY_F4_DOWN <- na.omit(vesicles_DEPs$HY_vs_F4_DOWN)

F0_F4_UP <- na.omit(vesicles_DEPs$F0_vs_F4_UP)
F0_F4_DOWN <- na.omit(vesicles_DEPs$F0_vs_F4_DOWN)


common_in_HY <- list(
  HY_vs_F0_UP   = intersect(common_healthy, HY_F0_UP),
  HY_vs_F0_DOWN = intersect(common_healthy, HY_F0_DOWN),
  
  HY_vs_F4_UP   = intersect(common_healthy, HY_F4_UP),
  HY_vs_F4_DOWN = intersect(common_healthy, HY_F4_DOWN)
)

common_in_HY

common_in_F0 <- list(
  F0_vs_F4_UP   = intersect(common_F0, F0_F4_UP),
  F0_vs_F4_DOWN = intersect(common_F0, F0_F4_DOWN),
  
  HY_vs_F0_UP   = intersect(common_F0, HY_F0_UP),
  HY_vs_F0_DOWN = intersect(common_F0, HY_F0_DOWN)
)

common_in_F0

common_in_F4 <- list(
  HY_vs_F4_UP   = intersect(common_F4, HY_F4_UP),
  HY_vs_F4_DOWN = intersect(common_F4, HY_F4_DOWN),
  
  F0_vs_F4_UP   = intersect(common_F4, F0_F4_UP),
  F0_vs_F4_DOWN = intersect(common_F4, F0_F4_DOWN)
)

common_in_F4


# fer grafica conjunta
library(dplyr)
library(tidyr)

# 1. Unir DEPs en un sol format llarg

get_dep_df <- function(dep_list, dataset_name) {
  bind_rows(
    data.frame(Protein = dep_list$HY_vs_F0_UP,   Comparison = "F0_vs_HY", Regulation = "UP"),
    data.frame(Protein = dep_list$HY_vs_F0_DOWN, Comparison = "F0_vs_HY", Regulation = "DOWN"),
    
    data.frame(Protein = dep_list$HY_vs_F4_UP,   Comparison = "F4_vs_HY", Regulation = "UP"),
    data.frame(Protein = dep_list$HY_vs_F4_DOWN, Comparison = "F4_vs_HY", Regulation = "DOWN"),
    
    data.frame(Protein = dep_list$F0_vs_F4_UP,   Comparison = "F4_vs_F0", Regulation = "UP"),
    data.frame(Protein = dep_list$F0_vs_F4_DOWN, Comparison = "F4_vs_F0", Regulation = "DOWN")
  ) %>%
    mutate(Dataset = dataset_name)
}

adipose_long  <- get_dep_df(adipose_DEPs, "Adipose")
vesicles_long <- get_dep_df(vesicles_DEPs, "Vesicles")

all_deps <- bind_rows(adipose_long, vesicles_long)

# 2. Unir proteïnes comunes

all_common <- unique(c(
  common_healthy,
  common_F0,
  common_F4
))


# 3. Filtrar només proteïnes comunes que són DEPs

common_dep_hits <- all_deps %>%
  filter(Protein %in% all_common)

nrow(all_deps)
nrow(common_dep_hits)


# 4. Resum final (agrupat)


final_table <- common_dep_hits %>%
  group_by(Protein, Comparison, Regulation) %>%
  summarise(
    Datasets = paste(unique(Dataset), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(Protein)

final_table
nrow(final_table)

write.csv(final_table, "common_proteins_DEPs.csv", row.names = FALSE)

library(ggplot2)

df_part1 <- final_table[1:110, ]
df_part2 <- final_table[111:232, ]
df_shared <- final_table %>%
  filter(grepl("Adipose, Vesicles", Datasets))

ggplot(df_part1, aes(x = Comparison, y = Protein)) +
  geom_point(aes(color = Regulation, shape = Datasets), size = 3) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6)) +
  labs(title = "Common DEPs across comparisons",
       x = "Comparison",
       y = "Protein")


# correlacio TA vs EVs
evs_hy_f0 <- read_excel("ATEVsn_resultados_proteómica_final.xlsx", sheet = "Healthy vs FO")
evs_hy_f0
evs_hy_f4 <- read_excel("ATEVsn_resultados_proteómica_final.xlsx", sheet = "Healthy vs F4")
evs_f0_f4 <- read_excel("ATEVsn_resultados_proteómica_final.xlsx", sheet = "F0 vs F4")

TA_hy_f0 <- read_excel("AT_REANALYSED.xlsx", sheet = "Healthy vs F0-F1")
TA_hy_f0
TA_hy_f4 <- read_excel("AT_REANALYSED.xlsx", sheet = "Healthy vs F3-F4")
TA_f0_f4 <- read_excel("AT_REANALYSED.xlsx", sheet = "F0-F1 vs F3-F4")

# Seleccionar columnes rellevants
adipose_sel <- TA_hy_f4 %>%
  dplyr::select(Protein = `T: Entry Name`, F0_Healthy = `F3-F4/Healthy`)

vesicles_sel <- evs_hy_f4 %>%
  dplyr::select(Protein = `T: Protein.Names`, F0_Healthy = `F4/Healthy`)

merged <- inner_join(adipose_sel, vesicles_sel, by = "Protein", suffix = c("_adipose", "_vesicles"))

merged

merged <- merged %>%
  mutate(
    log2_adipose = log2(F0_Healthy_adipose),
    log2_vesicles = log2(F0_Healthy_vesicles)
  )

cor_value  <-  cor(merged$log2_adipose, merged$log2_vesicles, method = "spearman")

ggplot(merged, aes(x = log2_adipose, y = log2_vesicles)) +
  geom_point(alpha = 0.3) +
  
  # destacar comunes
  geom_point(data = merged, size = 2) +
  
  geom_smooth(method = "lm") +
  
  theme_minimal() +
  labs(
    title = paste0("Correlation F0 vs Healthy (r = ", round(cor_value, 3), ")"),
    x = "log2(F0/Healthy) Adipose",
    y = "log2(F0/Healthy) Vesicles"
  )


library(ggplot2)

ggplot(merged, aes(x = log2_adipose, y = log2_vesicles)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal() +
  labs(
    title = paste0("Correlation F0 vs Healthy (r = ", round(cor_value, 3), ")"),
    x = "log2(F0/Healthy) Adipose",
    y = "log2(F0/Healthy) Vesicles"
  )
