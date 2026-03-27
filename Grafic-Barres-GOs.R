library(ggplot2)
library(dplyr)

# HY vs CTL
# df <- data.frame(
#   term = c(
#     "Multicellular organismal process (GO:0032501)",
#     "Extracellular matrix organization (GO:0043062)",
#     "Extracellular structure organization (GO:0045229)",
#     "Cell adhesion (GO:0007155)",
#     "L-amino acid transport (GO:0015807)",
#     "Collagen fibril organization (GO:0030199)",
#     "Regulation of cell communication (GO:0010646)",
#     "Regulation of cell migration (GO:0030334)",
# 
#     "Lipid metabolic process (GO:0006629)",
#     "Unsaturated fatty acid metabolic process (GO:0033559)",
#     "Estrogen metabolic process (GO:0008210)",
#     "Lipid transport (GO:0006869)",
#     "Digestion (GO:0007586)",
#     "Prostaglandin metabolic process (GO:0006693)",
#     "Monocarboxylic acid metabolic process (GO:0032787)",
#     "Isoprenoid metabolic process (GO:0006720)"
#   ),
# 
#   value = c(
#    4.6662,
#    4.2779,
#    4.2568,
#    2.4415,
#    2.4462,
#    2.3935,
#    2.83933,
#    2.2718,
#   -4.0725,
#   -3.0805,
#   -2.70143,
#   -1.4630,
#   -1.3456,
#   -3.4730,
#   -2.7156,
#   -1.8304
#   )
# )

# F0 vs CTL
# df <- data.frame(
#     term = c(
#       "Response to oxygen levels (GO:0070482)",
#       "Response to hypoxia (GO:0001666)",
#       "Regulation of programmed cell death (GO:0043067)",
#       "Apoptotic process (GO:0006915)",
#       "Mitotic cell cycle (GO:0000278)",
#       "Nuclear division (GO:0000280)",
#       "Cellular response to stimulus (GO:0051716)",
#       "Organic acid transmembrane transport (GO:1903825)",
#       "NADPH regeneration (GO:0006740)",
#       "NADP metabolic process (GO:0006739)",
#       "Cellular response to amino acid stimulus (GO:0071230)",
#       "Response to bacterium (GO:0009617)",
#       "Response to cytokine (GO:0034097)",
#       "Collagen fibril organization (GO:0030199)"
#     ),
# 
#     value = c(
#       -5.1333,
#       -4.3306,
#       -2.6324,
#       -1.5225,
#       -2.2465,
#       -1.3750,
#       -1.9629,
#       1.6797,
#       1.9829,
#       1.8763,
#       1.6737,
#       1.5994,
#       1.5553,
#       1.5617
#     )
#   )


# NASH vs CTL
# df <- data.frame(
#     term = c(
#       "Small molecule metabolic process (GO:0044281)",
#       "Response to xenobiotic stimulus (GO:0009410)",
#       "Alcohol metabolic process (GO:0006066)",
#       "Animal organ regeneration (GO:0031100)",
#       "Fatty acid transport (GO:0015908)",
#       "Regulation of primary metabolic process (GO:0080090)",
#       "Response to stress (GO:0006950)",
#       "Inflammatory response (GO:0006945)",
#       "Response to hypoxia (GO:0001666)",
#       "Cell migration (GO:0016477)",
#       "Cell population proliferation (GO:0008283)",
#       "Angiogenesis (GO:0001525)",
#       "Cholesterol homeostasis"
#     ),
# 
#     value = c(
#       -2.0408,
#       -1.6451,
#       -1.5999,
#       -1.3590,
#       4.0944,
#       4.2290,
#       9.3738,
#       7.0745,
#       5.5671,
#       2.7097,
#       4.9448,
#       3.7211,
#       2.3608
#     )
#   )

# F0 vs HY
# df <- data.frame(
#     term = c(
#       "System development (GO:0048731)",
#       "Cell migration (GO:0016477)",
#       "Multicellular organism development (GO:0007275)",
#       "Cell adhesion (GO:0007155)",
#       "Circulatory system development (GO:0072359)",
#       "Cell communication (GO:0007154)",
#       "Extracellular matrix organization (GO:0030198)",
#       "Extracellular structure organization (GO:0043062)",
#       "Cell surface receptor signaling pathway (G0:0007166)",
#       "Sodium Ion Transport (GO:0006814)",
#       "Pyruvate Metabolic Process (GO:0006090)",
#       "Lipid Hydroxylation (GO:0002933)",
#       "Organic Acid Metabolic Process (GO:0006082)",
#       "Epithelial Cell Development (GO:0002064)",
#       "Very-Low-Density Lipoprotein Particle Assembly (GO:0034379)"
#     ),
# 
#     value = c(
#      -6.7024,
#      -6.2118,
#      -5.1507,
#      -4.3345,
#      -3.7029,
#      -2.8327,
#      -2.3165,
#      -2.3038,
#      -2.0981,
#      2.0710,
#      1.7319,
#      1.5880,
#      1.4147,
#      0.9858,
#      1.6399
#     )
#   )

# F0 vs NASH
df <- data.frame(
  term = c(
    "Collagen fibril organization (GO:0030199)",
    "Monocarboxylic acid metabolic process (GO:0032787)",
    "Transmembrane transport (GO:0055085)",
    "Chemical homeostasis (GO:0048878)",
    "Lipid catabolic process (GO:0016042)",
    "Regulation of cell adhesion (GO:0030155)",
    "Glucose metabolic process (GO:0006006)",
    "Cellular developmental process (GO:0048869)",
    "cellular response to lipid (GO:0071396)",
    "regulation of fatty acid metabolic process (GO:0019217)"
  ),
  
  value = c(
    2.1899,
    2.3342,
    2.8827,
    3.48548,
    3.6314,
    -5.942,
    -5.4425,
    -4.9381,
    -4.756,
    -2.336
  )
)


# extra proteomica TA  f0 vs hy 
# df <- data.frame(
#   term = c(
#     "Immunoglobulin mediated immune response (GO:0016064)",
#     "Immune system process (GO:0002376)",
#     "Complement activation, classical pathway (GO:0006958)",
#     "Cellular oxidant detoxification (GO:0098869)",
#     "Adaptive immune response (GO:0002250)",
#     "Hydrogen peroxide catabolic process (GO:0042744)",
#     "Negative regulation of fibrinolysis (GO:0051918)",
#     "Lipoprotein metabolic process (GO:0042157)",
#     "Lipid transport (GO:0006869)",
#     "Branched-chain amino acid catabolic process (GO:0009083)",
#     "Fatty acid metabolic process (GO:0006631)",
#     "Lipid metabolic process (GO:0006629)"
#   ),
# 
#   value = c(
#     9.1524,
#     8.4749,
#     8.3979,
#     6.4388,
#     6.3152,
#     5.8632,
#     4.9545,
#     3.4798,
#     3.1129,
#     -8.1817,
#     -6.4271,
#     -5.9281
#   )
# )

# proteomica ta f4 vs hy 
# df <- data.frame(
#   term = c(
#    "Carbon dioxide transport (GO:0015670)",
#    "Oxygen transport (GO:0015671)",
#    "Hydrogen peroxide catabolic process (GO:0041744)",
#    "Cellular oxidant detoxification (GO:0098869)",
#    "Fatty acid metabolic process (GO:0006631)",
#    "Negative regulation of phosphatase activity (GO:0010923)"
#   ),
# 
#   value = c(
#     12.42,
#     11.46,
#     8.02,
#     6.82,
#    -4.42,
#    -3.96
#   )
# )

# proteomica ta+evs f0vshy
# df <- data.frame(
#   term = c(
#     "Adaptive immune response (GO:0002250)",
#     "Cell adhesion (GO:0007155)",
#     "Immune response (GO:0006955)",
#     "Complement activation (G0:0006956)",
#     "Blood coagulation (GO:0007596)",
#     "Inflammatory response (GO:0006954)",
#     "mRNA processing (GO:0006397)",
#     "rRNA processing (GO:0006364)",
#     "RNA splicing (GO:0008380)",
#     "Cellular response to insulin stimulus (GO:0032869)",
#     "Mitotic cell cycle (GO:0000278)"
#   ),
# 
#   value = c(
#   8.1524,
#   7.5017,
#   5.6308,
#   5.4622,
#   5.2526,
#   2.5768,
#   -4.91,
#   -4.91,
#   -4.32,
#   -3.97,
#   -3.80
# 
#   )
# )

# proteomica at+evs f4vshy
# df <- data.frame(
#   term = c(
#     "Adaptive immune response (GO:0002250)",
#     "Cell adhesion (GO:0007155)",
#     "Angiogenesis (GO:0001525)",
#     "Complement activation, classical pathway (G0:0006958)",
#     "Blood coagulation, fibrin clot formation (GO:0072378)",
#     "Acute-phase response (GO:0006953)",
#     "Cellular response to low-density lipoprotein particle (GO:0071404)",
#     "Mitochondrion organization (GO:0007005)",
#     "Negative regulation of translation (GO:0017148)",
#     "Maturation of 5.8S rRNA (GO:0000460)"
#   ),
# 
#   value = c(
#    25.66,
#    8.93,
#    4.16,
#    14.29,
#    5.58,
#    2.89,
#    2.52,
#    -4.39,
#    -4.05,
#    -3.85
#   )
# )

# proteomica ta-evs f4vs f0
# df <- data.frame(
#   term = c(
#     "Cell-matrix adhesion (GO:0007160)",
#     "Integrin-mediated signaling pathway (GO:0007229)",
#     "Small GTPase-mediated signal transduction (GO:0007264)",
#     "Positive regulation of telomerase RNA localization to Cajal body (GO:1904874)",
#     "Positive regulation of protein localization to Cajal body (GO:1904871)",
#     "Positive regulation of telomere maintenance via telomerase (GO:0032212)"
#   ),
# 
#   value = c(
#     4.56,
#     4.44,
#     3.02,
#     -7.61,
#     -7.45,
#     -6.29
#   )
# )



df$group <- ifelse(df$value > 0,"Up","Down")

# Ordena els termes segons value
df$term <- factor(df$term, levels = df$term[order(df$value)])

min_x <- min(df$value)  # punt inicial del gràfic

offset <- 0.1  # marge respecte al 0

ggplot(df, aes(x=value, y=term, fill=group)) +
  geom_col(width=0.7) +
  geom_vline(xintercept=0, size=0.6) +
  
  geom_text(
    aes(
      label = term,
      x = ifelse(group == "Up", 0 - offset, 0 + offset),  # Up: cap a l'esquerra, Down: cap a la dreta
    ),
    hjust = ifelse(df$group == "Up", 1, 0),  
    size = 3  # lletra una mica més petita
  ) +
  
  scale_fill_manual(values=c("Up"="red", "Down"="blue")) +
  labs(x="-log10(adj p-value)", y=NULL, title="F0 vs NASH") +
  theme_classic() +
  theme(
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.line.y = element_blank(),
    legend.position="none"
  ) 