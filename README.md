Aquest repositori conté els scripts utilitzats per a l’anàlisi de dades transcriptòmiques (RNA-seq) i proteòmiques en el context de l’estudi de la malaltia hepàtica associada a disfunció metabòlica (MASLD), amb especial focus en la caracterització de vesícules extracel·lulars (EVs) i teixit adipós (TA).

Estructura del repositori

Els scripts es poden agrupar en dues línies principals d’anàlisi:

1. Anàlisi de RNA-seq
Aquesta part inclou els scripts destinats a l’anàlisi de dades transcriptòmiques:

obtencióDEGs.R
Script utilitzat per a la identificació de gens diferencialment expressats (DEGs) entre diferents estadis de la malaltia.
Grafic-Barres-GOs.R
Script destinat a la visualització dels resultats d’enriquiment funcional mitjançant Gene Ontology (GO) aconseguits amb Enrichr, representant els termes més significatius en forma de gràfics de barres.

2. Anàlisi de proteòmica (TA i EVs)
Aquesta part inclou els scripts utilitzats per a l’anàlisi de dades proteòmiques, tant per a teixit adipós (TA) com per a vesícules extracel·lulars (EVs).

Caracteritzacio-Proteines.R
Script utilitzat per a la caracterització de les dades proteòmiques, incloent l’exploració de patrons d’expressió i anàlisi descriptiva.
proteomica TA+Evs.R
Script que analitza les dades proteòmiques de TA i EVs. Aquest script segueix una metodologia equivalent a la utilitzada en RNA-seq, adaptant únicament el conjunt de dades d’entrada, en aquest cas es troben els DEPs i després GOs.

3. Anàlisi exploratòria i reducció de dimensionalitat
pca.R
Script per a la realització d’anàlisi de components principals (PCA). Aquest script és adaptable a diferents conjunts de dades (RNA-seq, proteòmica de TA i proteòmica d’EVs), ja que la seva estructura permet modificar fàcilment el fitxer d’entrada.
