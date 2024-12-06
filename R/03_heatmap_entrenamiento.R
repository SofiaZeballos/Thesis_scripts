library(ggplot2)
library(dplyr)
library(reshape2)
library(pheatmap)
library(gplots)

# Leer el archivo CSV con los genes
genes <- read.csv("..\\Archivos\\Genes_entrenamiento\\Genes_cadena.csv")

# Filtrar los genes a usar basados en la columna "Usar.para.entrenar"
wbid_vector <- genes$WBID[genes$Usar.para.entrenar == 'si']

# ------- Heatmap para el dataset de sc ------- 
set.seed(0)
data <- read.csv("C:\\Users\\masoz\\Desktop\\TesisMaestria\\Grafo\\geneExpr_sinsacar0s.csv", sep = ',')
# Asignar la primera columna como nombres de fila
rownames(data) <- data[, 1]

# Eliminar la primera columna del data frame
data <- data[, -1]

# Guardar los nombres de fila originales
rownames_originales <- rownames(data)

# Normalizar los datos
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

data_normalized <- as.data.frame(lapply(data, normalize))
rownames(data_normalized) <- rownames_originales

# Seleccionar las columnas que coincidan con wbid_vector
selected_data <- data_normalized[, colnames(data_normalized) %in% wbid_vector]

# Crear un vector de correspondencia entre WBID y Gen.en.Celegans
wbid_to_gene <- setNames(genes$Gen.en.Celegans, genes$WBID)

# Renombrar las columnas de selected_data usando la correspondencia
colnames(selected_data) <- wbid_to_gene[colnames(selected_data)]

# Seleccionar aleatoriamente 10 genes (columnas)
random_genes <- sample(colnames(data_normalized), 10)

# Extraer los datos de estos genes seleccionados
random_data <- data_normalized[, random_genes]

# Convertir el data frame de genes aleatorios al formato largo
random_data_long <- melt(as.data.frame(random_data), variable.name = "Gene", value.name = "Expression")

# Asignar un nombre genérico para las células
random_data_long$Cell <- rownames(data_normalized)

# Asignar un Complejo genérico "Random" a estos genes
random_data_long$Complejo <- "R"

# Convertir el data frame original a formato largo
selected_data$Cell <- rownames(selected_data)
selected_data_long <- melt(selected_data, id.vars = "Cell", variable.name = "Gene", value.name = "Expression")

# Crear un vector de correspondencia entre los genes seleccionados y sus complejos
gene_complex_mapping <- setNames(genes$Complejo, genes$Gen.en.Celegans)

# Asignar el complejo correspondiente a cada gen en el formato largo
selected_data_long$Complejo <- gene_complex_mapping[as.character(selected_data_long$Gene)]

# Combinar los datos originales con los datos aleatorios
combined_data_long <- rbind(selected_data_long, random_data_long)

# Ordenar los datos para que "R" aparezca al final
combined_data_long$Complejo <- factor(combined_data_long$Complejo, 
                                      levels = c(sort(unique(selected_data_long$Complejo)), "R"))

# Crear el heatmap con el nuevo bloque de genes aleatorios
ggplot(combined_data_long, aes(x = Cell, y = Gene, fill = Expression)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", name = "Expresión\nnormalizada")+
  theme_minimal() +
  labs(title = "", x = "Células", y = "Genes") +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.5, "cm"), 
        legend.text = element_text(size = 6),     # Tamaño del texto de la leyenda
        legend.title = element_text(size = 7),    # Tamaño del título de la leyenda
        axis.text.y = element_text(size = 7),     # Tamaño del texto del eje y
        axis.title.x = element_text(size = 8),    # Tamaño del título del eje x
        axis.title.y = element_text(size = 8)) +  # Tamaño del título del eje y
  facet_grid(Complejo ~ ., scales = "free_y", space = "free_y")
