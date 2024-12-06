library(edgeR)

# Dado que los datos de expresión de Boeck et al. estaban en formato dcpm, 
# hay que convertirlos a TMM, TPM y CLRTPM
expression_table <- read.csv("..\\Archivos\\Tablas_conteo\\expresion_Boeck.csv", 
                             header = TRUE, 
                             row.names = 2, 
                             sep=';')

# Sacamos las dos primeras columans
expression_table <- expression_table[, -(1:2)]


d <- DGEList(counts = expression_table)

############## TMM #################

d <- calcNormFactors(d, method = "TMM")
normalized_expression <- cpm(d, normalized.lib.sizes = TRUE)
write.table(normalized_expression, file="..\\Archivos\\Datasets-normalizados\\Boeck_et_al_experiment_gene-count_normalized_TMM.tab", 
            quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)


############## TPM #################

norm_factors <- d$samples$norm.factors

tpm <- cpm(d, normalized.lib.sizes = TRUE) * 1e6 / colSums(cpm(d) * norm_factors)

write.table(tpm, file="..\\Archivos\\Datasets-normalizados\\Boeck_et_al_experiment_gene-count_normalized_TPM.tab", 
            quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

############## TPM + CLR #################

library(compositions)

# Apply CLR transformation to TPM-normalized data
clr_transformed_data <- clr(tpm)

write.table(clr_transformed_data, file="..\\Archivos\\Datasets-normalizados\\Boeck_et_al_experiment_gene-count_normalized_TPMCLR.tab", 
            quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

