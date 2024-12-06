# -- Archivos --

# Boeck

tmm_Boeck = read.table("..\\Archivos\\Datasets_normalizados\\Boeck_et_al_experiment_gene-count_normalized_TMM.tab", 
                       header = TRUE, 
                       sep = '\t')
tpm_Boeck = read.table("..\\Archivos\\Datasets_normalizados\\\\Boeck_et_al_experiment_gene-count_normalized_TPM.tab", 
                       header = TRUE, 
                       sep = '\t')
tpmclr_Boeck = read.table("..\\Archivos\\Datasets_normalizados\\Boeck_et_al_experiment_gene-count_normalized_TPMCLR.tab", 
                          header = TRUE, 
                          sep = '\t')

# Gomez

tmm_Gomez = read.table("..\\Archivos\\Datasets_normalizados\\gomez-oter_et._al\\gomez_et_al_experiment_gene-count_normalized_TMM.tab", 
                       header = TRUE, 
                       sep = '\t')
tpm_Gomez = read.table("..\\Archivos\\Datasets_normalizados\\gomez_et_al_experiment_gene-count_normalized_TPM.tab", 
                       header = TRUE, 
                       sep = '\t')
tpmclr_Gomez = read.table("..\\Archivos\\Datasets_normalizados\\gomez_et_al_experiment_gene-count_normalized_TPMCLR.tab", 
                          header = TRUE, 
                          sep = '\t')

# Kalesky

tmm_kalesky = read.table("..\\Archivos\\Datasets_normalizados\\kalesky_et_al_experiment_gene-count_normalized_TMM.tab", 
                         header = TRUE, 
                         sep = '\t')
tpm_kalesky = read.table("..\\Archivos\\Datasets_normalizados\\kalesky_et_al_experiment_gene-count_normalized_TPM.tab", 
                         header = TRUE, 
                         sep = '\t')
tpmclr_kalesky = read.table("..\\Archivos\\Datasets_normalizados\\kalesky_et_al_experiment_gene-count_normalized_TPMCLR.tab", 
                            header = TRUE, 
                            sep = '\t')

# Mirza

tmm_mirza = read.table("..\\Archivos\\Datasets_normalizados\\mirza_et_al_experiment_gene-count_normalized_TMM.tab", 
                       header = TRUE, 
                       sep = '\t')
tpm_mirza = read.table("..\\Archivos\\Datasets_normalizados\\mirza_et_al_experiment_gene-count_normalized_TPM.tab", 
                       header = TRUE, 
                       sep = '\t')
tpmclr_mirza = read.table("..\\Archivos\\Datasets_normalizados\\mirza_et_al_experiment_gene-count_normalized_TPMCLR.tab", 
                          header = TRUE, 
                          sep = '\t')


dim(tmm_Boeck)
dim(tmm_Gomez)
dim(tmm_kalesky)
dim(tmm_mirza)

head(tmm_Boeck,2)
head(tmm_Gomez,2)
head(tmm_kalesky,2)
head(tmm_mirza,2)


# Unir las tablas por los rownames
merged_data <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), 
                      list(tmm_Boeck, tmm_Gomez, tmm_kalesky, tmm_mirza))
# Si algún valor está ausente, reemplazar por NA
# merged_data[is.na(merged_data)] <- NaN

cat('Cantidad de genes con datos completos:', sum(complete.cases(merged_data)))

#Para fijarme el gen de Lau
gen_lau = subset(merged_data, merged_data$Gene == 'WBGene00021183')

# -- Para ver si los genes de entrenamiento están en los datasetes --

#Leo los genes
genes = read.csv("..\\Archivos\\Genes_entrenamiento\\Genes_cadena.csv", sep = ';', stringsAsFactors = FALSE)

# Me quedo con los wbid de genes que tengan la palabra "si" en la columna Usar.para.entrenar
genes_entrenamiento <- subset(genes, Usar.para.entrenar == "si")
wbid_entrenamiento <- genes_entrenamiento$WBID

#Me fijo si están todos en los datasets mergeados
genes_lista <- wbid_entrenamiento
genes_no_presentes <- genes_lista[!genes_lista %in% merged_data$Gene]
if (length(genes_no_presentes) > 0) {
  cat('Los siguientes genes no están presentes en merged_data:\n')
  cat(genes_no_presentes, sep = ', ')
} else {
  print('Todos los genes de la lista están presentes en merged_data.')
}

# Me fijo si la expresión de estos genes está completa en el dataset mergeado. 
if (any(is.na(merged_data$WBID[merged_data$WBID %in% wbid_entrenamiento]))) {
  print('Hay al menos un NA en los valores de WBID seleccionados para entrenamiento en merged_data.')
} else {
  print('No hay NA en los valores de WBID seleccionados para entrenamiento en merged_data.')
}

# -- Para hacer distintas combinaciones --

################# TMM ################# 
# Combinar boeck + gomez
boeck_gomez_tmm <- merge(tmm_Boeck, tmm_Gomez, by = "Gene", all = TRUE)

# Combinar boeck + kalesky
boeck_kalesky_tmm <- merge(tmm_Boeck, tmm_kalesky, by = "Gene", all = TRUE)

# Combinar boeck + mirza
boeck_mirza_tmm <- merge(tmm_Boeck, tmm_mirza, by = "Gene", all = TRUE)

# Combinar boeck + gomez + kalesky
boeck_gomez_kalesky_tmm <- merge(boeck_gomez_tmm, tmm_kalesky, by = "Gene", all = TRUE)

# Combinar boeck + gomez + mirza
boeck_gomez_mirza_tmm <- merge(boeck_gomez_tmm, tmm_mirza, by = "Gene", all = TRUE)

# Combinar boeck + kalesky + mirza
boeck_kalesky_mirza_tmm <- merge(boeck_kalesky_tmm, tmm_mirza, by = "Gene", all = TRUE)

# Combinar boeck + gomez + kalesky + mirza
boeck_gomez_kalesky_mirza_tmm <- merge(boeck_gomez_kalesky_tmm, tmm_mirza, by = "Gene", all = TRUE)

# Lista para almacenar los DataFrames
lista_combinaciones_tmm <- list()

# Función para llenar NaNs con ceros
fill_nans <- function(df) {
  df[is.na(df)] <- 0
  return(df)
}

# Combinar boeck + gomez
boeck_gomez_tmm <- fill_nans(boeck_gomez_tmm)
lista_combinaciones_tmm[["boeck_gomez_tmm"]] <- boeck_gomez_tmm

# Combinar boeck + kalesky
boeck_kalesky_tmm <- fill_nans(boeck_kalesky_tmm)
lista_combinaciones_tmm[["boeck_kalesky_tmm"]] <- boeck_kalesky_tmm

# Combinar boeck + mirza
boeck_mirza_tmm <- fill_nans(boeck_mirza_tmm)
lista_combinaciones_tmm[["boeck_mirza_tmm"]] <- boeck_mirza_tmm

# Combinar boeck + gomez + kalesky
boeck_gomez_kalesky_tmm <- fill_nans(boeck_gomez_kalesky_tmm)
lista_combinaciones_tmm[["boeck_gomez_kalesky_tmm"]] <- boeck_gomez_kalesky_tmm

# Combinar boeck + gomez + mirza
boeck_gomez_mirza_tmm <- fill_nans(boeck_gomez_mirza_tmm)
lista_combinaciones_tmm[["boeck_gomez_mirza_tmm"]] <- boeck_gomez_mirza_tmm

# Combinar boeck + kalesky + mirza
boeck_kalesky_mirza_tmm <- fill_nans(boeck_kalesky_mirza_tmm)
lista_combinaciones_tmm[["boeck_kalesky_mirza_tmm"]] <- boeck_kalesky_mirza_tmm

# Combinar boeck + gomez + kalesky + mirza
boeck_gomez_kalesky_mirza_tmm <- fill_nans(boeck_gomez_kalesky_mirza_tmm)
lista_combinaciones_tmm[["boeck_gomez_kalesky_mirza_tmm"]] <- boeck_gomez_kalesky_mirza_tmm


# -- Para guardarlo --
setwd("..\\Archivos\\Datasets_normalizados-combinados\\")

# Escribir la lista_combinaciones en un archivo CSV
lapply(names(lista_combinaciones_tmm), function(name) {
  write.csv(lista_combinaciones_tmm[[name]], file = paste0(name, ".csv"), row.names = FALSE)
})

################# TPM ################# 
# Combinar boeck + gomez
boeck_gomez_tpm <- merge(tpm_Boeck, tpm_Gomez, by = "Gene", all = TRUE)

# Combinar boeck + kalesky
boeck_kalesky_tpm <- merge(tpm_Boeck, tpm_kalesky, by = "Gene", all = TRUE)

# Combinar boeck + mirza
boeck_mirza_tpm <- merge(tpm_Boeck, tpm_mirza, by = "Gene", all = TRUE)

# Combinar boeck + gomez + kalesky
boeck_gomez_kalesky_tpm <- merge(boeck_gomez_tpm, tpm_kalesky, by = "Gene", all = TRUE)

# Combinar boeck + gomez + mirza
boeck_gomez_mirza_tpm <- merge(boeck_gomez_tpm, tpm_mirza, by = "Gene", all = TRUE)

# Combinar boeck + kalesky + mirza
boeck_kalesky_mirza_tpm <- merge(boeck_kalesky_tpm, tpm_mirza, by = "Gene", all = TRUE)

# Combinar boeck + gomez + kalesky + mirza
boeck_gomez_kalesky_mirza_tpm <- merge(boeck_gomez_kalesky_tpm, tpm_mirza, by = "Gene", all = TRUE)

# Lista para almacenar los DataFrames
lista_combinaciones_tpm <- list()

# Función para llenar NaNs con ceros
fill_nans <- function(df) {
  df[is.na(df)] <- 0
  return(df)
}

# Combinar boeck + gomez
boeck_gomez_tpm <- fill_nans(boeck_gomez_tpm)
lista_combinaciones_tpm[["boeck_gomez_tpm"]] <- boeck_gomez_tpm

# Combinar boeck + kalesky
boeck_kalesky_tpm <- fill_nans(boeck_kalesky_tpm)
lista_combinaciones_tpm[["boeck_kalesky_tpm"]] <- boeck_kalesky_tpm

# Combinar boeck + mirza
boeck_mirza_tpm <- fill_nans(boeck_mirza_tpm)
lista_combinaciones_tpm[["boeck_mirza_tpm"]] <- boeck_mirza_tpm

# Combinar boeck + gomez + kalesky
boeck_gomez_kalesky_tpm <- fill_nans(boeck_gomez_kalesky_tpm)
lista_combinaciones_tpm[["boeck_gomez_kalesky_tpm"]] <- boeck_gomez_kalesky_tpm

# Combinar boeck + gomez + mirza
boeck_gomez_mirza_tpm <- fill_nans(boeck_gomez_mirza_tpm)
lista_combinaciones_tpm[["boeck_gomez_mirza_tpm"]] <- boeck_gomez_mirza_tpm

# Combinar boeck + kalesky + mirza
boeck_kalesky_mirza_tpm <- fill_nans(boeck_kalesky_mirza_tpm)
lista_combinaciones_tpm[["boeck_kalesky_mirza_tpm"]] <- boeck_kalesky_mirza_tpm

# Combinar boeck + gomez + kalesky + mirza
boeck_gomez_kalesky_mirza_tpm <- fill_nans(boeck_gomez_kalesky_mirza_tpm)
lista_combinaciones_tpm[["boeck_gomez_kalesky_mirza_tpm"]] <- boeck_gomez_kalesky_mirza_tpm


# -- Para guardarlo --

# Escribir la lista_combinaciones en un archivo CSV
lapply(names(lista_combinaciones_tpm), function(name) {
  write.csv(lista_combinaciones_tpm[[name]], file = paste0(name, ".csv"), row.names = FALSE)
})

################# TPM + CLR ################# 
# Combinar boeck + gomez
boeck_gomez_tpmclr <- merge(tpmclr_Boeck, tpmclr_Gomez, by = "Gene", all = TRUE)

# Combinar boeck + kalesky
boeck_kalesky_tpmclr <- merge(tpmclr_Boeck, tpmclr_kalesky, by = "Gene", all = TRUE)

# Combinar boeck + mirza
boeck_mirza_tpmclr <- merge(tpmclr_Boeck, tpmclr_mirza, by = "Gene", all = TRUE)

# Combinar boeck + gomez + kalesky
boeck_gomez_kalesky_tpmclr <- merge(boeck_gomez_tpmclr, tpmclr_kalesky, by = "Gene", all = TRUE)

# Combinar boeck + gomez + mirza
boeck_gomez_mirza_tpmclr <- merge(boeck_gomez_tpmclr, tpmclr_mirza, by = "Gene", all = TRUE)

# Combinar boeck + kalesky + mirza
boeck_kalesky_mirza_tpmclr <- merge(boeck_kalesky_tpmclr, tpmclr_mirza, by = "Gene", all = TRUE)

# Combinar boeck + gomez + kalesky + mirza
boeck_gomez_kalesky_mirza_tpmclr <- merge(boeck_gomez_kalesky_tpmclr, tpmclr_mirza, by = "Gene", all = TRUE)

# Lista para almacenar los DataFrames
lista_combinaciones_tpmclr <- list()

# Combinar boeck + gomez
boeck_gomez_tpmclr <- fill_nans(boeck_gomez_tpmclr)
lista_combinaciones_tpmclr[["boeck_gomez_tpmclr"]] <- boeck_gomez_tpmclr

# Combinar boeck + kalesky
boeck_kalesky_tpmclr <- fill_nans(boeck_kalesky_tpmclr)
lista_combinaciones_tpmclr[["boeck_kalesky_tpmclr"]] <- boeck_kalesky_tpmclr

# Combinar boeck + mirza
boeck_mirza_tpmclr <- fill_nans(boeck_mirza_tpmclr)
lista_combinaciones_tpmclr[["boeck_mirza_tpmclr"]] <- boeck_mirza_tpmclr

# Combinar boeck + gomez + kalesky
boeck_gomez_kalesky_tpmclr <- fill_nans(boeck_gomez_kalesky_tpmclr)
lista_combinaciones_tpmclr[["boeck_gomez_kalesky_tpmclr"]] <- boeck_gomez_kalesky_tpmclr

# Combinar boeck + gomez + mirza
boeck_gomez_mirza_tpmclr <- fill_nans(boeck_gomez_mirza_tpmclr)
lista_combinaciones_tpmclr[["boeck_gomez_mirza_tpmclr"]] <- boeck_gomez_mirza_tpmclr

# Combinar boeck + kalesky + mirza
boeck_kalesky_mirza_tpmclr <- fill_nans(boeck_kalesky_mirza_tpmclr)
lista_combinaciones_tpmclr[["boeck_kalesky_mirza_tpmclr"]] <- boeck_kalesky_mirza_tpmclr

# Combinar boeck + gomez + kalesky + mirza
boeck_gomez_kalesky_mirza_tpmclr <- fill_nans(boeck_gomez_kalesky_mirza_tpmclr)
lista_combinaciones_tpmclr[["boeck_gomez_kalesky_mirza_tpmclr"]] <- boeck_gomez_kalesky_mirza_tpmclr


# -- Para guardarlo --

# Escribir la lista_combinaciones en un archivo CSV
lapply(names(lista_combinaciones_tpmclr), function(name) {
  write.csv(lista_combinaciones_tpmclr[[name]], file = paste0(name, ".csv"), row.names = FALSE)
})

