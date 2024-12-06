# Instalar paquetes
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(version = "3.18")
BiocManager::install("Rsubread")
BiocManager::install("edgeR", quietly=TRUE, ask=FALSE, update=FALSE)
BiocManager::install("biomaRt", quietly=TRUE, ask=FALSE, update=FALSE)
BiocManager::install("EnhancedVolcano", quietly=TRUE, ask=FALSE, update=FALSE)
BiocManager::install("GO.db", quietly=TRUE, ask=FALSE, update=FALSE)
BiocManager::install("org.Mm.eg.db", quietly=TRUE, ask=FALSE, update=FALSE)
#BiocManager::install("ClusterProfiler", quietly=TRUE, ask=FALSE, update=FALSE)

install.packages("gplots", quietly=TRUE)
install.packages("pheatmap", quietly=TRUE)
install.packages("fansi", quietly=TRUE)
install.packages("utf8", quietly=TRUE)

## Important function
## https://stackoverflow.com/questions/70025153/how-to-access-the-shell-in-google-colab-when-running-the-r-kernel
shell_call <- function(command, ...) {
  result <- system(command, intern = TRUE, ...)
  cat(paste0(result, collapse = "\n"))
}


# Importar paquetes
library(Rsubread)
library(ShortRead)
library(Rsamtools)
library(edgeR)
library(biomaRt)
library(ggplot2)
library(EnhancedVolcano)
library(gplots)
library(pheatmap)

# ----- 1. Mapeo de los reads al genoma -----

# NOTA: los reads no se encuentran subido a github debido al tamaño de los 
# archivos, pero se pueden desacargar de NCBI

#list fastq files of interest
path_reads = 'path\\to\\reads\\'

my.fastqs<- list.files(path= path_reads,
                       pattern=".+\\.gz$", 
                       full.names=TRUE)

# List for _pass_1.fastq.gz files
my.fastqs_1 <- my.fastqs[grep("pass_1.fastq.gz$", my.fastqs)]

# List for _pass_2.fastq.gz files
my.fastqs_2 <- my.fastqs[grep("pass_2.fastq.gz$", my.fastqs)]

# Alineamiento
align(index="..\\Archivos\\Tablas_conteo\\Genoma_indexado\\", 
      readfile1=my.fastqs_1,
      readfile2=my.fastqs_2,
      input_format="gzFASTQ", 
      nthreads=6)

#obtain filenames for bam files
my.bams <- list.files(path=path_reads, 
                      pattern="+\\.BAM$", 
                      full.names=TRUE)

# ----- 2. Conteo de reads -----
fc <- featureCounts(my.bams, 
                    annot.ext= "..\\Archivos\\Tablas_conteo\\Genoma_indexado\\Genomic.gtf", 
                    isGTFAnnotationFile=TRUE, 
                    nthreads=6,
                    GTF.featureType= 'gene',
                    countMultiMappingReads=TRUE,
                    isPairedEnd=TRUE	
)


#Tengo que sacarle el "GeneID:" de las filas para después poder usar esos identificadores
head(fc$counts,2)


# Para explorar la variable que se genera
class(fc)
names(fc)
head(fc$counts, 2)
head(fc$stat,3)
colnames(fc$stat)

# Para ponerle nombre a las columnas
my.samples <- c("Ecoli_15_.1",
                "Ecoli_15_.2",
                "Ecoli_15_.3",
                "Ecoli_20_.1",
                "Ecoli_20_.2",
                "Ecoli_20_.3",
                "Ecoli_25_.1",
                "Ecoli_25_.2",
                "Ecoli_25_.3",
                "Bsub_15_.1",
                "Bsub_15_.2",
                "Bsub_15_.3",
                "Bsub_20_.1",
                "Bsub_20_.2",
                "Bsub_20_.3",
                "Bsub_25_.1",
                "Bsub_25_.2",
                "Bsub_25_.3"
                )

colnames(fc$counts) <- my.samples

rownames(fc$counts) <- sub("^CELE_", "", rownames(fc$counts))


setwd('..\\Archivos\\Datasets_normalizados\\')
#save to an external file
write.table(fc$counts, file="Gomez-oter_experiment_gene-count.tab", 
            quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

#save to an RDS format object
saveRDS(fc, file="Gomez-oter_experiment_fc.Rds")

# ----- 3. Normalizacion -----

path_samples = '..\\Archivos\\Tablas_conteo\\Gomez-oter_info_samples.csv'
targets <- read.csv(path_samples, 
                    header=TRUE, 
                    sep = ";")

rownames(targets) <- targets$Samples

# ----- 3.1. Estandarizar nombres  -----
# Para ver los datasets disponibles https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html#step-2-choosing-a-dataset

listDatasets(useEnsembl(biomart = "genes"))

#selecting an Ensembl BioMart database and dataset
mart <- useMart("parasite_mart", 
                dataset = "wbps_gene", 
                host = "https://parasite.wormbase.org", 
                port = 443)

#listDatasets(useMart(biomart = "parasite_mart", 
#                     host = "https://parasite.wormbase.org", 
#                     port = 443))
#datasets <- listDatasets(mart)

genes.df <- getBM(attributes = c('wbps_gene_id',
                                 'external_gene_id',
                                 'description',
                                 'gene_biotype',
                                 'transcript_count', 
                                 'wormbase_gseq'),
                  filters = 'wormbase_gseqname', #este es el id que se usa en featurecounts
                  values = rownames(fc$counts),
                  mart = mart)
type(genes.df)
head(genes.df,3)
dim(genes.df)
dim(fc$counts)

#########################################
# NOTA: Acá tengo un problema y es que fc$counts es más grande que genes.df, 
# seguramente porque hay códigos que no encontró. Me voy a sacar de arriba los 
# genes de fc$count que no están en genes.df: 

# Crear un vector de entrezgene_id desde genes.df
wb_seqname <- genes.df$wormbase_gseq

# Filtrar las filas de fc$count basado en entrezgene_id
fc$counts_filtrado <- fc$counts[row.names(fc$counts) %in% wb_seqname, ]

rownames(genes.df) <- genes.df$wormbase_gseq

wb_seqname <- rownames(fc$counts_filtrado)
genes.df = genes.df[row.names(genes.df)%in% wb_seqname, ]
head(genes.df,2)

dim(genes.df)
dim(fc$counts_filtrado)


#########################################

# Ahora hago la DGEList
y <- DGEList(counts=fc$counts_filtrado, samples=targets, genes=genes.df, group=targets$Genotype)


# Para cambiar el nombre de las columnas por el WBID

# Extract the "counts_normalized" matrix and the relevant columns from the "genes" list
counts <- y$counts
genes_subset <- y$genes[, c("wbps_gene_id", "wormbase_gseq")]

# Match the row names in "counts_normalized" with "entrezgene_id" and replace with "ensembl_gene_id"
rownames(counts) <- genes_subset$wbps_gene_id[match(rownames(counts), 
                                                    genes_subset$wormbase_gseq)]
head(counts,2)
dim(counts)
sum(row.names(counts) %in% c("", NULL, NA))

# Update the "counts_normalized" list in the "y" object
y$counts <- counts
head(y$counts,2)


class(y)
dim(y)

# ----- 3.2. Filtrar conteos bajos  -----

10/min(y$samples$lib.size)*1000000

table(y$sample$group)
keep <- rowSums(cpm(y) > 0) >= 3
table(keep)

y <- y[keep, , keep.lib.sizes=FALSE]
dim(y)

# ----- 3.3. Normalizar por composition bias  -----


# Normalization by trimmed mean of M values (TMM) is performed by using the 
# `calcNormFactors` function, which returns the DGEList argument with only the 
# `norm.factors` changed. It calculates a set of normalization factors, one for 
# each sample, to eliminate composition biases between libraries. The product of 
# these factors and the library sizes defines the effective library size, which 
# replaces the original library size in all downstream analyses.

y <- calcNormFactors(y)
y$samples

## Visualización: Normalization for composition bias

# Chequeamos por tipo de alimentación 
group.c <- as.factor(y$samples$Groups)
group.c
# group.c <- relevel(group.c, ref="Ecoli_15")

pch <- c(15:23)
colors <- rep(c("blue", "red", "green", "purple", "orange", "pink"), each = 3)

plotMDS(y, col=colors, pch=16, cex=2)
legend("bottomleft", legend=levels(group.c), pch=16, col=unique(colors), ncol=2)

# ----- 3.4. Normalizar por TMM  -----

y <- calcNormFactors(y, method = "TMM")

## NOTA: calcNormFactors() crea una columna en y$samples con el factor de norm
# de cada sample. Cuando después se aplica cpm() sobre el objeto, toma en cuenta
# ese norm factor, entonces el resultado es en TMM https://www.biostars.org/p/317701/

# Extract normalized counts
y$counts_normalized_tmm <- cpm(y)

# View the normalized counts
head(y$counts_normalized_tmm)

## Para calcular los promedios. 

# Extract relevant data from the "samples" table
groups <- unique(y$samples$Groups)
counts_normalized <- y$counts_normalized_tmm

sum(row.names(counts_normalized) %in% c("", NULL, NA))

# Initialize an empty matrix for the averages
y$counts_norm_tmm_prom <- matrix(0, nrow = nrow(y$counts_normalized_tmm), ncol = length(groups))

# Assign column names to the new matrix
colnames(y$counts_norm_tmm_prom) <- groups
rownames(y$counts_norm_tmm_prom) <- rownames(y$counts_normalized_tmm)

head(y$counts_norm_tmm_prom)

# Calculate and fill in the averages for each group
for (i in seq_along(groups)) {
  group <- groups[i]
  group_cols <- grep(group, colnames(counts_normalized))
  y$counts_norm_tmm_prom[, i] <- rowMeans(counts_normalized[, group_cols], na.rm = TRUE)
}



head(y$counts_norm_tmm_prom)
dim(y$counts_norm_tmm_prom)

#save to an external file
write.table(y$counts_norm_tmm_prom, file="..\\Archivos\\Datasets_normalizados\\gomez_et_al_experiment_gene-count_normalized_TMM.tab", 
            quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)


# ----- 3.5. Normalizar por TPM  -----

# Para pasar a TPM: 
# Access the normalization factors
norm_factors <- y$samples$norm.factors

# Apply TPM normalization
tpm <- cpm(y, normalized.lib.sizes = TRUE) * 1e6 / colSums(cpm(y) * norm_factors)

## Para calcular los promedios. 

# Extract relevant data from the "samples" table
#groups <- unique(y$samples$Groups)
counts_normalized <- tpm

sum(row.names(counts_normalized) %in% c("", NULL, NA))

# Initialize an empty matrix for the averages
y$counts_norm_prom_tpm <- matrix(0, nrow = nrow(counts_normalized), ncol = length(groups))
rownames(y$counts_norm_prom_tpm) <- rownames(counts_normalized)
colnames(y$counts_norm_prom_tpm) <- groups

head(y$counts_norm_prom_tpm)


# Calculate and fill in the averages for each group
for (i in seq_along(groups)) {
  group <- groups[i]
  group_cols <- grep(group, colnames(counts_normalized))
  y$counts_norm_prom_tpm[, i] <- rowMeans(counts_normalized[, group_cols], na.rm = TRUE)
}

head(y$counts_norm_prom_tpm)
dim(y$counts_norm_prom_tpm)

write.table(y$counts_norm_prom_tpm, file="..\\Archivos\\Datasets_normalizados\\gomez_et_al_experiment_gene-count_normalized_TPM.tab", 
            quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

# ----- 3.6. Normalizar por TPM + CLR  -----

library(compositions)

# Apply CLR transformation to TPM-normalized data
clr_transformed_data <- clr(tpm)
y$counts_norm_tpm_clr = clr_transformed_data

head(y$counts_norm_tpm_clr) 
head(y$counts_normalized)

## Para calcular los promedios. 

# Extract relevant data from the "samples" table
#groups <- unique(y$samples$tissue)
counts_normalized <- y$counts_norm_tpm_clr

sum(row.names(counts_normalized) %in% c("", NULL, NA))

# Initialize an empty matrix for the averages
y$counts_norm_prom_tpm_clr <- matrix(0, nrow = nrow(counts_normalized), ncol = length(groups))
rownames(y$counts_norm_prom_tpm_clr) <- rownames(counts_normalized)
head(y$counts_norm_prom_tpm_clr)

# Calculate and fill in the averages for each group
for (i in seq_along(groups)) {
  group <- groups[i]
  group_cols <- grep(group, colnames(counts_normalized))
  y$counts_norm_prom_tpm_clr[, i] <- rowMeans(counts_normalized[, group_cols], na.rm = TRUE)
}

# Assign column names to the new matrix
colnames(y$counts_norm_prom_tpm_clr) <- groups

head(y$counts_norm_prom_tpm_clr)
dim(y$counts_norm_prom_tpm_clr)

#save to an external file
write.table(y$counts_norm_prom_tpm_clr, file="..\\Archivos\\Datasets_normalizados\\gomez_et_al_experiment_gene-count_normalized_TPMCLR.tab", 
            quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

#save to an RDS format object
saveRDS(y, file="..\\Archivos\\Tablas_conteo\\gomez_et_al_experiment_y.Rds")

