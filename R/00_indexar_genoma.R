# Antes de hacer las tablas de conteo, es necesario indexar el genoma.

library(Rsubread)

# Set the path to the directory containing the genome file
setwd('..\\Archivos\\tablas_conteo\\Genoma_indexado\\')
# Run index with basename "genome_celegansWBcel235_rsubread"
buildindex(basename = "genome_celegansWBcel235_rsubread", 
           reference = paste0(genome_directory, "GCF_000002985.6_WBcel235_genomic.fna"))

# Nota: los archivos resultado de este script no se encuentran en github ya que 
# cerca de 1.86Gb