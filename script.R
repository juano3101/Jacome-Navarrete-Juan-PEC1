# datos
# https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&DataMode=FactorsData&StudyID=ST001322&StudyType=MS&ResultType=1#DataTabs
# ST001322

# datos
# https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&DataMode=FactorsData&StudyID=ST001322&StudyType=MS&ResultType=1#DataTabs
# ST001322

#ion positive mod

library(metabolomicsWorkbenchR)
t(do_query('study','study_id','ST000957','summary'))

### descargar datos
library(SummarizedExperiment)
library(readr)      # para leer TSV
library(tibble)     # para asegurar nombres de filas
library(S4Vectors)  # para DataFrame


# Cambia los nombres de archivo si es necesario
assay <- read.delim("datos/assay.tsv", row.names = 1, check.names = FALSE)
coldata <- read.delim("datos/coldata.txt", row.names = 1, check.names = FALSE)
rowdata <- read.delim("datos/rowdata.tsv", row.names = 1, check.names = FALSE)

# Reordenar coldata para que coincida con el orden de las columnas de assay
coldata <- coldata[colnames(assay), , drop = FALSE]

se <- SummarizedExperiment(
  assays = list(count = assay),
  colData = coldata,
  rowData = rowdata
)