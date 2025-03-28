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
# Añadir columna de sujeto al colData
coldata$sujeto <- gsub("_.*", "", rownames(coldata))

se <- SummarizedExperiment(
  assays = list(count = assay),
  colData = coldata,
  rowData = rowdata
)

# Seleccionar solo las columnas deseadas
# Conservar solo las columnas que deseas
rowData(se) <- rowData(se)[, c("Standardized_name", "Formula", "Super_class", "Main_class", "Sub_class")]

dim(se)

## revisión del objeto

# Ver si hay NAs
sum(is.na(assay(se)))

# Obtener vector lógico de metabolitos con suma 0 de abundancia
metabolitos_cero <- rowSums(assay(se)) == 0
# Cuántos son
sum(metabolitos_cero)
# Ver nombres de esos metabolitos
rownames(se)[metabolitos_cero]
# quitar metabolitos que tienen 0
se <- se[!metabolitos_cero, ]

# Crear vector lógico con metabolitos que no estan dentro de la base de datos
# de Metabolomics Workbench
sin_nombre <- rowData(se)$Standardized_name == "-"
# Filtrar solo los que tienen nombre válido
se <- se[!sin_nombre, ]



dim(se)


se
dim(se)
assayNames(se)
head(rownames(se))
head(colnames(se))
colData(se)
rowData(se)


table(colData(se)$Milk_fraction)



library(dplyr)
library(tibble)
library(ggplot2)


# Extraer info y preparar
matriz <- assay(se)
metadata <- as.data.frame(colData(se)) %>%
  rownames_to_column("Muestra")

# Calcular metabolitos detectados (abundancia > 0) por muestra
# Calcular número de metabolitos detectados por muestra (2 es para cada col.) 
metabolitos_por_muestra <- apply(matriz, 2, function(x) sum(x > 0))

# Agregar al metadata
metadata$metabolitos_detectados <- metabolitos_por_muestra

# Agrupar por sujeto y fracción, y graficar
ggplot(metadata, aes(x = Milk_fraction, y = metabolitos_detectados, group = sujeto, color = sujeto)) +
  geom_line(size = 1.1) +
  geom_point(size = 2) +
  labs(title = "Metabolitos detectados por fracción de leche y sujeto",
       x = "Fracción de leche", y = "Número de metabolitos detectados") +
  theme_minimal()




library(dplyr)
library(ggplot2)


# Calcular media, SD y SEM por fracción
resumen_fraccion <- metadata %>%
  group_by(Milk_fraction) %>%
  summarise(
    media = mean(metabolitos_detectados),
    sd = sd(metabolitos_detectados),
    n = n(),
    sem = sd / sqrt(n)
  )

# Calcular media, SD y SEM por fracción
resumen_sujeto <- metadata %>%
  group_by(sujeto) %>%
  summarise(
    media = mean(metabolitos_detectados),
    sd = sd(metabolitos_detectados),
    n = n(),
    sem = sd / sqrt(n)
  )

resumen_combinado <- metadata %>%
  group_by(sujeto, Milk_fraction) %>%
  summarise(
    media = mean(metabolitos_detectados),
    n = n(),
    .groups = "drop"
  )







library(dplyr)
library(tibble)
library(ggplot2)
# Calcular suma de intensidades por muestra
suma_por_muestra <- colSums(assay(se))

metadata$intensidad_por_muestra <- suma_por_muestra


# Gráfico: barras por sujeto dentro de cada fracción
ggplot(metadata, aes(x = sujeto, y = intensidad_por_muestra, fill = Milk_fraction)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Intensidad total por sujeto y fracción de leche",
       x = "Sujeto", y = "Intensidad total") +
  theme_minimal()






library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)

# Convertir la matriz a data frame y pasar a formato largo
matriz_df <- as.data.frame(assay(se)) %>%
  rownames_to_column(var = "metabolito") %>%
  pivot_longer(-metabolito, names_to = "muestra", values_to = "intensidad")

# Añadir la fracción de leche desde colData
matriz_df <- matriz_df %>%
  left_join(as.data.frame(colData(se)) %>%
              rownames_to_column("muestra"),
            by = "muestra")

ggplot(matriz_df, aes(x = Milk_fraction, y = intensidad, fill = sujeto)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  scale_y_log10() +
  facet_wrap(~ sujeto) +
  labs(title = "Distribución log10 de intensidades por fracción y sujeto",
       x = "Fracción de leche", y = "Log10 Intensidad") +
  theme_minimal()


