# datos
# https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&DataMode=FactorsData&StudyID=ST001322&StudyType=MS&ResultType=1#DataTabs
# ST001322



library(metabolomicsWorkbenchR)

t(do_query('study','study_id','ST001322','summary'))

### descargar datos

SE <- do_query(
      context = 'study',
      input_item = 'analysis_id',
      input_value = 'AN002198',
      output_item = 'SummarizedExperiment'
    )


### como la estructura del query es diferente a la de los ejemplos
### procedo a extraer los datos manuales y volver a crear un SummarizedExperiment
library(SummarizedExperiment)

# Extraer la matriz de intensidades
matriz <- as.data.frame(SE@assays@data@listData[[1]])

# Extraer rowData y usar los nombres reales de metabolitos
rowData <- as.data.frame(SE@elementMetadata)
rownames(matriz) <- rowData$metabolite          # poner los nombres en la matriz
rownames(rowData) <- rowData$metabolite         # poner los nombres en rowData

# Extraer colData y alinear nombres
colData <- as.data.frame(SE@colData)
rownames(colData) <- colnames(matriz)           # los nombres de las muestras

# Crear el objeto SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(abundancias = matriz),
  rowData = rowData,
  colData = colData
)

assay(se)

# Eliminar metabolitos con nombre NA
tail(rownames(se), 10)
se <- se[!rownames(se) %in% c("NA", "NA.1"), ]


# Convertir toda la matriz a números
matriz <- assay(se)






head(colData(se))

str(colData(se))
summary(colData(se))

unique(se$treatment)
unique(se$group)
table(se$treatment)
table(se$group)


df <- as.data.frame(colData(se))
df$tiempo <- gsub(".* (2 week|2 month|6 month).*", "\\1", df$treatment)

colData(se)$tiempo <- df$tiempo



library(ggplot2)
df <- as.data.frame(colData(se))

ggplot(df, aes(x = tiempo, fill = group)) +
  geom_bar(position = "dodge") +
  labs(title = "Distribución de muestras por grupo y tiempo",
       x = "Tiempo", y = "Número de muestras") +
  theme_minimal()


library(ggplot2)
library(reshape2)

# Obtener matriz de abundancias
matriz <- assay(se)

# Convertir a formato largo para ggplot
df_long <- melt(matriz)
colnames(df_long) <- c("Metabolito", "Muestra", "Abundancia")

# Gráfico
ggplot(df_long, aes(x = Muestra, y = Abundancia)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Distribución de abundancias por muestra",
       x = "Muestras", y = "Abundancia (intensidad)")


library(ggplot2)
library(ggfortify)

# Transponer para tener muestras en filas

# 1. Extraer la matriz como data.frame para asegurar control total
matriz_df <- as.data.frame(assay(se), stringsAsFactors = FALSE)

# 2. Convertir cada columna (muestra) a numérico explícitamente
matriz_df <- data.frame(lapply(matriz_df, function(x) as.numeric(as.character(x))),
                        check.names = FALSE, row.names = rownames(se))

# 3. Convertir de nuevo a matriz
matriz_numerica <- as.matrix(matriz_df)

# 4. Confirmar que todo es numérico
str(matriz_numerica)
matriz_t <- t(matriz_numerica)

# Confirmar que sigue siendo numérica
str(matriz_t)



# PCA
pca <- prcomp(matriz_t, scale. = TRUE)

# Agregar metadatos
df_pca <- as.data.frame(pca$x)
df_pca$group <- colData(se)$group
df_pca$tiempo <- colData(se)$tiempo

# Gráfico PCA
ggplot(df_pca, aes(x = PC1, y = PC2, color = group, shape = tiempo)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA de abundancias de metabolitos",
       x = "PC1", y = "PC2")

