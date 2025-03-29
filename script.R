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

colData(se)[["metabolitos_detectados"]] <- metabolitos_por_muestra

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
colData(se)[["intensidad_por_muestra"]] <- suma_por_muestra

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
  pivot_longer(-metabolito, names_to = "Muestra", values_to = "intensidad")

# Añadir la fracción de leche desde colData
matriz_df <- matriz_df %>%
  left_join(as.data.frame(colData(se)) %>%
              rownames_to_column("Muestra"),
            by = "Muestra")

ggplot(matriz_df, aes(x = Milk_fraction, y = intensidad, fill = sujeto)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  scale_y_log10() +
  facet_wrap(~ sujeto) +
  labs(title = "Distribución log10 de intensidades por fracción y sujeto",
       x = "Fracción de leche", y = "Log10 Intensidad") +
  theme_minimal()






library(dplyr)
library(ggplot2)


# Abundancia absoluta promedio por fracción
matriz_promedio <- matriz_df %>%
  group_by(Milk_fraction, metabolito) %>%
  summarise(media = mean(intensidad), .groups = "drop")

ggplot(matriz_promedio, aes(x = Milk_fraction, y = media, fill = metabolito)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundancia absoluta promedio por fracción de leche",
       x = "Fracción de leche", y = "Media de intensidad") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


# calcular el porcentaje relativo por fracción
matriz_relativa_frac <- matriz_promedio %>%
  group_by(Milk_fraction) %>%
  mutate(porcentaje = media / sum(media) * 100) %>%
  ungroup()

# 3. Graficar
ggplot(matriz_relativa_frac, aes(x = Milk_fraction, y = porcentaje, fill = metabolito)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundancia relativa promedio por fracción de leche",
       x = "Fracción de leche", y = "Porcentaje (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")






# 1. Calcular la media por sujeto, fracción y metabolito
matriz_promedio_sujeto <- matriz_df %>%
  group_by(sujeto, Milk_fraction, metabolito) %>%
  summarise(media = mean(intensidad), .groups = "drop")

# 2. Graficar barras apiladas con facet por sujeto
ggplot(matriz_promedio_sujeto, aes(x = Milk_fraction, y = media, fill = metabolito)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ sujeto, scales = "free_y") +
  labs(title = "Abundancia absoluta promedio por fracción y sujeto",
       x = "Fracción de leche", y = "Media de intensidad") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


library(dplyr)
library(ggplot2)

# Calcular porcentaje relativo por fracción dentro de cada sujeto
matriz_relativa_sujeto <- matriz_promedio_sujeto %>%
  group_by(sujeto, Milk_fraction) %>%
  mutate(porcentaje = media / sum(media) * 100) %>%
  ungroup()

# Graficar barras apiladas con facetas por sujeto
ggplot(matriz_relativa_sujeto, aes(x = Milk_fraction, y = porcentaje, fill = metabolito)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ sujeto, scales = "free_y") +
  labs(title = "Abundancia relativa promedio por fracción y sujeto",
       x = "Fracción de leche", y = "Porcentaje (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


# Graficar barras apiladas con facetas por sujeto
ggplot(matriz_relativa_sujeto, aes(x = Milk_fraction, y = porcentaje, fill = metabolito)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ sujeto, scales = "free_y") +
  labs(title = "Abundancia relativa promedio por fracción y sujeto",
       x = "Fracción de leche", y = "Porcentaje (%)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 6),             # tamaño pequeño
    legend.key.size = unit(0.4, "lines"),             # tamaño de caja más pequeño
    legend.position = "right"                        
  ) +
  guides(fill = guide_legend(ncol = 2))                # 2 columnas en la leyenda






# 1. Función para crear SummarizedExperiment agrupado
agrupar_por_clase <- function(se, clase) {
  # Extraer matriz y metadatos
  matriz <- assay(se)
  rowdata <- as.data.frame(rowData(se))
  
  # Agregar la clase deseada como columna de la matriz
  matriz_df <- as.data.frame(matriz) %>%
    mutate(clase = rowdata[[clase]]) %>%
    group_by(clase) %>%
    summarise(across(everything(), sum, na.rm = TRUE)) %>%
    filter(!is.na(clase)) %>%
    column_to_rownames("clase")
  
  # Crear nuevo SE con la matriz agrupada
  SummarizedExperiment(
    assays = list(abundancia = as.matrix(matriz_df)),
    colData = colData(se),
    rowData = DataFrame(nombre_clase = rownames(matriz_df))
  )
}

# 2. Crear los tres objetos
se_superclass <- agrupar_por_clase(se, "Super_class")
se_mainclass  <- agrupar_por_clase(se, "Main_class")
se_subclass   <- agrupar_por_clase(se, "Sub_class")







# Convertir la matriz a data frame y pasar a formato largo
matriz_df_superclass <- as.data.frame(assay(se_superclass)) %>%
  rownames_to_column(var = "metabolito") %>%
  pivot_longer(-metabolito, names_to = "Muestra", values_to = "intensidad")

# Añadir la fracción de leche desde colData
matriz_df_superclass <- matriz_df_superclass %>%
  left_join(as.data.frame(colData(se_superclass)) %>%
              rownames_to_column("Muestra"),
            by = "Muestra")


# Abundancia absoluta promedio por fracción
matriz_promedio_superclass <- matriz_df_superclass %>%
  group_by(Milk_fraction, metabolito) %>%
  summarise(media = mean(intensidad), .groups = "drop")

ggplot(matriz_promedio_superclass, aes(x = Milk_fraction, y = media, fill = metabolito)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundancia absoluta promedio por fracción de leche",
       x = "Fracción de leche", y = "Media de intensidad") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# calcular el porcentaje relativo por fracción
matriz_relativa_frac_superclass <- matriz_promedio_superclass %>%
  group_by(Milk_fraction) %>%
  mutate(porcentaje = media / sum(media) * 100) %>%
  ungroup()

# 3. Graficar
ggplot(matriz_relativa_frac_superclass, aes(x = Milk_fraction, y = porcentaje, fill = metabolito)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundancia relativa promedio por fracción de leche",
       x = "Fracción de leche", y = "Porcentaje (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






# 1. Calcular la media por sujeto, fracción y metabolito
matriz_promedio_sujeto_superclass <- matriz_df_superclass %>%
  group_by(sujeto, Milk_fraction, metabolito) %>%
  summarise(media = mean(intensidad), .groups = "drop")

# 2. Graficar barras apiladas con facet por sujeto
ggplot(matriz_promedio_sujeto_superclass, aes(x = Milk_fraction, y = media, fill = metabolito)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ sujeto, scales = "free_y") +
  labs(title = "Abundancia absoluta promedio por fracción y sujeto",
       x = "Fracción de leche", y = "Media de intensidad") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


library(dplyr)
library(ggplot2)

# Calcular porcentaje relativo por fracción dentro de cada sujeto
matriz_relativa_sujeto_superclass <- matriz_promedio_sujeto_superclass %>%
  group_by(sujeto, Milk_fraction) %>%
  mutate(porcentaje = media / sum(media) * 100) %>%
  ungroup()

# Graficar barras apiladas con facetas por sujeto
ggplot(matriz_relativa_sujeto_superclass, aes(x = Milk_fraction, y = porcentaje, fill = metabolito)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ sujeto, scales = "free_y") +
  labs(title = "Abundancia relativa promedio por fracción y sujeto",
       x = "Fracción de leche", y = "Porcentaje (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))










# Convertir la matriz a data frame y pasar a formato largo
matriz_df_mainclass <- as.data.frame(assay(se_mainclass)) %>%
  rownames_to_column(var = "metabolito") %>%
  pivot_longer(-metabolito, names_to = "Muestra", values_to = "intensidad")

# Añadir la fracción de leche desde colData
matriz_df_mainclass <- matriz_df_mainclass %>%
  left_join(as.data.frame(colData(se_mainclass)) %>%
              rownames_to_column("Muestra"),
            by = "Muestra")


# Abundancia absoluta promedio por fracción
matriz_promedio_mainclass <- matriz_df_mainclass %>%
  group_by(Milk_fraction, metabolito) %>%
  summarise(media = mean(intensidad), .groups = "drop")

ggplot(matriz_promedio_mainclass, aes(x = Milk_fraction, y = media, fill = metabolito)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundancia absoluta promedio por fracción de leche",
       x = "Fracción de leche", y = "Media de intensidad") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# calcular el porcentaje relativo por fracción
matriz_relativa_frac_mainclass <- matriz_promedio_mainclass %>%
  group_by(Milk_fraction) %>%
  mutate(porcentaje = media / sum(media) * 100) %>%
  ungroup()

# 3. Graficar
ggplot(matriz_relativa_frac_mainclass, aes(x = Milk_fraction, y = porcentaje, fill = metabolito)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundancia relativa promedio por fracción de leche",
       x = "Fracción de leche", y = "Porcentaje (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






# Convertir la matriz a data frame y pasar a formato largo
matriz_df_mainclass <- as.data.frame(assay(se_mainclass)) %>%
  rownames_to_column(var = "metabolito") %>%
  pivot_longer(-metabolito, names_to = "Muestra", values_to = "intensidad")

# Añadir la fracción de leche desde colData
matriz_df_mainclass <- matriz_df_mainclass %>%
  left_join(as.data.frame(colData(se_mainclass)) %>%
              rownames_to_column("Muestra"),
            by = "Muestra")


# Abundancia absoluta promedio por fracción
matriz_promedio_mainclass <- matriz_df_mainclass %>%
  group_by(Milk_fraction, metabolito) %>%
  summarise(media = mean(intensidad), .groups = "drop")

ggplot(matriz_promedio_mainclass, aes(x = Milk_fraction, y = media, fill = metabolito)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundancia absoluta promedio por fracción de leche",
       x = "Fracción de leche", y = "Media de intensidad") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# calcular el porcentaje relativo por fracción
matriz_relativa_frac_mainclass <- matriz_promedio_mainclass %>%
  group_by(Milk_fraction) %>%
  mutate(porcentaje = media / sum(media) * 100) %>%
  ungroup()

# 3. Graficar
ggplot(matriz_relativa_frac_mainclass, aes(x = Milk_fraction, y = porcentaje, fill = metabolito)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundancia relativa promedio por fracción de leche",
       x = "Fracción de leche", y = "Porcentaje (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






# 1. Calcular la media por sujeto, fracción y metabolito
matriz_promedio_sujeto_mainclass <- matriz_df_mainclass %>%
  group_by(sujeto, Milk_fraction, metabolito) %>%
  summarise(media = mean(intensidad), .groups = "drop")

# 2. Graficar barras apiladas con facet por sujeto
ggplot(matriz_promedio_sujeto_mainclass, aes(x = Milk_fraction, y = media, fill = metabolito)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ sujeto, scales = "free_y") +
  labs(title = "Abundancia absoluta promedio por fracción y sujeto",
       x = "Fracción de leche", y = "Media de intensidad") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


library(dplyr)
library(ggplot2)

# Calcular porcentaje relativo por fracción dentro de cada sujeto
matriz_relativa_sujeto_mainclass <- matriz_promedio_sujeto_mainclass %>%
  group_by(sujeto, Milk_fraction) %>%
  mutate(porcentaje = media / sum(media) * 100) %>%
  ungroup()

# Graficar barras apiladas con facetas por sujeto
ggplot(matriz_relativa_sujeto_mainclass, aes(x = Milk_fraction, y = porcentaje, fill = metabolito)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ sujeto, scales = "free_y") +
  labs(title = "Abundancia relativa promedio por fracción y sujeto",
       x = "Fracción de leche", y = "Porcentaje (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))







# Convertir la matriz a data frame y pasar a formato largo
matriz_df_subclass  <- as.data.frame(assay(se_subclass )) %>%
  rownames_to_column(var = "metabolito") %>%
  pivot_longer(-metabolito, names_to = "Muestra", values_to = "intensidad")

# Añadir la fracción de leche desde colData
matriz_df_subclass  <- matriz_df_subclass  %>%
  left_join(as.data.frame(colData(se_subclass )) %>%
              rownames_to_column("Muestra"),
            by = "Muestra")


# Abundancia absoluta promedio por fracción
matriz_promedio_subclass  <- matriz_df_subclass  %>%
  group_by(Milk_fraction, metabolito) %>%
  summarise(media = mean(intensidad), .groups = "drop")

ggplot(matriz_promedio_subclass , aes(x = Milk_fraction, y = media, fill = metabolito)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundancia absoluta promedio por fracción de leche",
       x = "Fracción de leche", y = "Media de intensidad") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# calcular el porcentaje relativo por fracción
matriz_relativa_frac_subclass  <- matriz_promedio_subclass  %>%
  group_by(Milk_fraction) %>%
  mutate(porcentaje = media / sum(media) * 100) %>%
  ungroup()

# 3. Graficar
ggplot(matriz_relativa_frac_subclass , aes(x = Milk_fraction, y = porcentaje, fill = metabolito)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundancia relativa promedio por fracción de leche",
       x = "Fracción de leche", y = "Porcentaje (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






# 1. Calcular la media por sujeto, fracción y metabolito
matriz_promedio_sujeto_subclass  <- matriz_df_subclass  %>%
  group_by(sujeto, Milk_fraction, metabolito) %>%
  summarise(media = mean(intensidad), .groups = "drop")

# 2. Graficar barras apiladas con facet por sujeto
ggplot(matriz_promedio_sujeto_subclass , aes(x = Milk_fraction, y = media, fill = metabolito)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ sujeto, scales = "free_y") +
  labs(title = "Abundancia absoluta promedio por fracción y sujeto",
       x = "Fracción de leche", y = "Media de intensidad") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


library(dplyr)
library(ggplot2)

# Calcular porcentaje relativo por fracción dentro de cada sujeto
matriz_relativa_sujeto_subclass  <- matriz_promedio_sujeto_subclass  %>%
  group_by(sujeto, Milk_fraction) %>%
  mutate(porcentaje = media / sum(media) * 100) %>%
  ungroup()

# Graficar barras apiladas con facetas por sujeto
ggplot(matriz_relativa_sujeto_subclass , aes(x = Milk_fraction, y = porcentaje, fill = metabolito)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ sujeto, scales = "free_y") +
  labs(title = "Abundancia relativa promedio por fracción y sujeto",
       x = "Fracción de leche", y = "Porcentaje (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





seleccionar_top_y_otros <- function(se, top_n = 20) {
  matriz <- assay(se)
  
  # Sumar la abundancia total por fila (metabolito o clase)
  suma_total <- rowSums(matriz)
  
  # Seleccionar los top_n más abundantes
  top_elementos <- names(sort(suma_total, decreasing = TRUE))[1:top_n]
  
  # Dividir la matriz
  matriz_top <- matriz[top_elementos, , drop = FALSE]
  matriz_otros <- matriz[!rownames(matriz) %in% top_elementos, , drop = FALSE]
  
  # Sumar el resto como "Others"
  fila_otros <- matrix(colSums(matriz_otros), nrow = 1,
                       dimnames = list("Others", colnames(matriz)))
  
  # Combinar top + others
  nueva_matriz <- rbind(matriz_top, fila_otros)
  
  # Crear nuevo rowData con los nombres nuevos
  nueva_rowData <- DataFrame(nombre = rownames(nueva_matriz))
  
  # Crear el nuevo SummarizedExperiment
  SummarizedExperiment(
    assays = list(abundancia = nueva_matriz),
    colData = colData(se),
    rowData = nueva_rowData
  )
}




se_top <- seleccionar_top_y_otros(se, top_n = 10)


# Convertir la matriz a data frame y pasar a formato largo
matriz_df_top <- as.data.frame(assay(se_top)) %>%
  rownames_to_column(var = "metabolito") %>%
  pivot_longer(-metabolito, names_to = "Muestra", values_to = "intensidad")

# Añadir la fracción de leche desde colData
matriz_df_top <- matriz_df_top %>%
  left_join(as.data.frame(colData(se_top)) %>%
              rownames_to_column("Muestra"),
            by = "Muestra")


# Abundancia absoluta promedio por fracción
matriz_promedio_top <- matriz_df_top %>%
  group_by(Milk_fraction, metabolito) %>%
  summarise(media = mean(intensidad), .groups = "drop")


# calcular el porcentaje relativo por fracción
matriz_relativa_frac_top <- matriz_promedio_top %>%
  group_by(Milk_fraction) %>%
  mutate(porcentaje = media / sum(media) * 100) %>%
  ungroup()


# 1. Calcular la media por sujeto, fracción y metabolito
matriz_promedio_sujeto_top <- matriz_df_top %>%
  group_by(sujeto, Milk_fraction, metabolito) %>%
  summarise(media = mean(intensidad), .groups = "drop")




library(dplyr)
library(ggplot2)

# Calcular porcentaje relativo por fracción dentro de cada sujeto
matriz_relativa_sujeto_top <- matriz_promedio_sujeto_top %>%
  group_by(sujeto, Milk_fraction) %>%
  mutate(porcentaje = media / sum(media) * 100) %>%
  ungroup()

# Reordenar los niveles del factor 'metabolito'
orden_metabolitos <- matriz_df_top %>%
  distinct(metabolito) %>%
  pull(metabolito) %>%
  setdiff("Others") %>%
  c("Others",. )  # Poner "Others" al final

# Aplicar en todos los dataframes usados en gráficos
matriz_df_top$metabolito <- factor(matriz_df_top$metabolito, levels = orden_metabolitos)
matriz_promedio_top$metabolito <- factor(matriz_promedio_top$metabolito, levels = orden_metabolitos)
matriz_relativa_frac_top$metabolito <- factor(matriz_relativa_frac_top$metabolito, levels = orden_metabolitos)
matriz_promedio_sujeto_top$metabolito <- factor(matriz_promedio_sujeto_top$metabolito, levels = orden_metabolitos)
matriz_relativa_sujeto_top$metabolito <- factor(matriz_relativa_sujeto_top$metabolito, levels = orden_metabolitos)

library(scales)

n_metab <- length(orden_metabolitos) - 1  # sin contar "Others"
colores <- c("black", hue_pal()(n_metab))  # "Others" negro
names(colores) <- orden_metabolitos

ggplot(matriz_promedio_top, aes(x = Milk_fraction, y = media, fill = metabolito)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundancia absoluta promedio por fracción de leche",
       x = "Fracción de leche", y = "Media de intensidad") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values = colores)


# 3. Graficar
ggplot(matriz_relativa_frac_top, aes(x = Milk_fraction, y = porcentaje, fill = metabolito)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundancia relativa promedio por fracción de leche",
       x = "Fracción de leche", y = "Porcentaje (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values = colores)

# 2. Graficar barras apiladas con facet por sujeto
ggplot(matriz_promedio_sujeto_top, aes(x = Milk_fraction, y = media, fill = metabolito)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ sujeto, scales = "free_y") +
  labs(title = "Abundancia absoluta promedio por fracción y sujeto",
       x = "Fracción de leche", y = "Media de intensidad") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values = colores)

# Graficar barras apiladas con facetas por sujeto
ggplot(matriz_relativa_sujeto_top, aes(x = Milk_fraction, y = porcentaje, fill = metabolito)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ sujeto, scales = "free_y") +
  labs(title = "Abundancia relativa promedio por fracción y sujeto",
       x = "Fracción de leche", y = "Porcentaje (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_manual(values = colores)
 
