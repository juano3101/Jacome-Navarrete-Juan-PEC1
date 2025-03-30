
## EL ESTUDIO ESCOGIDO FUE DE LA BASE DE DATOS DE METABOLOMICS WORKBENCH ST000957
#ion positive mod

### CARGAR LIBRERIAS
library(metabolomicsWorkbenchR)
library(POMA)
library(SummarizedExperiment)
library(readr)    
library(tibble)    
library(S4Vectors)  
library(dplyr)
library(ggplot2)
library(patchwork)

# RESUMEN DEL ESTUDIO
t(do_query('study','study_id','ST000957','summary'))


# CARGAR MATRICES DE DATOS
assay <- read.delim("datos/assay.tsv", row.names = 1, check.names = FALSE)
coldata <- read.delim("datos/coldata.txt", row.names = 1, check.names = FALSE)
rowdata <- read.delim("datos/rowdata.tsv", row.names = 1, check.names = FALSE)

# REORDENAR EL COLDATA PARA QUE COINCIDA CON LAS FILAS DE LA MATRIZ
coldata <- coldata[colnames(assay), , drop = FALSE]

# AÑADIR LA COLUMNA DE SUJETO EN COLDATA
coldata$sujeto <- gsub("_.*", "", rownames(coldata))


# CREAR EL OBJETO SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(count = assay),
  colData = coldata,
  rowData = rowdata
)


# CONSERVAR SOLO LAS COLUMAS QUE SE DESAN
colData(se) <- colData(se)[, c("Milk_fraction", "sujeto")]
rowData(se) <- rowData(se)[, c("Standardized_name", "Formula", "Super_class", "Main_class", "Sub_class")]


# CREAR UN VECTOR CON METABOLITOS QUE NO ESTAN DENTRO DE LA BASE
# DE Metabolomics Workbench
sin_nombre_estandarizado <- rowData(se)$Standardized_name == "-"
# ELIMINAR ESOS METABOLITOS
se <- se[!sin_nombre_estandarizado, ]

# CONVERTIR EN FACTORES A LAS VARIABLES CATEGÓRICAS
colData(se)$Milk_fraction <- as.factor(colData(se)$Milk_fraction)
colData(se)$sujeto <- as.factor(colData(se)$sujeto)


# SE REALIZA UNA IMPUTACIÓN CON EL METODO KNN, SIN REMOVER MUESTRAS
se_imputed <- PomaImpute(se, method = "knn", zeros_as_na = TRUE, remove_na = FALSE)

# POR ALGUNA RAZÓN ELIMINA LAS VARIABLES DE LOS METABOLITOS, SE LAS VUELVO A AÑADIR
rowData(se_imputed) <- rowData(se)


# REVISIÓN DEL OBJETO
se_imputed
dim(se_imputed)
assayNames(se_imputed)
head(rownames(se_imputed))
head(colnames(se_imputed))
colData(se_imputed)
rowData(se_imputed)
# Ver si hay NAs
sum(is.na(assay(se_imputed)))
# Obtener vector lógico de metabolitos con suma 0 de abundancia
metabolitos_cero <- rowSums(assay(se_imputed)) == 0
sum(metabolitos_cero) # cauntos 0 tiene


#################################################################
### ANALISIS DE ABUNDANCIA ABSOLUTA Y RELATIVA DE LOS METABOLITOS
#################################################################


### ABUNDANCIA ABSOLUTA PROMEDIO POR FRACCIÓN DE LECHE

# Función para calcular abundancias promedio de un objeto SummarizedExperiment
calcular_abundancias_promedio <- function(se_objeto) {
  # Convertir a formato largo
  matriz_df <- as.data.frame(assay(se_objeto)) %>%
    rownames_to_column(var = "metabolito") %>%
    pivot_longer(-metabolito, names_to = "Muestra", values_to = "intensidad") %>%
    left_join(as.data.frame(colData(se_objeto)) %>% # Añadir la fracción de leche desde colData
                rownames_to_column("Muestra"),
              by = "Muestra")
  
  # Calcular media y porcentaje en un solo paso
  matriz_resumen <- matriz_df %>%
    group_by(Milk_fraction, metabolito) %>%
    summarise(
      media = mean(intensidad), #calcular la media
      .groups = "drop_last" # no mostrar el último grupo
    ) %>%
    mutate(porcentaje = media / sum(media) * 100) %>% #calcular el porcentaje
    ungroup() # desagrupar
  
  return(matriz_resumen)
}


resumen_se_imputed <- calcular_abundancias_promedio(se_imputed)

# Grafico de abudancia promedio por fracción
mean_ab_frac <- ggplot(resumen_se_imputed, aes(x = Milk_fraction, y = media, fill = metabolito)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundancia absoluta promedio por fracción de leche",
       x = "Fracción de leche", y = "Media de intensidad") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") #+
#  theme_minimal() +
#  theme(
#    axis.text.x = element_text(angle = 45, hjust = 1),
#    legend.text = element_text(size = 6),             # tamaño pequeño
#    legend.key.size = unit(0.4, "lines"),             # tamaño de caja más pequeño
#    legend.position = "right"                        
#  ) +
#  guides(fill = guide_legend(ncol = 2))                # 2 columnas en la leyenda

# Grafico de abudancia relativa promedio por fracción
mean_rel_frac <- ggplot(resumen_se_imputed, aes(x = Milk_fraction, y = porcentaje, fill = metabolito)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundancia relativa promedio por fracción de leche",
       x = "Fracción de leche", y = "Porcentaje (%)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 8),             # tamaño pequeño
    legend.key.size = unit(0.5, "lines"),             # tamaño de caja más pequeño
    legend.position = "right"                        
  ) +
  guides(fill = guide_legend(ncol = 2))                # 2 columnas en la leyenda

# usar el paquete patchwork para unir los dos gráficos
# plot_layout(guides = "collect") usa una de las leyendas para los dos gráficos
# (mean_ab_frac)/(mean_rel_frac) pone un gráficos debajo del otro
abundancia_frac <- (mean_ab_frac)/(mean_rel_frac) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') # para poner etiquetas en los gráficos




### ABUNDANCIA ABSOLUTA PROMEDIO POR SUJETO

calcular_abundancia_por_sujeto <- function(se_objeto) {
  # Convertir a formato largo
  matriz_df <- as.data.frame(assay(se_objeto)) %>% 
    rownames_to_column(var = "metabolito") %>% 
    pivot_longer(-metabolito, names_to = "Muestra", values_to = "intensidad") %>% 
    left_join(as.data.frame(colData(se_objeto)) %>% # Añadir la fracción de leche desde colData
                rownames_to_column("Muestra"),
              by = "Muestra")
  
  # Calcular media por sujeto y fracción
  matriz_resumen <- matriz_df %>% 
    group_by(sujeto, Milk_fraction, metabolito) %>% # Agrupar por sujeto, fracción y metabolito
    summarise(media = mean(intensidad), .groups = "drop") # Calcular media
  
  # Calcular porcentaje relativo dentro de cada fracción y sujeto
  matriz_resumen <- matriz_resumen %>%
    group_by(sujeto, Milk_fraction) %>% # Agrupar por sujeto y fracción
    mutate(porcentaje = media / sum(media) * 100) %>% # Calcular porcentaje
    ungroup()
  
  return(matriz_resumen)
}

matriz_promedio_sujeto <- calcular_abundancia_por_sujeto(se_imputed)

# Graficar barras apiladas con facet por sujeto
mean_ab_sujeto <- ggplot(matriz_promedio_sujeto, aes(x = Milk_fraction, y = media, fill = metabolito)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ sujeto, scales = "free_y") +
  labs(title = "Abundancia absoluta promedio por fracción y sujeto",
       x = "Fracción de leche", y = "Media de intensidad") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Graficar barras apiladas con facetas por sujeto
mean_rel_sujeto <- ggplot(matriz_promedio_sujeto, aes(x = Milk_fraction, y = porcentaje, fill = metabolito)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ sujeto, scales = "free_y") +
  labs(title = "Abundancia relativa promedio por fracción y sujeto",
       x = "Fracción de leche", y = "Porcentaje (%)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 8),             # tamaño pequeño
    legend.key.size = unit(0.5, "lines"),             # tamaño de caja más pequeño
    legend.position = "right"                        
  ) +
  guides(fill = guide_legend(ncol = 2))                # 2 columnas en la leyenda

# usar el paquete patchwork para unir los dos gráficos
abundancia_sujeto <- (mean_ab_sujeto)/(mean_rel_sujeto) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') # para poner etiquetas en los gráficos



### ABUNDANCIA ABSOLUTA PROMEDIO POR CLASE

# Función para crear SummarizedExperiment agrupado
agrupar_por_clase <- function(se, clase) {
  # Extraer matriz y metadatos
  matriz <- assay(se)
  rowdata <- as.data.frame(rowData(se)) 
  
  # Agregar la clase deseada como columna de la matriz
  matriz_df <- as.data.frame(matriz) %>%
    mutate(clase = rowdata[[clase]]) %>% # Agregar la clase
    group_by(clase) %>% # Agrupar por clase
    summarise(across(everything(), sum, na.rm = TRUE)) %>% # Sumar por clase
    filter(!is.na(clase)) %>% # Eliminar NA
    column_to_rownames("clase") # Usar la clase como fila
  
  # Crear nuevo SE con la matriz agrupada
  SummarizedExperiment(
    assays = list(abundancia = as.matrix(matriz_df)),
    colData = colData(se),
    rowData = DataFrame(nombre_clase = rownames(matriz_df))
  )
}

#Crear los tres objetos
se_superclass <- agrupar_por_clase(se_imputed, "Super_class")
se_mainclass  <- agrupar_por_clase(se_imputed, "Main_class")
se_subclass   <- agrupar_por_clase(se_imputed, "Sub_class")

# Ejecutar para cada nivel con la funcion calcular_abundancias_promedio creada anteriormente
resumen_superclass <- calcular_abundancias_promedio(se_superclass)
resumen_mainclass  <- calcular_abundancias_promedio(se_mainclass)
resumen_subclass   <- calcular_abundancias_promedio(se_subclass)


ggplot(res_superclass, aes(x = Milk_fraction, y = media, fill = metabolito)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundancia absoluta promedio por fracción de leche",
       x = "Fracción de leche", y = "Media de intensidad") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# 3. Graficar
ggplot(res_superclass, aes(x = Milk_fraction, y = porcentaje, fill = metabolito)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundancia relativa promedio por fracción de leche",
       x = "Fracción de leche", y = "Porcentaje (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(res_mainclass, aes(x = Milk_fraction, y = media, fill = metabolito)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundancia absoluta promedio por fracción de leche",
       x = "Fracción de leche", y = "Media de intensidad") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3. Graficar
ggplot(res_mainclass, aes(x = Milk_fraction, y = porcentaje, fill = metabolito)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundancia relativa promedio por fracción de leche",
       x = "Fracción de leche", y = "Porcentaje (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(res_subclass , aes(x = Milk_fraction, y = media, fill = metabolito)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundancia absoluta promedio por fracción de leche",
       x = "Fracción de leche", y = "Media de intensidad") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(res_subclass , aes(x = Milk_fraction, y = porcentaje, fill = metabolito)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundancia relativa promedio por fracción de leche",
       x = "Fracción de leche", y = "Porcentaje (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

















# Ejecutar para cada nivel
matriz_promedio_sujeto_superclass <- calcular_abundancia_por_sujeto(se_superclass)
matriz_promedio_sujeto_mainclass  <- calcular_abundancia_por_sujeto(se_mainclass)
matriz_promedio_sujeto_subclass   <- calcular_abundancia_por_sujeto(se_subclass)





# 2. Graficar barras apiladas con facet por sujeto
ggplot(matriz_promedio_sujeto_superclass, aes(x = Milk_fraction, y = media, fill = metabolito)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ sujeto, scales = "free_y") +
  labs(title = "Abundancia absoluta promedio por fracción y sujeto",
       x = "Fracción de leche", y = "Media de intensidad") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

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




se_top <- seleccionar_top_y_otros(se, top_n = 20)


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

























se_normalized <- se_imputed %>% 
  PomaNorm(method = "log_pareto")

rowData(se_normalized) <- rowData(se)


dim(se_normalized)


se_normalized
dim(se_normalized)
assayNames(se_normalized)
head(rownames(se_normalized))
head(colnames(se_normalized))
colData(se_normalized)
rowData(se_normalized)


table(colData(se_normalized)$Milk_fraction)


PomaBoxplots(se_imputed, x = "samples") # data before normalization
PomaBoxplots(se_normalized, x = "samples") # data after normalization


PomaDensity(se_imputed, x = "features") + # data after normalization  +
  theme(legend.position = "none")

PomaDensity(se_normalized, x = "features") + # data after normalization  +
  theme(legend.position = "none")


# no hay outliers
se_normalized
PomaOutliers(se_normalized)$data

library(dplyr)
library(tibble)
library(ggplot2)


# Extraer info y preparar
matriz <- assay(se_normalized)
metadata <- as.data.frame(colData(se_normalized)) %>%
  rownames_to_column("Muestra")

# Calcular metabolitos detectados (abundancia > 0) por muestra
# Calcular número de metabolitos detectados por muestra (2 es para cada col.) 
metabolitos_por_muestra <- apply(matriz, 2, function(x) sum(x > 0))

# Agregar al metadata
metadata$metabolitos_detectados <- metabolitos_por_muestra

colData(se_normalized)[["metabolitos_detectados"]] <- metabolitos_por_muestra

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
suma_por_muestra <- colSums(assay(se_normalized))

metadata$intensidad_por_muestra <- suma_por_muestra
colData(se_normalized)[["intensidad_por_muestra"]] <- suma_por_muestra

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
matriz_df_normalized <- as.data.frame(assay(se_normalized)) %>%
  rownames_to_column(var = "metabolito") %>%
  pivot_longer(-metabolito, names_to = "Muestra", values_to = "intensidad")

# Añadir la fracción de leche desde colData
matriz_df_normalized <- matriz_df_normalized %>%
  left_join(as.data.frame(colData(se_normalized)) %>%
              rownames_to_column("Muestra"),
            by = "Muestra")

ggplot(matriz_df, aes(x = Milk_fraction, y = intensidad, fill = sujeto)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  facet_wrap(~ sujeto) +
  labs(
    title = "Distribución de intensidades normalizadas por fracción y sujeto",
    y = "Intensidad (log pareto normalizada)",
    x = "Fracción de leche"
  ) +
  theme_minimal()   +
  theme(legend.position = "none")











library(pheatmap)
library(RColorBrewer)

# 1. Extraer matriz de intensidades
matriz <- assay(se)

# 2. Calcular varianzas y seleccionar los top 30 más variables (con varianza > 0)
varianzas <- apply(matriz, 1, var, na.rm = TRUE)
varianzas_filtradas <- varianzas[varianzas > 0]
top_n <- min(30, length(varianzas_filtradas))
top_metabolitos <- names(sort(varianzas_filtradas, decreasing = TRUE))[1:top_n]
matriz_top_var <- matriz[top_metabolitos, ]

# 3. Escalar por fila (z-score)
matriz_escalada <- t(scale(t(matriz_top_var)))

# 4. Anotaciones: sujeto y fracción
anotaciones <- as.data.frame(colData(se_top)) %>%
  select(Milk_fraction, sujeto)
rownames(anotaciones) <- colnames(se_top)

# 5. Colores automáticos para anotaciones
# Generar paletas de colores distintas para cada variable categórica
colors_fractions <- setNames(brewer.pal(length(unique(anotaciones$Milk_fraction)), "Pastel1"),
                             unique(anotaciones$Milk_fraction))

colors_sujetos <- setNames(brewer.pal(n = length(unique(anotaciones$sujeto)), name = "Set2"),
                           unique(anotaciones$sujeto))

ann_colors <- list(
  Milk_fraction = colors_fractions,
  sujeto = colors_sujetos
)

# 6. Heatmap final
pheatmap(matriz_escalada,
         annotation_col = anotaciones,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Heatmap de los 30 metabolitos más variables",
         fontsize_row = 7)





library(ggplot2)
library(ggrepel)

# 1. Obtener matriz de intensidades (muestras en filas)
matriz <- t(assay(se))

# 2. Realizar PCA
pca <- prcomp(matriz, scale. = TRUE)

# 3. Extraer coordenadas del PCA
pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("Muestra") %>%
  left_join(as.data.frame(colData(se_top)) %>%
              tibble::rownames_to_column("Muestra"),
            by = "Muestra")

# 4. Graficar PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = sujeto, shape = Milk_fraction, label = Muestra)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3.5, max.overlaps = 100) +
  labs(title = "PCA de muestras basado en intensidades",
       subtitle = "Colores: sujeto | Forma: fracción | Etiqueta: ID muestra",
       x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "% var)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "% var)")) +
  theme_minimal() +
  theme(legend.position = "right")











# 1. Obtener matriz: metabolitos en filas
matriz <- assay(se)

# 2. PCA con escalado
pca_metab <- prcomp(matriz, scale. = TRUE)

# 3. Coordenadas de metabolitos (puntos)
metab_coords <- as.data.frame(pca_metab$x) %>%
  rownames_to_column("metabolito")

# 4. Coordenadas de muestras (vectores)
muestra_coords <- as.data.frame(pca_metab$rotation) %>%
  rownames_to_column("Muestra") %>%
  left_join(as.data.frame(colData(se_top)) %>%
              rownames_to_column("Muestra"),
            by = "Muestra")

# 5. Gráfico PCA biplot
ggplot() +
  # Puntos: metabolitos
  geom_point(data = metab_coords, aes(x = PC1, y = PC2), color = "black", alpha = 0.7, size = 2) +
  geom_text_repel(data = metab_coords, aes(x = PC1, y = PC2, label = metabolito), size = 3, max.overlaps = 50) +
  
  # Flechas: muestras
  geom_segment(data = muestra_coords,
               aes(x = 0, y = 0, xend = PC1 * 5, yend = PC2 * 5, color = sujeto),
               arrow = arrow(length = unit(0.2, "cm")), alpha = 0.8) +
  geom_label_repel(data = muestra_coords,
                   aes(x = PC1 * 5, y = PC2 * 5, label = Muestra, fill = Milk_fraction),
                   size = 3.2, color = "white", show.legend = TRUE) +
  
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Pastel1") +
  labs(title = "PCA biplot: Metabolitos vs. Muestras",
       subtitle = "Metabolitos como puntos (negro), muestras como flechas (color: sujeto, relleno: fracción)",
       x = paste0("PC1 (", round(summary(pca_metab)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_metab)$importance[2, 2] * 100, 1), "%)"),
       color = "Sujeto", fill = "Fracción de leche") +
  theme_minimal()

