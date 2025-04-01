
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
library(scales)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(knitr)
library(kableExtra)

# Función para calcular abundancias promedio por sujeto de un objeto SummarizedExperiment
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


# funcion para seleccionar los top n metabolitos y agrupar el resto como "otros"
seleccionar_top_y_otros <- function(se_imputed, top_n = 20) {
  matriz <- assay(se_imputed)
  
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

colores_fraccion <- c(
  "fat"   = "#1f77b4",  # azul
  "skim"  = "#e377c2",  # rosa
  "whole" = "#ff7f0e"   # naranja
)


colores_sujeto <- c(
  "BLS001A" = "#d62728",  # rojo
  "BLS002A" = "#2ca02c",  # verde
  "BLS003A" = "#9467bd",  # púrpura
  "BLS010A" = "#17becf"   # celeste
)


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


#save(se, file = "SummarizedExperiment_dataset.Rda")

# SE REALIZA UNA IMPUTACIÓN CON EL METODO KNN, SIN REMOVER MUESTRAS
se_imputed <- PomaImpute(se, method = "knn", zeros_as_na = TRUE, remove_na = FALSE)

# POR ALGUNA RAZÓN ELIMINA LAS VARIABLES DE LOS METABOLITOS, SE LAS VUELVO A AÑADIR
rowData(se_imputed) <- rowData(se)


# REVISIÓN DEL OBJETO
se_imputed
head(rownames(se_imputed))
head(colnames(se_imputed))
colData(se_imputed)
rowData(se_imputed)
# Ver si hay NAs
sum(is.na(assay(se_imputed)))
# Obtener vector lógico de metabolitos con suma 0 de abundancia
metabolitos_cero <- rowSums(assay(se_imputed)) == 0
sum(metabolitos_cero) # cauntos 0 tiene

table(colData(se_imputed)$Milk_fraction)

table(colData(se_imputed)$sujeto)

#################################################################
### ANALISIS DE ABUNDANCIA ABSOLUTA Y RELATIVA DE LOS METABOLITOS
#################################################################

##ABUNDANCIA TOP 50#######
# aplicar la función
se_top <- seleccionar_top_y_otros(se, top_n = 50)
matriz_promedio_sujeto_top <- calcular_abundancia_por_sujeto(se_top)

# Reordenar los niveles del factor 'metabolito'
orden_metabolitos <- matriz_promedio_sujeto_top %>%
  distinct(metabolito) %>%
  pull(metabolito) %>%
  setdiff("Others") %>%
  c("Others",. )  # Poner "Others" al final

# Aplicar en todos los dataframes usados en gráficos
matriz_promedio_sujeto_top$metabolito <- factor(matriz_promedio_sujeto_top$metabolito, levels = orden_metabolitos)

# Paleta de colores
n_metab <- length(orden_metabolitos) - 1  # sin contar "Others"
colores <- c("black", hue_pal()(n_metab))  # "Others" negro
names(colores) <- orden_metabolitos

# Abundancia absoluta por sujeto y fracción (top metabolitos)
ab_suj_top <- ggplot(matriz_promedio_sujeto_top, 
                     aes(x = Milk_fraction, y = media, fill = metabolito)) +
  geom_bar(stat = "identity") + facet_wrap(~ sujeto, scales = "free_y") +
  labs(title = "Abundancia absoluta por fracción y sujeto (Top metabolitos)",
       x = "Fracción de leche", y = "Media de intensidad") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.5, "lines"), legend.position = "right") +
  guides(fill = guide_legend(ncol = 2))  +
  scale_fill_manual(values = colores)

### ABUNDANCIA PROMEDIO POR CLASE
#Crear los tres objetos
se_mainclass  <- agrupar_por_clase(se_imputed, "Main_class")
matriz_promedio_sujeto_mainclass  <- calcular_abundancia_por_sujeto(se_mainclass)

# graficos para MAIN CLASE
ab_main_suj <- ggplot(matriz_promedio_sujeto_mainclass, 
                      aes(x = Milk_fraction, y = media, fill = metabolito)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ sujeto, scales = "free_y") +
  labs(title = "Abundancia absoluta promedio por fracción y sujeto (Main class)",
       x = "Fracción de leche", y = "Media de intensidad") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
                            legend.key.size = unit(0.5, "lines"), 
                            legend.position = "right") +
  guides(fill = guide_legend(ncol = 2))


ab_suj_top / ab_main_suj + plot_annotation(tag_levels = 'A') 

### NORMALIZAR DATOS

se_normalized <- se_imputed %>% 
  PomaNorm(method = "log_pareto")

rowData(se_normalized) <- rowData(se)

se_normalized
dim(se_normalized)
head(rownames(se_normalized))
head(colnames(se_normalized))
colData(se_normalized)
rowData(se_normalized)


# grafico de densidad de datos
# Densidad antes de la normalización
den_no_normalizado <- PomaDensity(se_imputed, x = "features") +
  labs(title = "Densidad de intensidades antes de la normalización") +
  theme(legend.position = "none")

# Densidad después de la normalización
den_normalizado <- PomaDensity(se_normalized, x = "features") +
  labs(title = "Densidad de intensidades después de la normalización") +
  theme(legend.position = "none")

# Unir gráficos con títulos
den_no_normalizado / den_normalizado +
  plot_annotation(tag_levels = 'A') # para poner etiquetas en los gráficos


# Convertir a formato largo
df_long <- as.data.frame(assay(se_normalized)) %>%
  tibble::rownames_to_column("metabolito") %>%
  pivot_longer(-metabolito, names_to = "muestra", values_to = "valor") %>%
  left_join(as.data.frame(colData(se_normalized)) %>%
              tibble::rownames_to_column("muestra"),
            by = "muestra")

# Gráfico de intensidad por sujeto y fracción
intensidad_normalizada_sujeto <- ggplot(df_long, aes(x = Milk_fraction, y = valor, fill = Milk_fraction)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ sujeto, scales = "free_y") +
  labs(
    title = "Distribución normalizada de intensidades por fracción y muestra",
    x = "Fracción de leche",
    y = "Intensidad normalizada"
  ) +
  theme_minimal() + theme(legend.position = "none") + scale_fill_manual(values = colores_fraccion)


intensidad_normalizada_sujeto



# no hay outliers
se_normalized
PomaOutliers(se_imputed)$data

# Gráfico de outliers
outliers_plot_no_normalizados <- PomaOutliers(se_imputed)$polygon_plot + scale_fill_manual(values = colores_fraccion) +
  theme_minimal() + theme(legend.position = "none")

outliers_plot_normalizados <- PomaOutliers(se_normalized)$polygon_plot + scale_fill_manual(values = colores_fraccion) +
  theme_minimal()

outliers_plot_no_normalizados + outliers_plot_normalizados + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') 

  

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
num_metabolitos_sujeto <- ggplot(metadata, aes(x = Milk_fraction, y = metabolitos_detectados, group = sujeto, color = sujeto)) +
    geom_line(size = 1.1) +
    geom_point(size = 2) +
    labs(title = "Metabolitos detectados por fracción de leche y sujeto",
         x = "Fracción de leche", y = "Número de metabolitos detectados") +
    theme_minimal() + scale_color_manual(values = colores_sujeto)

intensidad_normalizada_sujeto + num_metabolitos_sujeto

# Calcular media, SD y SEM por fracción
resumen_fraccion <- metadata %>%
  group_by(Milk_fraction) %>%
  summarise(
    media = mean(metabolitos_detectados), #calcular la media
    sd = sd(metabolitos_detectados), #calcular la desviación estándar
    n = n(),
    sem = sd / sqrt(n) #calcular el error estándar de la media
  )

#Calcular varianzas y seleccionar los top más variables
varianzas <- apply(matriz, 1, var, na.rm = TRUE) # Calcular varianza por fila
varianzas_filtradas <- varianzas[varianzas > 0] # Filtrar las que no son 0
top_n <- min(30, length(varianzas_filtradas)) # Seleccionar los 30 más variables
top_metabolitos <- names(sort(varianzas_filtradas, decreasing = TRUE))[1:top_n] # Nombres de los top 30
matriz_top_var <- matriz[top_metabolitos, ] # Filtrar la matriz para quedarme con los top 30

# Anotaciones (Milk_fraction y sujeto)
anotaciones <- as.data.frame(colData(se_normalized)) %>% # Extraer metadata
  select(Milk_fraction, sujeto) # Seleccionar las columnas de interés
rownames(anotaciones) <- colnames(se_normalized) # Asignar nombres de filas

# colores para para la fraccion y los sujetos, que se definieron antes
colores_para_heatmap <- list(
  Milk_fraction = colores_fraccion,
  sujeto = colores_sujeto
)

# Heatmap final
heatmap <- pheatmap(matriz_top_var,
           annotation_col = anotaciones,
           annotation_colors = colores_para_heatmap,
           show_rownames = TRUE,
           show_colnames = FALSE,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           main = "Heatmap de los 30 metabolitos más variables",
           fontsize_row = 7)

heatmap



# realizar pca

matriz_pca <- t(assay(se_normalized)) # Transponer la matriz para calcular pca para muestras
pca <- prcomp(matriz_pca, scale. = TRUE)

# Extraer coordenadas del PCA
pca_df <- as.data.frame(pca$x) %>% # Convertir a dataframe
  rownames_to_column("Muestra") %>% # Añadir nombres de filas
  left_join(as.data.frame(colData(se_normalized)) %>% # Añadir metadata
              tibble::rownames_to_column("Muestra"), 
            by = "Muestra") # Unir por la columna Muestra


# Graficar PCA de las muestras
biplot_pca_muestras  <-  ggplot(pca_df, aes(x = PC1, y = PC2, color = sujeto, shape = Milk_fraction, label = Muestra)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3.5, max.overlaps = 100) +
  labs(title = "PCA de muestras basado en intensidades",
       x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "% var)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "% var)")) +
  theme_minimal() +
  theme(legend.position = "right") + scale_color_manual(values = colores_sujeto)

biplot_pca_muestras


