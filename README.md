# Informe PEC1 - Análisis de datos metabolómicos
## Descripción
Este repositorio contiene el análisis exploratorio del estudio ST000957, basado en perfiles metabolómicos de leche humana. 

Se realizó como parte de la asignatura Análisis de Datos Ómicos del máster en Bioestadística y Bioinformática UOC, 2025.

Acceso al estudio original:
https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST000957

## Contenido del repositorio
- `informe.qmd`: Informe en formato Quarto.
- `informe.pdf`: Informe final generado.
- `SummarizedExperiment_dataset.Rda`: Objeto `SummarizedExperiment` con datos y metadatos integrados.
- `script.R`: Script con el análisis principal.
- `metadatos_descripcion.Rmd`: Descripción de variables del estudio.
- `metadatos_descripcion.html`: Versión HTML generada.
- `datos/`: Carpeta con los datos originales en formato texto.
- `informe_files/`: Carpeta con figuras generadas automáticamente por Quarto.
- `references.bib`: Bibliografía utilizada en formato BibTeX.
- `apa.csl` / `vancouver.csl`: Estilos de citación. Se uso APA para el informe final
- `script_pruebas.R`: Script de pruebas, usado para hacer ensayos de código y más analisis que no se incluyeron en el informe final

## Requisitos
- R (versión ≥ 4.2)
- RStudio

## Cómo reproducir el análisis
1. Clona el repositorio:
```bash
git clone https://github.com/tu_usuario/nombre_del_repo.git
