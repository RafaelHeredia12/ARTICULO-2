# Análisis Probabilístico de la Enfermedad Renal Crónica mediante Redes Bayesianas Gaussianas

Este repositorio contiene el código fuente, preprocesamiento de datos, modelos probabilísticos y el artículo científico desarrollado para evaluar la interdependencia de biomarcadores renales, comorbilidades y la degradación de la función renal mediante **Redes Bayesianas Gaussianas (GBN)**, **Redes Condicionales Gaussianas (CLG)** y **Modelos Aditivos Generalizados (GAM/Splines)**.

---

## 📌 Estructura del Repositorio

```text
.
├── data/                                  # Gestión e intermediación de datos del proyecto
│   ├── processed/                         # Objetos serializados R y datos procesados
│   │   ├── bic_por_nodo.rds               # Descomposición de la métrica BIC por cada nodo de la red
│   │   ├── ckd_clean.csv                  # Dataset limpio y depurado (formato CSV)
│   │   ├── ckd_clean.rds                  # Dataset limpio preservando tipos de datos R (formato RDS)
│   │   ├── ckd_gbn.csv                    # Dataset filtrado con las 11 variables continuas para GBN
│   │   ├── ckd_gbn.rds                    # Dataset continuo para GBN serializado
│   │   ├── comparacion_bic_aic.rds        # Tabla comparativa de scores BIC y AIC entre estructuras
│   │   ├── cv_gam.rds                     # Resultados de validación cruzada para los modelos GAM
│   │   ├── edf_gam.rds                    # Grados de libertad efectivos (EDF) estimados por el GAM
│   │   ├── mejor_dag_nombre.rds           # Identificador serializado de la estructura óptima (DAG 1)
│   │   ├── mejor_gbn.rds                  # Objeto bn.fit con los parámetros ajustados del modelo óptimo
│   │   ├── mejor_gbn_mixto.rds            # Objeto de la red condicional gaussiana (variables mixtas)
│   │   ├── resultados_queries.rds         # Resultados de las inferencias probabilísticas (Q1 a Q6)
│   │   ├── scores_mixto.rds               # Evaluaciones de ajuste y verosimilitud de la red mixta (CLG)
│   │   └── scores_noparametrico.rds       # Métricas de rendimiento y ajuste para la extensión GAM
│   └── raw/                               # Microdatos originales sin alterar
│       └── ckd.csv                        # Base de datos clínica cruda de Enfermedad Renal Crónica
├── notebooks/                             # Publicación y reporte técnico en Quarto
│   ├── ARTICULO2_RBG.html                 # Reporte final compilado e interactivo para visualización web
│   └── ARTICULO2_RBG.qmd                  # Código fuente del artículo científico en Quarto
├── output/                                # Artefactos y gráficos generados
│   └── figures/                           # Exportación de gráficos y diagramas causales
│       ├── dag1_ckd.png                   # Diagrama del DAG 1 (Convergencia en creatinina sérica)
│       ├── dag2_ckd.png                   # Diagrama del DAG 2 (Filtro glomerular primario)
│       ├── dag3_ckd.png                   # Diagrama del DAG 3 (Cascada hematológica independiente)
│       ├── dag_mixto.png                  # Estructura del modelo condicional-gaussiano (CLG)
│       ├── efectos_parciales_gam.png      # Curvas de efectos no lineales ajustadas por splines
│       └── histogramas_transformacion.png # Evaluaciones de normalidad y distribución de variables
├── scripts/                               # Código modular en R para el pipeline de análisis
│   ├── 01_limpieza_ckd.R                  # Limpieza, imputación, filtrado e inspección de variables
│   ├── 02_comparacion_dags.R              # Definición de estructuras teóricas y evaluación BIC/AIC
│   ├── 03_queries.R                       # Inferencia probabilística por likelihood weighting (Q1 a Q6)
│   ├── 04_categoricas.R                   # Incorporación de comorbilidades discretas (Redes CLG)
│   ├── 05_no_parametrico.R                # Ajuste de splines penados mediante GAM para no linealidades
│   ├── dag1.R                             # Definición estructural e hipótesis clínica del DAG 1
│   ├── dag2.R                             # Definición estructural e hipótesis clínica del DAG 2
│   ├── dag3.R                             # Definición estructural e hipótesis clínica del DAG 3
│   └── queries.R                          # Definición funcional de las preguntas de inferencia
├── .gitignore                             # Archivos y patrones excluidos del control de versiones
├── ARTICULO-2.Rproj                       # Configuración del entorno de trabajo en RStudio
└── README.md                              # Documentación principal del repositorio
```