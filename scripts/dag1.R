# ============================================================
# dag1.R
# médico uno (mamá de grethel)
# ============================================================


library(bnlearn)
library(dplyr)
library(Rgraphviz)

# 1. Cargar y preparar datos
df <- read.csv("data/processed/ckd_clean.csv")
datos_bn <- na.omit(df)
datos_bn[] <- lapply(datos_bn, function(x) if (is.factor(x)) factor(x) else x)

# 2. Definir DAG (biomarcadores ERC)
dag <- empty.graph(names(datos_bn))
arcs(dag) <- matrix(c(
  "age", "bp", "age", "bgr",
  "bp", "sc", "bgr", "sc",
  "sc", "bu", "sc", "sod", "sc", "pot", "sc", "hemo",
  "hemo", "pcv", "hemo", "rc"
), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from", "to")))

# 3. Guardar imagen de la DAG
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

png("output/figures/dag1_ckd.png", width = 1200, height = 900, res = 150)
graphviz.plot(dag, layout = "dot", shape = "ellipse")
dev.off()

