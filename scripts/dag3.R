# ============================================================
# dag1.R
# médico tres (prima de Daniel)
# ============================================================


library(bnlearn)
library(dplyr)
library(Rgraphviz)

# 1. Cargar y preparar datos
df <- read.csv("data/processed/ckd_gbn.csv")
datos_bn <- na.omit(df)
datos_bn[] <- lapply(datos_bn, function(x) if (is.factor(x)) factor(x) else as.numeric(x))

# 2. Definir DAG (biomarcadores ERC)
dag <- empty.graph(names(datos_bn))
arcs(dag) <- matrix(c(
  "age",  "bgr",
  "age",  "bp",
  "bgr",  "bp",
  "bgr",  "wc",
  "bp",   "sc",
  "bp",   "pot",
  "sc",   "bu",
  "pot",  "sod",
  "bu",   "hemo",
  "sod",  "hemo",
  "wc",   "hemo",
  "hemo", "pcv",
  "pcv",  "rc"
), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from", "to")))

# 3. Guardar imagen de la DAG
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

png("output/figures/dag_fisiopatologica_ckd.png", width = 1200, height = 900, res = 150)
graphviz.plot(dag, layout = "dot", shape = "ellipse")
dev.off()

score(dag, data = datos_bn, type = "bic-g")
score(dag, data = datos_bn, type = "aic-g")