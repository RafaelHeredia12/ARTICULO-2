# ============================================================
# 04_categoricas.R
# Modelo mixto: extension de la mejor GBN (DAG 1) con
# variables categoricas dm (diabetes) y htn (hipertension)
# Proyecto: Redes Bayesianas Gaussianas - ERC
# ============================================================
# Las variables categoricas solo pueden ser nodos raiz o
# padres de nodos continuos en bnlearn (redes condicional-
# gaussianas). Nunca hijos de nodos continuos.
#
# Arcos nuevos justificados clinicamente:
#   dm  -> bgr  (diabetes eleva glucosa en sangre)
#   dm  -> sc   (diabetes daña riñones, eleva creatinina)
#   htn -> bp   (hipertension es por definicion bp elevada)
#   htn -> sc   (hipertension daña riñones igual que diabetes)
#
# Output:
#   output/figures/dag_mixto.png       -> imagen de la DAG mixta
#   data/processed/mejor_gbn_mixto.rds -> bn.fit del modelo mixto
#   data/processed/scores_mixto.rds    -> BIC/AIC gaussiano vs mixto
# ============================================================

library(bnlearn)
library(Rgraphviz)

# ------------------------------------------------------------
# 1. Cargar datos
# ------------------------------------------------------------
df_full <- readRDS("data/processed/ckd_clean_full.rds")

# Seleccionar las 11 continuas + dm + htn
vars_modelo <- c("age", "bgr", "wc", "bp", "sc", "bu",
                 "pot", "sod", "hemo", "pcv", "rc",
                 "dm", "htn")
dat_mixto <- df_full[, vars_modelo]

# Convertir categoricas a factor
dat_mixto$dm  <- factor(dat_mixto$dm,  levels = c("no", "yes"))
dat_mixto$htn <- factor(dat_mixto$htn, levels = c("no", "yes"))

# Casos completos (mismo criterio que en los scripts anteriores)
dat_mixto <- na.omit(dat_mixto)
cat("Casos completos para el modelo mixto:", nrow(dat_mixto), "\n")

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Definir la DAG mixta
# DAG 1 extendida con dm y htn como nodos raiz
# ------------------------------------------------------------
dag_mixto <- empty.graph(names(dat_mixto))
arcs(dag_mixto) <- matrix(c(
  # Arcos originales DAG 1
  "age",  "bp",
  "age",  "bgr",
  "bp",   "sc",
  "bgr",  "sc",
  "sc",   "bu",
  "sc",   "sod",
  "sc",   "pot",
  "sc",   "hemo",
  "sc",   "wc",
  "hemo", "pcv",
  "hemo", "rc",
  # Arcos nuevos: categoricas -> continuas
  "dm",   "bgr",
  "dm",   "sc",
  "htn",  "bp",
  "htn",  "sc"
), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from", "to")))

# Verificar que es un DAG valido
stopifnot(acyclic(dag_mixto))
cat("DAG mixta valida (aciclica): OK\n")

# ------------------------------------------------------------
# 3. Imagen de la DAG mixta
# ------------------------------------------------------------
# Colorear nodos: categoricos en gris, continuos en blanco
nodos_cat  <- c("dm", "htn")
nodos_cont <- setdiff(names(dat_mixto), nodos_cat)

node_attrs <- list(
  fillcolor = setNames(
    c(rep("lightgrey", length(nodos_cat)),
      rep("white", length(nodos_cont))),
    c(nodos_cat, nodos_cont)
  )
)

png("output/figures/dag_mixto.png", width = 1400, height = 1000, res = 150)
graphviz.plot(dag_mixto, layout = "dot", shape = "ellipse",
              attrs = list(node = list(style = "filled")),
              highlight = list(nodes = nodos_cat,
                               fill  = "lightgrey",
                               col   = "black"))
dev.off()
cat("Imagen guardada: output/figures/dag_mixto.png\n")

# ------------------------------------------------------------
# 4. Ajustar el modelo mixto
# bnlearn ajusta automaticamente como condicional-gaussiana
# cuando detecta nodos factor con hijos continuos
# ------------------------------------------------------------
fit_mixto <- bn.fit(dag_mixto, dat_mixto)
cat("Tipo de nodo dm  :", class(fit_mixto$dm),  "\n")
cat("Tipo de nodo htn :", class(fit_mixto$htn), "\n")
cat("Tipo de nodo bgr :", class(fit_mixto$bgr), "\n")

# ------------------------------------------------------------
# 5. Comparar BIC/AIC: gaussiano puro vs. modelo mixto
# Usamos los mismos casos completos para comparacion justa
# ------------------------------------------------------------

# Subconjunto de casos que tienen dm y htn (para el gaussiano puro)
dat_gauss <- dat_mixto[, c("age","bgr","wc","bp","sc","bu",
                           "pot","sod","hemo","pcv","rc")]
dat_gauss[] <- lapply(dat_gauss, as.numeric)

dag_gauss <- empty.graph(names(dat_gauss))
arcs(dag_gauss) <- matrix(c(
  "age","bp", "age","bgr",
  "bp","sc",  "bgr","sc",
  "sc","bu",  "sc","sod",
  "sc","pot", "sc","hemo",
  "sc","wc",  "hemo","pcv",
  "hemo","rc"
), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from","to")))

bic_gauss <- score(dag_gauss,  data = dat_gauss,  type = "bic-g")
aic_gauss <- score(dag_gauss,  data = dat_gauss,  type = "aic-g")
bic_mixto <- score(dag_mixto,  data = dat_mixto,  type = "bic-cg")
aic_mixto <- score(dag_mixto,  data = dat_mixto,  type = "aic-cg")

tabla_scores <- data.frame(
  Modelo = c("GBN Gaussiana pura (DAG 1)",
             "GBN Condicional-Gaussiana (DAG 1 + dm + htn)"),
  Nodos  = c(11, 13),
  Arcos  = c(nrow(arcs(dag_gauss)), nrow(arcs(dag_mixto))),
  BIC    = round(c(bic_gauss, bic_mixto), 2),
  AIC    = round(c(aic_gauss, aic_mixto), 2)
)

cat("\n=== Comparacion: GBN pura vs. modelo mixto ===\n")
print(tabla_scores)
cat("(mayor = mejor en bnlearn)\n\n")

# ------------------------------------------------------------
# 6. Inferencia de ejemplo con el modelo mixto
# ¿Cual es la probabilidad de sc > 1.5 dado dm = yes y htn = yes?
# Muestra como las categoricas enriquecen la inferencia
# ------------------------------------------------------------
set.seed(123)
cat("Ejemplo de inferencia con modelo mixto:\n")
cat("P(sc > 1.5 | dm = yes, htn = yes):\n")
p_mixto <- cpquery(fit_mixto,
                   event    = (sc > 1.5),
                   evidence = list(dm = "yes", htn = "yes"),
                   method   = "lw", n = 1e6)
cat("Resultado:", round(p_mixto, 4), "\n\n")

cat("P(sc > 1.5 | dm = no, htn = no):\n")
p_sano <- cpquery(fit_mixto,
                  event    = (sc > 1.5),
                  evidence = list(dm = "no", htn = "no"),
                  method   = "lw", n = 1e6)
cat("Resultado:", round(p_sano, 4), "\n")
cat("-> La diferencia entre ambos ilustra el efecto de dm y htn sobre sc.\n\n")

# ------------------------------------------------------------
# 7. Guardar outputs
# ------------------------------------------------------------
saveRDS(fit_mixto,    "data/processed/mejor_gbn_mixto.rds")
saveRDS(tabla_scores, "data/processed/scores_mixto.rds")
cat("Guardado: mejor_gbn_mixto.rds, scores_mixto.rds\n")