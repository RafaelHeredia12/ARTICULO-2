# ============================================================
# 02_comparacion_dags.R
# Compara las 3 DAGs propuestas por especialistas
# Proyecto: Redes Bayesianas Gaussianas - ERC
# ============================================================
# Output:
#   output/figures/dag{1,2,3}_ckd.png      -> imagenes de las 3 DAGs
#   data/processed/comparacion_bic_aic.rds -> tabla de scores
#   data/processed/mejor_gbn.rds           -> bn.fit de la mejor estructura
#   data/processed/mejor_dag_nombre.rds    -> string con el nombre de la mejor
# ============================================================

library(bnlearn)
library(Rgraphviz)

# ------------------------------------------------------------
# 1. Cargar datos
# ------------------------------------------------------------
df  <- readRDS("data/processed/ckd_gbn.rds")
dat <- na.omit(df)
dat[] <- lapply(dat, as.numeric)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Definir las 3 DAGs de especialistas
# ------------------------------------------------------------

# --- DAG 1: mamá de Grethel ---
dag1 <- empty.graph(names(dat))
arcs(dag1) <- matrix(c(
  "age", "bp",
  "age", "bgr",
  "bp",  "sc",
  "bgr", "sc",
  "sc",  "bu",
  "sc",  "sod",
  "sc",  "pot",
  "sc",  "hemo",
  "sc",  "wc",
  "hemo","pcv",
  "hemo","rc"
), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from","to")))

# --- DAG 2: tía de Matías ---
dag2 <- empty.graph(names(dat))
arcs(dag2) <- matrix(c(
  "age", "bp",
  "bp",  "sc",
  "bgr", "sc",
  "wc",  "sc",
  "sc",  "bu",
  "sc",  "sod",
  "sc",  "pot",
  "sc",  "hemo",
  "hemo","pcv",
  "hemo","rc"
), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from","to")))

# --- DAG 3: prima de Daniel ---
dag3 <- empty.graph(names(dat))
arcs(dag3) <- matrix(c(
  "age", "bgr",
  "age", "bp",
  "bgr", "bp",
  "bgr", "wc",
  "bp",  "sc",
  "bp",  "pot",
  "sc",  "bu",
  "pot", "sod",
  "bu",  "hemo",
  "sod", "hemo",
  "wc",  "hemo",
  "hemo","pcv",
  "pcv", "rc"
), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from","to")))

# ------------------------------------------------------------
# 3. Imagenes de las 3 DAGs
# ------------------------------------------------------------
dags           <- list(dag1 = dag1, dag2 = dag2, dag3 = dag3)
names_friendly <- c("DAG 1 (Mamá de Grethel)",
                    "DAG 2 (Tía de Matías)",
                    "DAG 3 (Prima de Daniel)")

for (nombre in names(dags)) {
  png(paste0("output/figures/", nombre, "_ckd.png"),
      width = 1200, height = 900, res = 150)
  graphviz.plot(dags[[nombre]], layout = "dot", shape = "ellipse")
  dev.off()
}

# ------------------------------------------------------------
# 4. Calcular BIC-g y AIC-g (mismo n para las 3)
# ------------------------------------------------------------
bic <- sapply(dags, function(d) score(d, data = dat, type = "bic-g"))
aic <- sapply(dags, function(d) score(d, data = dat, type = "aic-g"))

tabla <- data.frame(
  DAG   = names_friendly,
  Arcos = sapply(dags, function(d) nrow(arcs(d))),
  BIC   = round(bic, 2),
  AIC   = round(aic, 2),
  row.names = NULL
)

tabla <- tabla[order(tabla$BIC, decreasing = TRUE), ]

cat("\n=== Comparación BIC / AIC (mayor = mejor en bnlearn) ===\n")
print(tabla)

# ------------------------------------------------------------
# 5. Identificar la mejor estructura y ajustar bn.fit
# ------------------------------------------------------------
mejor_nombre <- names(dags)[which.max(bic)]
mejor_dag    <- dags[[mejor_nombre]]
mejor_fit    <- bn.fit(mejor_dag, dat)

cat("\nMejor estructura:", mejor_nombre, "\n")
cat("BIC:", round(max(bic), 2), " | AIC:", round(aic[which.max(bic)], 2), "\n")

# ------------------------------------------------------------
# 6. Guardar outputs para los scripts siguientes y el articulo
# ------------------------------------------------------------
saveRDS(tabla,        "data/processed/comparacion_bic_aic.rds")
saveRDS(mejor_fit,    "data/processed/mejor_gbn.rds")
saveRDS(mejor_nombre, "data/processed/mejor_dag_nombre.rds")

cat("\nGuardado: comparacion_bic_aic.rds, mejor_gbn.rds, mejor_dag_nombre.rds\n")


