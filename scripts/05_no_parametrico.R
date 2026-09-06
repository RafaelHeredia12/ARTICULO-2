# ============================================================
# 05_no_parametrico.R
# Comparacion: GBN lineal vs. modelo no parametrico (GAM/splines)
# Proyecto: Redes Bayesianas Gaussianas - ERC
# ============================================================
# Una GBN asume que cada nodo hijo es una combinacion LINEAL de
# sus padres. El modelo no parametrico relaja ese supuesto:
# para cada nodo hijo y con padres x_1, ..., x_p, se modela
#
#     y = beta_0 + f_1(x_1) + ... + f_p(x_p) + error
#
# donde cada f_j es una funcion suave (spline) estimada de los
# datos en lugar de un coeficiente lineal fijo. Esto se ajusta
# con gam() del paquete mgcv, tal como se vio en clase.
#
# El BIC/AIC de la red completa es la SUMA de los BIC/AIC de
# las distribuciones locales (un modelo por nodo), ya que la
# log-verosimilitud conjunta se factoriza segun la DAG.
#
# Output:
#   output/figures/efectos_parciales_gam.png
#   data/processed/scores_noparametrico.rds
# ============================================================

library(bnlearn)
library(mgcv)

# ------------------------------------------------------------
# 1. Cargar datos (mismos 214 casos completos que los scripts
#    anteriores, para que los scores sean comparables)
# ------------------------------------------------------------
df  <- readRDS("data/processed/ckd_gbn.rds")
dat <- na.omit(df)
dat[] <- lapply(dat, as.numeric)
cat("Casos completos:", nrow(dat), "\n\n")

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Estructura de la mejor DAG (DAG 1)
#    Se declara como lista de padres por nodo, que es lo que
#    necesitamos para construir las formulas de cada modelo local.
# ------------------------------------------------------------
padres <- list(
  age  = character(0),          # nodo raiz
  bp   = c("age"),
  bgr  = c("age"),
  sc   = c("bp", "bgr"),
  bu   = c("sc"),
  sod  = c("sc"),
  pot  = c("sc"),
  wc   = c("sc"),
  hemo = c("sc"),
  pcv  = c("hemo"),
  rc   = c("hemo")
)

# ------------------------------------------------------------
# 3. Ajustar las distribuciones locales
#    - Nodos con padres  -> gam() con spline s() en cada padre
#    - Nodos raiz        -> lm() con solo intercepto
#
#    NOTA sobre k: s() usa por defecto una base de k = 10. Si una
#    variable tiene menos valores unicos que eso, mgcv falla con
#    "fewer unique covariate combinations than specified maximum
#    degrees of freedom". En este dataset bp solo toma 9 valores
#    distintos (medicion clinica redondeada a decenas), asi que
#    ajustamos k al numero de valores unicos disponibles.
# ------------------------------------------------------------
modelos_gam <- list()
modelos_lin <- list()

# k maximo admisible por variable = (valores unicos - 1), acotado a 10
k_por_var <- sapply(dat, function(x) min(10, length(unique(x)) - 1))
cat("Base k utilizada por variable:\n")
print(k_por_var)
cat("\n")

# Construye el termino spline con k adecuado para esa variable
termino_spline <- function(v) {
  k <- k_por_var[[v]]
  if (k < 10) paste0("s(", v, ", k = ", k, ")") else paste0("s(", v, ")")
}

for (nodo in names(padres)) {
  pa <- padres[[nodo]]
  
  if (length(pa) == 0) {
    # Nodo raiz: sin padres, mismo modelo en ambos enfoques
    f_lin <- as.formula(paste(nodo, "~ 1"))
    modelos_lin[[nodo]] <- lm(f_lin, data = dat)
    modelos_gam[[nodo]] <- lm(f_lin, data = dat)
    cat("Nodo raiz (lm ~ 1):", nodo, "\n")
  } else {
    # Modelo lineal: y ~ x1 + x2 + ...
    f_lin <- as.formula(paste(nodo, "~", paste(pa, collapse = " + ")))
    modelos_lin[[nodo]] <- lm(f_lin, data = dat)
    
    # Modelo no parametrico: y ~ s(x1) + s(x2) + ...
    f_gam <- as.formula(paste(nodo, "~",
                              paste(sapply(pa, termino_spline), collapse = " + ")))
    modelos_gam[[nodo]] <- gam(f_gam, data = dat, method = "REML")
    
    cat("Nodo:", nodo, "| padres:", paste(pa, collapse = ", "),
        "| formula:", deparse(f_gam), "\n")
  }
}

# ------------------------------------------------------------
# 4. BIC y AIC de la red completa = suma sobre nodos
# ------------------------------------------------------------
bic_lineal <- sum(sapply(modelos_lin, BIC))
aic_lineal <- sum(sapply(modelos_lin, AIC))
bic_gam    <- sum(sapply(modelos_gam, BIC))
aic_gam    <- sum(sapply(modelos_gam, AIC))

tabla_np <- data.frame(
  Modelo = c("GBN lineal (DAG 1)",
             "Modelo no parametrico GAM/splines (DAG 1)"),
  BIC    = round(c(bic_lineal, bic_gam), 2),
  AIC    = round(c(aic_lineal, aic_gam), 2)
)

cat("\n=== Comparacion: GBN lineal vs. GAM no parametrico ===\n")
print(tabla_np)
cat("\nNOTA: aqui BIC/AIC se calculan con las funciones base de R,\n")
cat("donde MENOR es mejor (convencion opuesta a score() de bnlearn).\n\n")

if (bic_gam < bic_lineal) {
  cat("El modelo no parametrico MEJORA el BIC.\n")
  cat("Diferencia:", round(bic_lineal - bic_gam, 2), "puntos.\n")
} else {
  cat("El modelo no parametrico NO mejora el BIC.\n")
  cat("Diferencia:", round(bic_gam - bic_lineal, 2), "puntos a favor del lineal.\n")
}

if (aic_gam < aic_lineal) {
  cat("El modelo no parametrico MEJORA el AIC.\n")
} else {
  cat("El modelo no parametrico NO mejora el AIC.\n")
}

# ------------------------------------------------------------
# 5. Grados de libertad efectivos (edf) por nodo
#    edf ~ 1 indica relacion practicamente lineal;
#    edf > 2 sugiere no linealidad relevante.
#    Esto justifica CUALES relaciones se benefician del spline.
# ------------------------------------------------------------
cat("\n=== Grados de libertad efectivos (edf) por termino spline ===\n")
cat("edf cercano a 1 => la relacion es esencialmente lineal\n\n")

filas_edf <- list()
for (nodo in names(modelos_gam)) {
  m <- modelos_gam[[nodo]]
  if (inherits(m, "gam") && length(m$smooth) > 0) {
    s_table <- summary(m)$s.table
    for (i in seq_len(nrow(s_table))) {
      filas_edf[[length(filas_edf) + 1]] <- data.frame(
        Nodo    = nodo,
        Termino = rownames(s_table)[i],
        edf     = round(s_table[i, "edf"], 3),
        p_valor = round(s_table[i, "p-value"], 4)
      )
    }
  }
}
tabla_edf <- do.call(rbind, filas_edf)
print(tabla_edf, row.names = FALSE)

# ------------------------------------------------------------
# 6. Figura: efectos parciales de los splines
# ------------------------------------------------------------
nodos_con_spline <- names(modelos_gam)[
  sapply(modelos_gam, function(m) inherits(m, "gam"))
]

png("output/figures/efectos_parciales_gam.png",
    width = 1600, height = 1200, res = 150)
par(mfrow = c(3, 4), mar = c(4, 4, 3, 1))
for (nodo in nodos_con_spline) {
  plot(modelos_gam[[nodo]], residuals = TRUE, shade = TRUE,
       main = paste("Nodo:", nodo), cex.main = 1)
}
dev.off()
cat("\nImagen guardada: output/figures/efectos_parciales_gam.png\n")

# ------------------------------------------------------------
# 7. Guardar outputs
# ------------------------------------------------------------
saveRDS(tabla_np,  "data/processed/scores_noparametrico.rds")
saveRDS(tabla_edf, "data/processed/edf_gam.rds")
cat("Guardado: scores_noparametrico.rds, edf_gam.rds\n")