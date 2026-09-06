# ============================================================
# 05_no_parametrico.R
# Comparacion: GBN gaussiana vs. GBN con variables transformadas
# Proyecto: Redes Bayesianas Gaussianas - ERC
# ============================================================
# Las GBNs asumen que cada nodo sigue una distribucion normal
# condicional. Cuando las variables no cumplen este supuesto,
# una alternativa no parametrica es transformar las variables
# para acercarse a la normalidad (transformacion logaritmica),
# ajustar la GBN sobre los datos transformados, y comparar
# el BIC/AIC contra el modelo gaussiano original.
#
# Output:
#   output/figures/histogramas_transformacion.png
#   data/processed/scores_noparametrico.rds
# ============================================================

library(bnlearn)

# ------------------------------------------------------------
# 1. Cargar datos (mismos casos completos que scripts anteriores)
# ------------------------------------------------------------
df  <- readRDS("data/processed/ckd_gbn.rds")
dat <- na.omit(df)
dat[] <- lapply(dat, as.numeric)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Explorar normalidad de las variables originales
# Prueba de Shapiro-Wilk (H0: la variable sigue distribucion normal)
# ------------------------------------------------------------
cat("=== Prueba de Shapiro-Wilk (p < 0.05 rechaza normalidad) ===\n")
shapiro_original <- sapply(dat, function(x) shapiro.test(x)$p.value)
print(round(sort(shapiro_original), 4))

vars_no_normales <- names(shapiro_original[shapiro_original < 0.05])
cat("\nVariables que rechazan normalidad:", paste(vars_no_normales, collapse = ", "), "\n\n")

# ------------------------------------------------------------
# 3. Transformacion logaritmica en variables no normales
# con valores positivos (log(x) no esta definido para x <= 0)
# ------------------------------------------------------------
dat_trans <- dat

# Verificar minimos para decidir si aplicar log o log1p
mins <- sapply(dat[vars_no_normales], min, na.rm = TRUE)
cat("Minimos de variables no normales:\n")
print(round(mins, 3))

for (v in vars_no_normales) {
  if (min(dat[[v]], na.rm = TRUE) > 0) {
    dat_trans[[v]] <- log(dat[[v]])
    cat("log() aplicado a:", v, "\n")
  } else {
    dat_trans[[v]] <- log1p(dat[[v]])  # log(1+x) para valores >= 0
    cat("log1p() aplicado a:", v, "\n")
  }
}

# ------------------------------------------------------------
# 4. Verificar mejora en normalidad post-transformacion
# ------------------------------------------------------------
cat("\n=== Shapiro-Wilk post-transformacion ===\n")
shapiro_trans <- sapply(dat_trans[vars_no_normales],
                        function(x) shapiro.test(x)$p.value)
comparacion_shapiro <- data.frame(
  Variable     = vars_no_normales,
  p_original   = round(shapiro_original[vars_no_normales], 4),
  p_trans      = round(shapiro_trans, 4),
  Mejora       = ifelse(shapiro_trans > shapiro_original[vars_no_normales],
                        "Si", "No")
)
print(comparacion_shapiro)

# ------------------------------------------------------------
# 5. Histogramas: original vs. transformado (para el articulo)
# ------------------------------------------------------------
vars_plot <- vars_no_normales[1:min(4, length(vars_no_normales))]
png("output/figures/histogramas_transformacion.png",
    width = 1600, height = 800, res = 150)
par(mfrow = c(2, length(vars_plot)), mar = c(4, 4, 3, 1))
for (v in vars_plot) {
  hist(dat[[v]],      main = paste("Original:", v),      xlab = v,
       col = "steelblue", border = "white")
  hist(dat_trans[[v]], main = paste("Transformada:", v), xlab = paste("log(", v, ")"),
       col = "darkorange", border = "white")
}
dev.off()
cat("\nImagen guardada: output/figures/histogramas_transformacion.png\n")

# ------------------------------------------------------------
# 6. Ajustar la misma DAG 1 sobre datos transformados
# y comparar BIC/AIC con el modelo gaussiano original
#
# IMPORTANTE: las verosimilitudes de modelos ajustados sobre
# escalas distintas NO son directamente comparables. Al aplicar
# Y = log(X), la densidad cambia segun el jacobiano de la
# transformacion:
#     f_Y(y) = f_X(x) * |dx/dy|,  con dx/dy = x  para y = log(x)
# Por lo tanto, para comparar en la escala ORIGINAL hay que
# restar sum(log(x)) a la log-verosimilitud del modelo
# transformado. Sin esta correccion, el BIC del modelo
# transformado parece artificialmente mejor solo por el cambio
# de unidades (los valores se comprimen y la densidad sube).
# ------------------------------------------------------------
dag1 <- empty.graph(names(dat))
arcs(dag1) <- matrix(c(
  "age", "bp",  "age", "bgr",
  "bp",  "sc",  "bgr", "sc",
  "sc",  "bu",  "sc",  "sod",
  "sc",  "pot", "sc",  "hemo",
  "sc",  "wc",  "hemo","pcv",
  "hemo","rc"
), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from","to")))

bic_original <- score(dag1, data = dat,       type = "bic-g")
aic_original <- score(dag1, data = dat,       type = "aic-g")

# Scores del modelo transformado EN SU PROPIA ESCALA (no comparables)
bic_trans_escala <- score(dag1, data = dat_trans, type = "bic-g")
aic_trans_escala <- score(dag1, data = dat_trans, type = "aic-g")

# --- Correccion jacobiana ---
# Para cada variable transformada con log(), el termino de
# correccion por observacion es -log(x). Sumamos sobre todas las
# variables transformadas y todas las observaciones.
log_jacobiano <- sum(sapply(vars_no_normales, function(v) sum(log(dat[[v]]))))

bic_trans_corregido <- bic_trans_escala - log_jacobiano
aic_trans_corregido <- aic_trans_escala - log_jacobiano

cat("\n=== Correccion jacobiana ===\n")
cat("Termino de correccion sum(log(x)):", round(log_jacobiano, 2), "\n")
cat("BIC transformado sin corregir (escala log):", round(bic_trans_escala, 2), "\n")
cat("BIC transformado CORREGIDO (escala original):", round(bic_trans_corregido, 2), "\n\n")

tabla_np <- data.frame(
  Modelo = c("GBN gaussiana (datos originales)",
             "GBN gaussiana (log, sin corregir - NO comparable)",
             "GBN gaussiana (log, con correccion jacobiana)"),
  BIC    = round(c(bic_original, bic_trans_escala, bic_trans_corregido), 2),
  AIC    = round(c(aic_original, aic_trans_escala, aic_trans_corregido), 2),
  Comparable = c("Si", "No", "Si")
)

cat("=== Comparacion BIC/AIC ===\n")
print(tabla_np)
cat("(mayor = mejor en bnlearn; comparar solo las filas marcadas 'Si')\n\n")

if (bic_trans_corregido > bic_original) {
  cat("La transformacion logaritmica MEJORA el ajuste (BIC corregido).\n")
  cat("Interpretacion: al relajar la desviacion de normalidad, el modelo\n")
  cat("captura mejor la estructura de los datos aun penalizando complejidad.\n")
} else {
  cat("La transformacion logaritmica NO mejora el ajuste (BIC corregido).\n")
  cat("Interpretacion: una vez corregido el cambio de escala, la GBN sobre\n")
  cat("datos originales sigue siendo preferible. La mejora aparente en el\n")
  cat("BIC sin corregir se debia unicamente al cambio de unidades.\n")
}

# ------------------------------------------------------------
# 7. Guardar outputs
# ------------------------------------------------------------
saveRDS(tabla_np, "data/processed/scores_noparametrico.rds")
cat("\nGuardado: scores_noparametrico.rds\n")