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
# Los scores se convierten a la escala de bnlearn
#     BIC_bnlearn = -BIC_R / 2
# para que sean directamente comparables con los reportados en
# los scripts 02 y 04 (donde MAYOR es mejor).
#
# Output:
#   output/figures/efectos_parciales_gam.png
#   data/processed/scores_noparametrico.rds
#   data/processed/edf_gam.rds
#   data/processed/bic_por_nodo.rds
#   data/processed/cv_gam.rds
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
#
#    Las funciones BIC()/AIC() de R usan la convencion
#        BIC_R = -2*loglik + k*log(n)      (MENOR es mejor)
#    mientras que score() de bnlearn usa
#        BIC_bnlearn = loglik - (k/2)*log(n)   (MAYOR es mejor)
#    Ambas se relacionan por  BIC_bnlearn = -BIC_R / 2.
#    Convertimos para reportar todo en la escala de bnlearn, que es
#    la misma usada en las tablas de los scripts 02 y 04.
# ------------------------------------------------------------
bic_lineal_R <- sum(sapply(modelos_lin, BIC))
aic_lineal_R <- sum(sapply(modelos_lin, AIC))
bic_gam_R    <- sum(sapply(modelos_gam, BIC))
aic_gam_R    <- sum(sapply(modelos_gam, AIC))

# Conversion a la escala de bnlearn
bic_lineal <- -bic_lineal_R / 2
aic_lineal <- -aic_lineal_R / 2
bic_gam    <- -bic_gam_R / 2
aic_gam    <- -aic_gam_R / 2

tabla_np <- data.frame(
  Modelo = c("GBN lineal (DAG 1)",
             "Modelo no parametrico GAM/splines (DAG 1)"),
  BIC    = round(c(bic_lineal, bic_gam), 2),
  AIC    = round(c(aic_lineal, aic_gam), 2)
)

cat("\n=== Comparacion: GBN lineal vs. GAM no parametrico ===\n")
cat("(escala bnlearn: MAYOR es mejor, comparable con scripts 02 y 04)\n\n")
print(tabla_np)

# Verificacion: el BIC lineal calculado aqui debe coincidir con el
# que reporto score(dag1, dat, type = "bic-g") en el script 02.
dag1 <- empty.graph(names(dat))
arcs(dag1) <- matrix(c(
  "age", "bp",  "age", "bgr",
  "bp",  "sc",  "bgr", "sc",
  "sc",  "bu",  "sc",  "sod",
  "sc",  "pot", "sc",  "hemo",
  "sc",  "wc",  "hemo","pcv",
  "hemo","rc"
), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from","to")))

bic_bnlearn <- score(dag1, data = dat, type = "bic-g")
cat("\nVerificacion de la conversion:\n")
cat("  BIC lineal via BIC() convertido :", round(bic_lineal, 2), "\n")
cat("  BIC lineal via score() bnlearn  :", round(bic_bnlearn, 2), "\n")
cat("  Diferencia                      :", round(abs(bic_lineal - bic_bnlearn), 4), "\n\n")

if (bic_gam > bic_lineal) {
  cat("El modelo no parametrico MEJORA el BIC.\n")
  cat("Diferencia:", round(bic_gam - bic_lineal, 2), "puntos.\n")
} else {
  cat("El modelo no parametrico NO mejora el BIC.\n")
  cat("Diferencia:", round(bic_lineal - bic_gam, 2), "puntos a favor del lineal.\n")
}

if (aic_gam > aic_lineal) {
  cat("El modelo no parametrico MEJORA el AIC.\n")
  cat("Diferencia:", round(aic_gam - aic_lineal, 2), "puntos.\n")
} else {
  cat("El modelo no parametrico NO mejora el AIC.\n")
  cat("Diferencia:", round(aic_lineal - aic_gam, 2), "puntos a favor del lineal.\n")
}

# ------------------------------------------------------------
# 5. Descomposicion de la mejora por nodo
#    La mejora agregada puede ocultar que el spline ayuda en unos
#    nodos y estorba en otros. Descomponemos BIC por distribucion
#    local para ver de donde viene realmente la ganancia.
# ------------------------------------------------------------
filas_nodo <- list()
for (nodo in names(padres)) {
  if (length(padres[[nodo]]) == 0) next   # nodo raiz: mismo modelo
  bic_l <- -BIC(modelos_lin[[nodo]]) / 2
  bic_g <- -BIC(modelos_gam[[nodo]]) / 2
  filas_nodo[[length(filas_nodo) + 1]] <- data.frame(
    Nodo       = nodo,
    BIC_lineal = round(bic_l, 2),
    BIC_GAM    = round(bic_g, 2),
    Mejora     = round(bic_g - bic_l, 2)
  )
}
tabla_nodo <- do.call(rbind, filas_nodo)
tabla_nodo <- tabla_nodo[order(-tabla_nodo$Mejora), ]

cat("\n=== Mejora del BIC por nodo (escala bnlearn) ===\n")
cat("Mejora > 0 => el spline ayuda en ese nodo\n\n")
print(tabla_nodo, row.names = FALSE)

n_mejora <- sum(tabla_nodo$Mejora > 0)
cat("\nNodos donde el spline mejora el BIC:", n_mejora, "de", nrow(tabla_nodo), "\n")
cat("La mejora agregada NO esta repartida de forma uniforme.\n")

# ------------------------------------------------------------
# 6. Validacion cruzada 10-fold
#    El BIC penaliza complejidad pero se calcula dentro de muestra.
#    Comprobamos si la flexibilidad del spline se traduce en mejor
#    prediccion fuera de muestra o si es sobreajuste.
# ------------------------------------------------------------
set.seed(42)
n_obs <- nrow(dat)
folds <- sample(rep(1:10, length.out = n_obs))

filas_cv <- list()
for (nodo in names(padres)) {
  pa <- padres[[nodo]]
  if (length(pa) == 0) next
  
  f_lin <- as.formula(paste(nodo, "~", paste(pa, collapse = " + ")))
  err_lin <- err_gam <- numeric(0)
  
  for (k in 1:10) {
    entrena <- dat[folds != k, ]
    prueba  <- dat[folds == k, ]
    
    # k del spline recalculado dentro del fold de entrenamiento
    ts_cv <- function(v) {
      kk <- min(10, length(unique(entrena[[v]])) - 1)
      paste0("s(", v, ", k = ", kk, ")")
    }
    f_gam_cv <- as.formula(paste(nodo, "~",
                                 paste(sapply(pa, ts_cv), collapse = " + ")))
    
    m_lin <- lm(f_lin, data = entrena)
    m_gam <- suppressWarnings(gam(f_gam_cv, data = entrena, method = "REML"))
    
    err_lin <- c(err_lin, prueba[[nodo]] - predict(m_lin, prueba))
    err_gam <- c(err_gam, prueba[[nodo]] - predict(m_gam, prueba))
  }
  
  rmse_lin <- sqrt(mean(err_lin^2))
  rmse_gam <- sqrt(mean(err_gam^2))
  filas_cv[[length(filas_cv) + 1]] <- data.frame(
    Nodo       = nodo,
    RMSE_lin   = round(rmse_lin, 3),
    RMSE_GAM   = round(rmse_gam, 3),
    Mejora_pct = round(100 * (rmse_lin - rmse_gam) / rmse_lin, 1)
  )
}
tabla_cv <- do.call(rbind, filas_cv)
tabla_cv <- tabla_cv[order(-tabla_cv$Mejora_pct), ]

cat("\n=== Validacion cruzada 10-fold (error fuera de muestra) ===\n")
cat("Mejora_pct > 0 => el spline predice mejor en datos no vistos\n\n")
print(tabla_cv, row.names = FALSE)
cat("\nNodos donde el spline predice mejor:",
    sum(tabla_cv$Mejora_pct > 0), "de", nrow(tabla_cv), "\n")

# ------------------------------------------------------------
# 7. Grados de libertad efectivos (edf) por nodo
#    edf ~ 1 indica relacion practicamente lineal;
#    edf > 2 sugiere no linealidad relevante.
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
# 8. Figura: efectos parciales de los splines
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
# 9. Guardar outputs
# ------------------------------------------------------------
saveRDS(tabla_np,   "data/processed/scores_noparametrico.rds")
saveRDS(tabla_edf,  "data/processed/edf_gam.rds")
saveRDS(tabla_nodo, "data/processed/bic_por_nodo.rds")
saveRDS(tabla_cv,   "data/processed/cv_gam.rds")
cat("Guardado: scores_noparametrico.rds, edf_gam.rds,\n")
cat("          bic_por_nodo.rds, cv_gam.rds\n")