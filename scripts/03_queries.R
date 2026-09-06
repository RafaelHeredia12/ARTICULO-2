# ============================================================
# 03_queries.R
# Inferencia con la mejor GBN (DAG 1)
# Proyecto: Redes Bayesianas Gaussianas - ERC
# ============================================================
# Carga la mejor red ya ajustada por 02_comparacion_dags.R y
# responde las 6 consultas propuestas por las especialistas.
#
# Estructura DAG 1 (referencia para el analisis de d-separacion):
#   age  -> bp,  age -> bgr
#   bp   -> sc,  bgr -> sc
#   sc   -> bu, sod, pot, hemo, wc
#   hemo -> pcv, rc
#
# ============================================================

library(bnlearn)

# ------------------------------------------------------------
# 1. Cargar la mejor red ajustada
# ------------------------------------------------------------
ckd_fit <- readRDS("data/processed/mejor_gbn.rds")

set.seed(123)   # cpquery con method="lw" es estocastico
N <- 1e6        # numero de muestras por consulta

# ------------------------------------------------------------
# Q1 - Especialista 1
# ¿Probabilidad de creatinina serica > 1.5 mg/dL dado que la
# presion arterial es 100 mmHg y la glucosa 150 mg/dL?
# Camino: bp -> sc <- bgr (padres directos)
# ------------------------------------------------------------
q1 <- cpquery(ckd_fit, event = (sc > 1.5),
              evidence = list(bp = 100, bgr = 150),
              method = "lw", n = N)
cat("Q1  P(sc > 1.5 | bp = 100, bgr = 150) =", round(q1, 4), "\n")

# ------------------------------------------------------------
# Q2 - Especialista 1
# ¿Probabilidad de hemoglobina < 10 g/dL dado 65 años?
# Camino: age -> bp -> sc -> hemo  y  age -> bgr -> sc -> hemo
# ------------------------------------------------------------
q2 <- cpquery(ckd_fit, event = (hemo < 10),
              evidence = list(age = 65),
              method = "lw", n = N)
cat("Q2  P(hemo < 10 | age = 65) =", round(q2, 4), "\n")

# ------------------------------------------------------------
# Q3 - Especialista 2
# ¿Probabilidad de leucocitos > 11,000 dado hematocrito bajo?
# Camino ACTIVO: wc <- sc -> hemo -> pcv
# El enunciado clinico dice "pcv < 30"; cpquery con method="lw"
# requiere evidencia puntual para nodos continuos, por lo que se
# usa pcv = 28 como valor representativo de ese rango. Esta
# aproximacion debe declararse en el articulo.
# ------------------------------------------------------------
q3 <- cpquery(ckd_fit, event = (wc > 11000),
              evidence = list(pcv = 28),
              method = "lw", n = N)
cat("Q3  P(wc > 11000 | pcv = 28) =", round(q3, 4), "\n")

# ------------------------------------------------------------
# Q4 - Especialista 2
# ¿Probabilidad de potasio > 5.5 mEq/L dado sc = 2.2 y bu = 45?
# sc -> pot es el camino informativo. La evidencia en bu es
# REDUNDANTE dado sc (ambos hijos de sc: condicionar en el padre
# comun los d-separa). Se calculan ambas versiones para
# documentar empiricamente esa redundancia.
# ------------------------------------------------------------
q4 <- cpquery(ckd_fit, event = (pot > 5.5),
              evidence = list(sc = 2.2, bu = 45),
              method = "lw", n = N)
q4_solo_sc <- cpquery(ckd_fit, event = (pot > 5.5),
                      evidence = list(sc = 2.2),
                      method = "lw", n = N)
cat("Q4  P(pot > 5.5 | sc = 2.2, bu = 45) =", round(q4, 4), "\n")
cat("    P(pot > 5.5 | sc = 2.2)          =", round(q4_solo_sc, 4),
    "  <- practicamente identica: bu es redundante dado sc\n")

# ------------------------------------------------------------
# Q5 - Especialista 3
# ¿Probabilidad de hematocrito < 33% dado bu, sod y wc?
# Caminos ACTIVOS: cada uno de bu, sod y wc se conecta con pcv
# a traves de  X <- sc -> hemo -> pcv
# ------------------------------------------------------------
q5 <- cpquery(ckd_fit, event = (pcv < 33),
              evidence = list(bu = 60, sod = 130, wc = 11000),
              method = "lw", n = N)
cat("Q5  P(pcv < 33 | bu = 60, sod = 130, wc = 11000) =", round(q5, 4), "\n")

# ------------------------------------------------------------
# Q6 - Especialista 3
# ¿Probabilidad de globulos rojos < 3.8 dado hemoglobina 9.5?
# Camino: hemo -> rc (padre directo)
# ------------------------------------------------------------
q6 <- cpquery(ckd_fit, event = (rc < 3.8),
              evidence = list(hemo = 9.5),
              method = "lw", n = N)
cat("Q6  P(rc < 3.8 | hemo = 9.5) =", round(q6, 4), "\n")

# ------------------------------------------------------------
# 2. Probabilidades marginales de referencia
#    Permiten cuantificar cuanto desplaza la evidencia a la
#    probabilidad del evento: si condicional ~ marginal, la
#    evidencia aporta poca informacion aunque el camino
#    este activo.
# ------------------------------------------------------------
cat("\n=== Marginales de referencia (sin evidencia) ===\n")
m_sc   <- cpquery(ckd_fit, event = (sc > 1.5),    evidence = TRUE, n = N)
m_hemo <- cpquery(ckd_fit, event = (hemo < 10),   evidence = TRUE, n = N)
m_wc   <- cpquery(ckd_fit, event = (wc > 11000),  evidence = TRUE, n = N)
m_pot  <- cpquery(ckd_fit, event = (pot > 5.5),   evidence = TRUE, n = N)
m_pcv  <- cpquery(ckd_fit, event = (pcv < 33),    evidence = TRUE, n = N)
m_rc   <- cpquery(ckd_fit, event = (rc < 3.8),    evidence = TRUE, n = N)

marginales <- c(m_sc, m_hemo, m_wc, m_pot, m_pcv, m_rc)
condicionales <- c(q1, q2, q3, q4, q5, q6)

comparacion <- data.frame(
  Consulta   = paste0("Q", 1:6),
  Marginal   = round(marginales, 4),
  Condicional= round(condicionales, 4),
  Desplaza   = round(condicionales - marginales, 4)
)
print(comparacion, row.names = FALSE)

# ------------------------------------------------------------
# 3. Guardar resultados para el articulo
# ------------------------------------------------------------
resultados <- list(
  q1 = list(enunciado = "P(sc > 1.5 | bp = 100, bgr = 150)",
            resultado = q1, marginal = m_sc),
  q2 = list(enunciado = "P(hemo < 10 | age = 65)",
            resultado = q2, marginal = m_hemo),
  q3 = list(enunciado = "P(wc > 11000 | pcv = 28)",
            resultado = q3, marginal = m_wc,
            nota = "enunciado original pcv < 30; se usa 28 como valor puntual"),
  q4 = list(enunciado = "P(pot > 5.5 | sc = 2.2, bu = 45)",
            resultado = q4, marginal = m_pot,
            nota = "bu redundante dado sc; sin bu: P = ", q4_solo_sc),
  q5 = list(enunciado = "P(pcv < 33 | bu = 60, sod = 130, wc = 11000)",
            resultado = q5, marginal = m_pcv),
  q6 = list(enunciado = "P(rc < 3.8 | hemo = 9.5)",
            resultado = q6, marginal = m_rc)
)

saveRDS(resultados,  "data/processed/resultados_queries.rds")
saveRDS(comparacion, "data/processed/comparacion_marginales.rds")
cat("\nGuardado: resultados_queries.rds, comparacion_marginales.rds\n")