# ============================================================
# Preprocesamiento CKD-IIJ para Redes Bayesianas Gaussianas
# Grupo 101 - Análisis de métodos de razonamiento e incertidumbre
# AD 2026 - Dr. Javier E. Garrido Guillén
# ============================================================
# Este script:
#   1. Carga los datos crudos (ckd.csv)
#   2. Limpia contaminación de texto (tabs, espacios, "?") que
#      corrompe columnas que deberían ser numéricas o categóricas
#   3. Corrige valores fisiológicamente imposibles en sod y pot
#   4. Construye el dataset con el universo continuo final
#      acordado (11 variables) para las 4 GBN, más classification
#      y dm de apoyo
#   5. Guarda un dataset limpio en data/processed/
# ============================================================

# ============================================================
# 01_preprocesamiento.R
# ============================================================
library(dplyr)

raw <- read.csv("data/raw/ckd.csv", stringsAsFactors = FALSE, na.strings = c("", "NA"))

clean_text <- function(x) {
  x <- trimws(x)
  x[x %in% c("?", "")] <- NA
  x
}

df <- raw %>% select(-id)

# Universo continuo final acordado por el equipo (se excluyen sg, al,
# su por no considerarse totalmente continuas; ver Métodos del artículo)
vars_continuas <- c("age", "bgr", "wc", "bp", "sc", "bu",
                    "pot", "sod", "hemo", "pcv", "rc")

# El csv trae columnas de texto contaminadas con tabs/espacios que
# rompen la lectura numérica de pcv, wc, rc y ensucian dm, cad,
# classification (valores como "\t43", "\t?", "\tno", " yes", "ckd\t")
es_texto <- sapply(df, is.character)
for (col in names(df)[es_texto]) {
  df[[col]] <- clean_text(df[[col]])
}

# Forzamos numérico en las continuas (incluye sg/al/su para dejarlas
# usables en el dataset completo, aunque no entren al set final)
for (col in c(vars_continuas, "sg", "al", "su")) {
  if (!is.numeric(df[[col]])) df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
}

df$classification <- clean_text(df$classification)
df$classification <- factor(df$classification, levels = c("notckd", "ckd"))

bin_cols <- c("rbc", "pc", "pcc", "ba", "htn", "dm", "cad", "appet", "pe", "ane")
for (col in bin_cols) { df[[col]] <- clean_text(df[[col]]) }

# Valores fisiológicamente imposibles (sod=4.5, pot=39/47 mEq/L no son
# compatibles con un paciente vivo; probable error de captura/decimal)
sospechosos <- df[which((df$pot > 10 & !is.na(df$pot)) |
                          (df$sod < 50 & !is.na(df$sod))),
                  c("age", "bp", "bu", "sc", "sod", "pot", "hemo", "classification")]
cat("Filas con sod/pot fisiológicamente implausibles (revisar manualmente):\n")
print(sospechosos)

df$sod[!is.na(df$sod) & df$sod < 50] <- NA
df$pot[!is.na(df$pot) & df$pot > 15] <- NA

# Diagnóstico de faltantes sobre el universo final, para citar en Métodos
faltantes <- sapply(df[vars_continuas], function(x) mean(is.na(x)) * 100)
cat("\n% de NA por variable (universo continuo final, 11 vars):\n")
print(round(sort(faltantes, decreasing = TRUE), 1))

cc <- sum(complete.cases(df[vars_continuas]))
cat(sprintf("\nCasos completos con las 11 variables juntas: %d de %d (%.0f%%)\n",
            cc, nrow(df), 100 * cc / nrow(df)))
cat("Usen la misma submuestra (mismas filas) para ajustar las 4 GBN,\n")
cat("para que BIC/AIC sean comparables entre estructuras (mismo n).\n")

# Guardar
df_final <- df %>% select(all_of(vars_continuas), classification, dm)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
saveRDS(df_final, "data/processed/ckd_clean.rds")
write.csv(df_final, "data/processed/ckd_clean.csv", row.names = FALSE, na = "")

# Version SOLO para ajuste de GBN: unicamente las 11 continuas, nada
# de classification/dm. 
df_gbn <- df %>% select(all_of(vars_continuas))
saveRDS(df_gbn, "data/processed/ckd_gbn.rds")
write.csv(df_gbn, "data/processed/ckd_gbn.csv", row.names = FALSE, na = "")


# También guardamos la versión completa (incluye sg, al, su, y las
# categóricas de texto) por si se usan en la sección de variables
# categóricas/discretas del artículo
# saveRDS(df, "data/processed/ckd_clean_full.rds")
# write.csv(df, "data/processed/ckd_clean_full.csv", row.names = FALSE, na = "")
