library(bnlearn)

df <- read.csv("data/processed/ckd_gbn.csv")
datos_bn <- na.omit(df)
datos_bn[] <- lapply(datos_bn, as.numeric)

# Ajustar la red gaussiana con los datos ya preparados (dag1.R)
ckd_fit <- bn.fit(dag, datos_bn)

# Query 1: ¿Cuál es la probabilidad de que la creatinina sérica (sc) sea
# mayor a 1.5 mg/dL, dado que la presión arterial (bp) es de 100 mmHg y
# la glucosa en sangre (bgr) es de 150 mg/dL?
# P(sc > 1.5 | bp = 100, bgr = 150)
cpquery(ckd_fit, event = (sc > 1.5), evidence = list(bp = 100, bgr = 150), method = "lw")
# 0.7733688 

# Query 2: ¿Cuál es la probabilidad de que la hemoglobina (hemo) sea menor
# a 10 g/dL (anemia), dado que el paciente tiene 65 años de edad (age)?
# P(hemo < 10 | age = 65)
cpquery(ckd_fit, event = (hemo < 10), evidence = list(age = 65), method = "lw")
#0.1414545


