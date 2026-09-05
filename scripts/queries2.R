library(bnlearn)

# 1. Cargar la DAG fisiopatológica redefinida
source("scripts/dag_fisiopatologica.R")

# 2. Cargar y preparar los datos
df <- read.csv("data/processed/ckd_gbn.csv")
datos_bn <- na.omit(df)
datos_bn[] <- lapply(datos_bn, as.numeric)

# 3. Ajustar la red bayesiana gaussiana
ckd_fit <- bn.fit(dag, datos_bn)


# Query 1: ¿Cuál es la probabilidad de presentar presión arterial elevada
# (bp > 140 mmHg), dado que la glucosa en sangre (bgr) es de 180 mg/dL
# y el paciente tiene 60 años de edad (age)?
# P(bp > 140 | bgr = 180, age = 60)
cpquery(ckd_fit, event = (bp > 140), evidence = list(bgr = 180, age = 60), method = "lw")

# Query 2: ¿Cuál es la probabilidad de registrar valores altos de urea en sangre
# (bu > 50 mg/dL), dado que la creatinina sérica (sc) está en 2.5 mg/dL?
# P(bu > 50 | sc = 2.5)
cpquery(ckd_fit, event = (bu > 50), evidence = list(sc = 2.5), method = "lw")

# Query 3: ¿Cuál es la probabilidad de tener leucocitosis o glóbulos blancos elevados
# (wc > 10000 cells/cumm), dado que la glucosa en sangre (bgr) es de 200 mg/dL?
# P(wc > 10000 | bgr = 200)
cpquery(ckd_fit, event = (wc > 10000), evidence = list(bgr = 200), method = "lw")

# Query 4: ¿Cuál es la probabilidad de que el hematocrito (pcv) sea menor
# a 33%, dado un perfil con urea (bu) de 60 mg/dL, sodio (sod) de 130 mEq/L
# y glóbulos blancos (wc) de 11,000 cells/cumm?
# P(pcv < 33 | bu = 60, sod = 130, wc = 11000)
cpquery(ckd_fit, event = (pcv < 33), evidence = list(bu = 60, sod = 130, wc = 11000), method = "lw")

# Query 5: ¿Cuál es la probabilidad de tener un conteo bajo de glóbulos rojos
# (rc < 3.8 M/mcL), dado que la hemoglobina (hemo) ha caído a 9.5 g/dL?
# P(rc < 3.8 | hemo = 9.5)
cpquery(ckd_fit, event = (rc < 3.8), evidence = list(hemo = 9.5), method = "lw")