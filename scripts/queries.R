library(bnlearn)

# Load DAG from dag1.R
source("scripts/dag1.R")

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

# Query 3: ¿Cuál es la probabilidad de que el conteo de glóbulos blancos (wc)
# sea mayor a 11,000 cells/cumm, dado que el hematocrito (pcv) del paciente
# es menor a 30% (anemia)?
# P(wc > 11000 | pcv < 30)
cpquery(ckd_fit, event = (wc > 11000), evidence = list(pcv = 28), method = "lw")
#0.2159386

# Query 4: ¿Cuál es la probabilidad de que el potasio (pot) sea mayor a
# 5.5 mEq/L (hipercalemia), dado que la creatinina sérica (sc) es de
# 2.2 mg/dL y la urea en sangre (bu) es de 45 mg/dL?
# P(pot > 5.5 | sc = 2.2, bu = 45)
cpquery(ckd_fit, event = (pot > 5.5), evidence = list(sc = 2.2, bu = 45), method = "lw")
#0.06177862