import numpy as np
from itertools import combinations

# Parámetros
input_pdb = "SG_DSB.pdb" # this pdb only has the SG atoms, after extrating visually in pymol 
output_file = "puentes_disulfuro.txt"

# Paso 1: Leer los SG de residuos CYX
# Guardamos: (residue_number, coordenadas)
residuos = []
with open(input_pdb) as f:
    for line in f:
        if line.startswith("ATOM") and " SG " in line and "CYX" in line:
            resnum = int(line[22:26])  # <- Número del residuo
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            residuos.append((resnum, np.array([x, y, z])))

# Paso 2: Calcular distancias entre todos los pares posibles
pares = list(combinations(residuos, 2))
distancias = []
for (r1, p1), (r2, p2) in pares:
    dist = np.linalg.norm(p1 - p2)
    distancias.append((dist, r1, r2))

# Paso 3: Emparejar los más cercanos sin reutilizar
distancias.sort()
usados = set()
resultado = []

for dist, r1, r2 in distancias:
    if r1 not in usados and r2 not in usados:
        usados.add(r1)
        usados.add(r2)
        resultado.append((r1, r2, dist))

# Paso 4: Escribir resultado en archivo
with open(output_file, "w") as f:
    for r1, r2, d in resultado:
        f.write(f"{r1}\t{r2}\t{d:.3f}\n")

print(f"Archivo guardado como {output_file}")
