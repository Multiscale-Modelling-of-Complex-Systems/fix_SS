import os
import glob
from Bio.PDB import PDBParser

# === Leer pares disulfuro desde archivo ===
def leer_pares_disulfuro(nombre_archivo):
    pares = []
    with open(nombre_archivo, 'r') as f:
        for linea in f:
            if linea.strip():
                partes = linea.split()
                pares.append((int(partes[0]), int(partes[1])))
    return pares

# === Buscar serial_number de átomos SG para cada residuo CYS ===
def encontrar_sg_seriales(pdb_file, pares_disulfuro):
    parser = PDBParser(QUIET=True)
    estructura = parser.get_structure("modelo", pdb_file)
    sg_seriales = {}  # resid -> serial_number del átomo SG

    for model in estructura:
        for chain in model:
            for residue in chain:
                if residue.resname == 'CYS':
                    resid = residue.id[1]
                    for atom in residue:
                        if atom.name == 'SG':
                            sg_seriales[(chain.id, resid)] = atom.serial_number

    conect_pairs = []
    for r1, r2 in pares_disulfuro:
        # Buscar SG sin importar cadena
        sg1 = next((v for (cid, resid), v in sg_seriales.items() if resid == r1), None)
        sg2 = next((v for (cid, resid), v in sg_seriales.items() if resid == r2), None)

        if sg1 is not None and sg2 is not None:
            conect_pairs.append((sg1, sg2))
        else:
            print(f"WARNING: SG no encontrado para residuo {r1} o {r2} en {pdb_file}")

    return conect_pairs

# === Generar líneas CONECT ===
def generar_conect_lines(pares_seriales):
    lines = []
    for s1, s2 in pares_seriales:
        lines.append(f"CONECT{s1:>5}{s2:>5}\n")
    return lines

# === Agregar CONECT al PDB ===
def agregar_conect_al_pdb(pdb_file, conect_lines, output_file):
    with open(pdb_file, 'r') as original, open(output_file, 'w') as nuevo:
        lines = original.readlines()
        # Elimina líneas CONECT o END antiguas
        lines = [l for l in lines if not l.startswith('CONECT') and not l.startswith('END')]
        nuevo.writelines(lines)
        nuevo.write("".join(conect_lines))
        nuevo.write("END\n")

# === Proceso principal ===
pares_disulfuro = leer_pares_disulfuro("puentes_disulfuro.txt")
pdb_files = glob.glob("*noHG1.pdb")

for pdb_file in pdb_files:
    print(f"\nProcesando: {pdb_file}")
    conect_pairs = encontrar_sg_seriales(pdb_file, pares_disulfuro)
    conect_lines = generar_conect_lines(conect_pairs)
    output_pdb = pdb_file.replace(".pdb", "_conectado.pdb")
    agregar_conect_al_pdb(pdb_file, conect_lines, output_pdb)
    print(f"Archivo generado: {output_pdb}")

