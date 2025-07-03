"""
Script para eliminar átomos HG1 de residuos listados en un archivo TXT dentro de un archivo PDB.

Uso:
    python eliminar_HG1.py --input archivo.pdb --list residuos.txt --output archivo_modificado.pdb

Flags:
    --input   : Ruta al archivo PDB de entrada.
    --list    : Ruta al archivo TXT con la lista de residuos (dos columnas separadas por tabulador).
    --output  : Ruta donde se guardará el archivo PDB modificado.

Ejemplo:
    python eliminar_HG1.py --input mi_estructura.pdb --list residuos.txt --output mi_estructura_mod.pdb
"""

import argparse
from Bio.PDB import PDBParser, PDBIO

def cargar_residuos_a_eliminar(ruta_txt):
    residuos = set()
    with open(ruta_txt, 'r') as f:
        for linea in f:
            if linea.strip():
                partes = linea.strip().split('\t')
                residuos.add(int(partes[0]))
                residuos.add(int(partes[1]))
    return residuos

def eliminar_HG1_de_residuos(pdb_input, pdb_output, lista_residuos_txt):
    parser = PDBParser(QUIET=True)
    estructura = parser.get_structure('estructura', pdb_input)

    residuos_a_eliminar = cargar_residuos_a_eliminar(lista_residuos_txt)

    for modelo in estructura:
        for cadena in modelo:
            for residuo in list(cadena):
                res_id = residuo.get_id()
                resseq = res_id[1]

                if resseq in residuos_a_eliminar:
                    for atomo in list(residuo):
                        if atomo.get_id() == 'HG1':
                            residuo.detach_child(atomo.get_id())

    io = PDBIO()
    io.set_structure(estructura)
    io.save(pdb_output)

def main():
    parser = argparse.ArgumentParser(description="Eliminar átomos HG1 de residuos específicos en un archivo PDB")
    parser.add_argument('--input', required=True, help="Archivo PDB de entrada")
    parser.add_argument('--list', required=True, help="Archivo TXT con lista de residuos (dos columnas tabuladas)")
    parser.add_argument('--output', required=True, help="Archivo PDB de salida modificado")

    args = parser.parse_args()

    eliminar_HG1_de_residuos(args.input, args.output, args.list)

if __name__ == '__main__':
    main()

