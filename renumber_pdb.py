import glob
import os
import pymol
from pymol import cmd

def procesar_pdb(filepath):
    filename = os.path.basename(filepath)
    object_name = "mol"

    # Cargar y modificar
    cmd.load(filepath, object_name)
    cmd.alter("chain B", "resi=str(int(resi)+1121)")
    cmd.alter("chain C", "resi=str(int(resi)+2242)")
    cmd.sort()  # ordena átomos si es necesario

    # Guardar y limpiar
    cmd.save(filepath, object_name)
    cmd.delete("all")

def main():
    pymol.finish_launching(['pymol', '-cq'])  # headless
    pdb_files = glob.glob("*_min.pdb")
    
    for pdb in pdb_files:
        print(f"Procesando: {pdb}")
        procesar_pdb(pdb)

if __name__ == "__main__":
    main()
