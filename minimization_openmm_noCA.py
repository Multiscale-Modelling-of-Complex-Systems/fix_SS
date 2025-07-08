import os
import glob
import time
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

# === Selección de GPU si está disponible ===
def seleccionar_gpu():
    for nombre in ['CUDA', 'OpenCL']:
        try:
            platform = Platform.getPlatformByName(nombre)
            print(f"Usando plataforma GPU: {platform.getName()}")
            return platform
        except Exception:
            continue
    print("No se encontró GPU, usando CPU.")
    return Platform.getPlatformByName("CPU")

# === Extraer enlaces CONECT (índices de átomos) desde archivo PDB ===
def extraer_conect_pairs(pdb_filename):
    conect_pairs = []
    with open(pdb_filename, 'r') as f:
        for line in f:
            if line.startswith("CONECT"):
                partes = line.split()
                if len(partes) >= 3:
                    a1 = int(partes[1])
                    a2 = int(partes[2])
                    conect_pairs.append((a1 - 1, a2 - 1))  # Convertimos a 0-based
    return conect_pairs

# === Proceso principal ===
gpu_platform = seleccionar_gpu()
pdb_files = glob.glob("*_conectado.pdb")

for pdb_file in pdb_files:
    print(f"\nMinimizando: {pdb_file}")

    conect_pairs = extraer_conect_pairs(pdb_file)

    # Cargar PDB en OpenMM
    pdb = PDBFile(pdb_file)
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

    system = forcefield.createSystem(pdb.topology,
                                     nonbondedMethod=NoCutoff,
                                     constraints=HBonds)

    # Agregar fuerza de enlace para disulfuros
    bond_length = 0.2 * nanometers
    bond_k = 500 * kilojoules_per_mole / nanometers**2
    disulfide_bonds = HarmonicBondForce()
    for a1, a2 in conect_pairs:
        disulfide_bonds.addBond(a1, a2, bond_length, bond_k)
    system.addForce(disulfide_bonds)

    integrator = LangevinIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
    simulation = Simulation(pdb.topology, system, integrator, platform=gpu_platform)
    simulation.context.setPositions(pdb.positions)

    # Imprimir energía inicial
    energia_inicial = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    print(f"Energía inicial: {energia_inicial}")

    # Minimización
    start = time.time()
    simulation.minimizeEnergy(maxIterations=50000)
    end = time.time()

    energia_final = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    print(f"Energía final: {energia_final}")
    print(f"Duración: {end - start:.2f} s")

    # Guardar nuevo archivo minimizado
    output_file = pdb_file.replace("_conectado.pdb", "_min.pdb")
    with open(output_file, 'w') as f:
        positions = simulation.context.getState(getPositions=True).getPositions()
        PDBFile.writeFile(simulation.topology, positions, f)

    print(f"Archivo guardado: {output_file}")
