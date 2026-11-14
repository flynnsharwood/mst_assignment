#!/usr/bin/env python3
"""Build an AB-stacked bilayer graphene wall that is periodic in x and y.
nx x ny rectangular 4-atom cells (a = 0.246 nm, b = 0.426 nm).
Default 25 x 15 cells -> 6.15 x 6.39 nm, 1500 C per layer.

Outputs
  graphene.gro / graphene.itp     (GROMACS; neutral C, Steele LJ, position restraints)
  graphene_lammps.data            (LAMMPS, atom_style full, with empty slit of width H)
Usage: python3 build_graphene.py [nx ny box_z_nm]
"""
import sys
acc = 0.142                    # C-C bond, nm
a = 3**0.5*acc                 # 0.2460 nm
b = 3*acc                      # 0.4260 nm
nx, ny = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (25, 15)
boxz = float(sys.argv[3]) if len(sys.argv) > 3 else 7.0   # nm; NPT will fix it
dz = 0.335                     # interlayer spacing, nm
basis = [(0, 0), (0, acc), (a/2, 1.5*acc), (a/2, 2.5*acc)]
Lx, Ly = nx*a, ny*b
atoms = []
for layer in range(2):
    shift = 0.0 if layer == 0 else acc        # AB stacking: shift along y by one bond
    z = 0.05 + layer*dz
    for i in range(nx):
        for j in range(ny):
            for bx, by in basis:
                atoms.append((round(i*a+bx, 5), round((j*b+by+shift) % Ly, 5), z))
n = len(atoms)

with open("graphene.gro", "w") as f:
    f.write(f"AB bilayer graphene {nx}x{ny} cells\n{n}\n")
    for k, (x, y, z) in enumerate(atoms, 1):
        f.write(f"{1:5d}{'GRA':<5s}{'C':>5s}{k%100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
    f.write(f"{Lx:10.5f}{Ly:10.5f}{boxz:10.5f}\n")

with open("graphene.itp", "w") as f:
    f.write("; Neutral graphene carbon, Steele LJ: sigma 0.340 nm, eps 0.2328 kJ/mol\n")
    f.write("; Put this line in [ atomtypes ] of your top file:\n")
    f.write(";  CG  CG  6  12.011  0.000  A  0.340  0.2328\n\n")
    f.write("[ moleculetype ]\nGRA  1\n\n[ atoms ]\n")
    for k in range(1, n+1):
        f.write(f"{k:6d}  CG  1  GRA  C  {k:6d}  0.000  12.011\n")
    f.write("\n#ifndef NO_GRA_POSRES\n[ position_restraints ]\n")
    for k in range(1, n+1):
        f.write(f"{k:6d}  1  100000 100000 100000\n")
    f.write("#endif\n")

# LAMMPS data (Angstrom). Types: 1 C_gra, 2 C_ch4, 3 H_ch4, 4 C_co2, 5 O_co2
H = 20.0   # slit width wall-to-wall in Angstrom for Part B (edit as needed)
zhi = 0.5 + dz*10 + H      # wall at z~0.5-3.85 A, its periodic image closes the slit
with open("graphene_lammps.data", "w") as f:
    f.write(f"Bilayer graphene wall, slit H = {H} A (Part B)\n\n{n} atoms\n5 atom types\n\n")
    f.write(f"0.0 {Lx*10:.5f} xlo xhi\n0.0 {Ly*10:.5f} ylo yhi\n0.0 {zhi:.5f} zlo zhi\n\n")
    f.write("Masses\n\n1 12.011\n2 12.011\n3 1.008\n4 12.011\n5 15.999\n\nAtoms # full\n\n")
    for k, (x, y, z) in enumerate(atoms, 1):
        f.write(f"{k} 1 1 0.0 {x*10:.4f} {y*10:.4f} {z*10:.4f}\n")
print(f"{n} C atoms, box {Lx:.4f} x {Ly:.4f} nm, LAMMPS zhi = {zhi:.2f} A")
