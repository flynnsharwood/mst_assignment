# Starter kit: graphene slit, CH4/CO2(/H2O)

Templates only - nothing here was run through GROMACS, Packmol or LAMMPS. Check every file before use.

common/build_graphene.py   AB bilayer graphene, 25x15 cells (6.15 x 6.39 nm, 3000 C):
                           writes graphene.gro/.itp (GROMACS) and graphene_lammps.data (LAMMPS, slit H = 20 A)
common/pr_fugacity.py      Peng-Robinson densities and fugacity coefficients for the 50/50 mixture
partA/setup_partA.sh       runs Part A for one system + temperature (bash setup_partA.sh slit 275)
partA/                     GROMACS: topology, molecule itps/pdbs, Packmol inputs, em/nvt/npt/prod mdp files
partB/                     LAMMPS: molecule templates, graphene data file, hybrid GCMC/MD input

Fill in before running:
- (done) CO2 bending constant and graphene LJ now taken from the paper SI, Tables S1 and S3
- ref-t in every mdp for 295 and 315 K runs; gen-temp in nvt.mdp
- partB: -var T, Patm, phiCH4, phiCO2 from the table in the plan
# mst_assignment
# mst_assignment
