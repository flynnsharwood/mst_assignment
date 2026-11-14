#!/bin/bash
# Part A pipeline for ONE system at ONE temperature.  RESUMABLE.
# Usage:  bash setup_partA.sh slit 275     or     bash setup_partA.sh bulk 315
# Stop any time with Ctrl+C (once). Run the same command again to resume:
# finished stages are skipped, an interrupted stage continues from its checkpoint.
set -e
export HWLOC_COMPONENTS=-gl,-opencl     # avoids the hwloc X-display hang on this machine
SYS=$1; T=$2; NT=${NT:-16}
D=${SYS}_${T}; TOP=topol_${SYS}.top
[ "$SYS" = slit ] && NPT=npt_slit.mdp || NPT=npt_bulk.mdp

if [ ! -d $D ]; then
  mkdir -p $D; cp *.itp *.pdb *.mdp $TOP pack_${SYS}.inp graphene.gro $D/
  for f in $D/*.mdp; do sed -i -e "s/^ref-t .*/ref-t           = $T/" -e "s/^gen-temp .*/gen-temp        = $T/" $f; done
  echo "seed $RANDOM" >> $D/pack_${SYS}.inp
fi
cd $D

# run_stage NAME MDP COORD [extra grompp args]
run_stage () {
  local N=$1 MDP=$2 C=$3; shift 3
  if [ -f $N.done ]; then echo ">>> $N already finished, skipping"; return; fi
  [ -f $N.tpr ] || gmx grompp -f $MDP -c $C -r $C -p $TOP -n index.ndx -o $N.tpr -maxwarn 1 "$@"
  gmx mdrun -deffnm $N -cpi $N.cpt -ntmpi 1 -ntomp $NT || true
  # complete = no stop signal after the last start of mdrun in the log, and a final .gro exists
  local s=$(grep -n "Started mdrun" $N.log | tail -1 | cut -d: -f1)
  local k=$(grep -n "Received the .* signal" $N.log | tail -1 | cut -d: -f1)
  if [ ! -f $N.gro ] || { [ -n "$k" ] && [ "$k" -gt "${s:-0}" ]; }; then
    echo ">>> $N was interrupted. Run the same command again to resume from $N.cpt."; exit 1
  fi
  touch $N.done; echo ">>> $N finished"
}

# 1-2. build coordinates and index once
if [ ! -f start.gro ]; then
  packmol < pack_${SYS}.inp
  if [ "$SYS" = slit ]; then
    gmx editconf -f fluid_slit.pdb -o fluid.gro -box 6.14878 6.39000 8.00000
    n1=$(sed -n 2p graphene.gro); n2=$(sed -n 2p fluid.gro)
    { echo "graphene slit"; echo $((n1+n2)); sed '1,2d;$d' graphene.gro; sed '1,2d;$d' fluid.gro
      echo "   6.14878   6.39000   8.00000"; } > start.tmp && mv start.tmp start.gro
  else
    gmx editconf -f fluid_bulk.pdb -o start.gro -box 6.5 6.5 6.5
    sed -i 's/^compressed-x-grps.*//' prod.mdp
  fi
fi
if [ ! -f index.ndx ]; then
  printf 'q\n' | gmx make_ndx -f start.gro -o index_tmp.ndx
  if [ "$SYS" = slit ]; then
    gmx select -s start.gro -select '"Fluid" not resname GRA' -on fluid.ndx; cat fluid.ndx >> index_tmp.ndx
  fi
  mv index_tmp.ndx index.ndx
fi

# 3. minimise, heat, equilibrate density
if [ ! -f em.done ]; then      # EM has no checkpoint; it just reruns (about a minute)
  gmx grompp -f em.mdp -c start.gro -r start.gro -p $TOP -n index.ndx -o em.tpr -maxwarn 1
  gmx mdrun -deffnm em -ntmpi 1 -ntomp $NT
  [ -f em.gro ] && touch em.done
fi
run_stage nvt nvt.mdp em.gro
run_stage npt $NPT nvt.gro -t nvt.cpt

# 4. box = mean of last 2 ns of NPT (velocities kept), then production
if [ ! -f prod_start.gro ]; then
  [ "$SYS" = slit ] && TERM=Box-Z || TERM=Volume
  echo $TERM | gmx energy -f npt.edr -b 3000 -o box.xvg > box.log 2>&1
  MEAN=$(grep "^$TERM" box.log | tail -1 | awk '{print $2}'); echo ">>> mean $TERM: $MEAN"
  if [ "$SYS" = slit ]; then
    LZ=$(tail -1 npt.gro | awk '{print $3}'); S=$(python3 -c "print($MEAN/$LZ)")
    gmx editconf -f npt.gro -o prod_start.tmp.gro -scale 1 1 $S
  else
    V=$(tail -1 npt.gro | awk '{print $1*$2*$3}'); S=$(python3 -c "print(($MEAN/$V)**(1/3))")
    gmx editconf -f npt.gro -o prod_start.tmp.gro -scale $S $S $S
  fi
  mv prod_start.tmp.gro prod_start.gro
fi
run_stage prod prod.mdp prod_start.gro
echo ">>> all stages done: $D"
