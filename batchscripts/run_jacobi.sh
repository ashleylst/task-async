#!/bin/sh
#SBATCH --time=00:40:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192

export NX=8192
export NY=8192

export OMP_DISPLAY_ENV=true
export OMP_DISPLAY_AFFINITY=true
export OMP_AFFINITY_FORMAT="Thread Affinity: %0.3L %.8n %.15{thread_affinity} %.12H"
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OMP_MAX_TASK_PRIORITY=20

for t in 1 2 4 8 16 32 64 128 192; do
  echo "== pf threads=$t"
  OMP_NUM_THREADS=$t ./jacobi_pf --n=${NB} --N=${NX} --nt=100
done

unset OMP_NUM_THREADS
for r in 1 2 4 8 16 32 64 128 256 512 1024; do
  nb=$((NX / r))
  echo "== ompss NX/nb=$r nb=$nb"
  taskset -c 0-191 ./jacobi_ompss --n=$nb --NX=$NX --NY=$NY --nt=100
done
