
rm *.png *.nc layermodel2D
./build/cmake_linuxlaptop.sh $22D $3 $4 0
make
mpirun.mpich -n $1 ./layermodel2D input-layermodel2D.txt
python3 plot_layermodel2D.py $2 $3 $4

# ./run_layermodel2D.sh NPROCS EQNSET(swe/tswe) NTRACERS NTRACERS_FCT
# ./run_layermodel2D.sh 4 tswe 2 2

# ./run_layermodel2D.sh NPROCS EQNSET(ce) NTRACERS (IGNORED) NTRACERS_FCT (IGNORED) THERMO (idealgasentropy, idealgaspottemp)
# ./run_layermodel2D.sh NPROCS EQNSET(mcerho/mcerhod) NTRACERS (IGNORED) NTRACERS_FCT (IGNORED) THERMO (constkappavirtpottemp)
