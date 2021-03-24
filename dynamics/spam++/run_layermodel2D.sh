
rm *.png *.nc layermodel2D
./build/cmake_linuxlaptop.sh $22D $3 $4 0
make
mpirun.mpich -n $1 ./layermodel2D input-layermodel2D.txt
python3 plot_layermodel2D.py $2 $3 $4

# ./run_layermodel2D.sh NPROCS EQNSET NTRACERS NTRACERS_FCT
# ./run_layermodel2D.sh 4 tswe 2 2
