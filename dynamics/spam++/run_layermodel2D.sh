
rm *.png *.nc layermodel2D
./build/cmake_linuxlaptop.sh $22D $3
make
mpirun.mpich -n $1 ./layermodel2D input-layermodel2D.txt
python3 plot_layermodel2D.py $2 $3

# ./run_layeredmodel2D.sh NPROCS EQNSET NTRACERS
