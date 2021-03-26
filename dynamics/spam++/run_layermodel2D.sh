
rm *.png *.nc layermodel2D
./build/cmake_linuxlaptop.sh $22D
make
mpirun.mpich -n $1 ./layermodel2D input-layermodel2D.txt
python3 plot_layermodel2D.py $2

# ./run_layermodel2D.sh NPROCS EQNSET(swe/tswe) 
# ./run_layermodel2D.sh 4 tswe
# ./run_layermodel2D.sh 4 ce
