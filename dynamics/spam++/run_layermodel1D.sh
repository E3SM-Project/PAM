
rm *.png *.nc layermodel1D
./build/cmake_linuxlaptop.sh $21D
make
mpirun.mpich -n $1 ./layermodel1D input-layermodel1D.txt
python3 plot_layermodel1D.py $2

# ./run_layermodel1D.sh NPROCS EQNSET(swe/tswe) 
# ./run_layermodel1D.sh 4 tswe
