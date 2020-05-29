
rm *.png *.nc layermodel
./build/cmake_linuxlaptop.sh $2 $3
make
mpirun.mpich -n $1 ./layermodel input-layermodel.txt
python3 plot_layermodel.py $2 $3
