
rm *.png *.nc layermodel1D
./build/cmake_linuxlaptop.sh $21D $3
make
mpirun.mpich -n $1 ./layermodel1D input-layermodel1D.txt
python3 plot_layermodel1D.py $2 $3
