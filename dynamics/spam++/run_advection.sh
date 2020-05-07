
rm *.png *.nc advection
./build/cmake_linuxlaptop.sh advection
make
mpirun.mpich -n $1 ./advection input-advection.txt
python3 plot_advection.py
