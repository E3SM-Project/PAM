
rm *.png *.nc advection
./build/cmake_linuxlaptop.sh
make
mpirun -n $1 ./advection input.txt
python3 plot_advection.py
