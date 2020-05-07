
rm *.png *.nc dualadvection
./build/cmake_linuxlaptop.sh dualadvection
make
mpirun.mpich -n $1 ./dualadvection input-advection.txt
python3 plot_dualadvection.py
