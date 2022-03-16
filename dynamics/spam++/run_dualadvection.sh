
rm *.png *.nc dualadvection
./build/cmake_linuxlaptop.sh dualadvection $2
make
mpirun.mpich -n $1 ./dualadvection input-advection.txt
python3 plot_dualadvection.py $2
