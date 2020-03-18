
rm *.png *.nc dualadvection
./build/cmake_linuxlaptop.sh
make
mpirun -n $1 ./dualadvection input.txt
python3 plot_dualadvection.py
