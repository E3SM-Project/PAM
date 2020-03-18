
rm *.png *.nc swe
./build/cmake_linuxlaptop.sh
make
mpirun -n $1 ./swe input.txt
python3 plot_swe.py
