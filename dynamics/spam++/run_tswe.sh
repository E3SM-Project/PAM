
rm *.png *.nc tswe
./build/cmake_linuxlaptop.sh tswe
make
mpirun.mpich -n $1 ./tswe input-tswe.txt
python3 plot_tswe.py
