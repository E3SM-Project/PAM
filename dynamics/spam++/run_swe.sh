
rm *.png *.nc swe
./build/cmake_linuxlaptop.sh swe
make
mpirun.mpich -n $1 ./swe input-swe.txt
python3 plot_swe.py
