
rm *.png *.nc extrudedadvection2D
./build/cmake_linuxlaptop.sh extrudedadvection2D
make
mpirun.mpich -n $1 ./extrudedadvection2D input-extrudedadvection2D.txt
python3 plot_extrudedadvection2D.py
