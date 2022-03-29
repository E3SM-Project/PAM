
rm *.png *.nc extrudedmodel2D
./build/cmake_linuxlaptop.sh $22Dext
make
mpirun.mpich -n $1 ./extrudedmodel2D input-extrudedmodel2D.txt
python3 plot_extrudedmodel2D.py $2

# ./run_extrudedmodel2D.sh NPROCS EQNSET(swe/tswe/ce/mce) 
# ./run_extrudedmodel2D.sh 4 tswe
# ./run_extrudedmodel2D.sh 4 ce
