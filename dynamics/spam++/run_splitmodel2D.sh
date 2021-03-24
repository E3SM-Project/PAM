
rm *.png *.nc splitmodel2D
# settings NTRACERS and NTRACERS_FCT to 0 for now...
./build/cmake_linuxlaptop.sh split$22D 0 0
make
mpirun.mpich -n $1 ./splitmodel2D input-splitmodel2D.txt
python3 plot_splitmodel2D.py $2 $3 $4

# ./run_splitmodel2D.sh NPROCS EQNSET
