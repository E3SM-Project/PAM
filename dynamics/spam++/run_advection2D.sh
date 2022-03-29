
rm *.png *.nc advection
./build/cmake_linuxlaptop.sh advection $2 $3 $4
make
mpirun.mpich -n $1 ./advection input-advection.txt
python3 plot_advection.py $2 $3 $4

# ./run_advection.sh NPROCS NTRACERS NTRACERS_FCT NTRACERS_Q
