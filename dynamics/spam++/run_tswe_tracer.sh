
rm *.png *.nc tswe_tracer
./build/cmake_linuxlaptop.sh tswe_tracer
make
mpirun.mpich -n $1 ./tswe_tracer input-tswe.txt
python3 plot_tswe_tracer.py
