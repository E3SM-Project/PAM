
rm *.png *.nc advection
make
mpirun -n $1 ./advection input.txt
python3 plot_advection.py
