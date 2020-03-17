
rm *.png *.nc dualadvection
make
mpirun -n $1 ./dualadvection input.txt
python3 plot_dualadvection.py
