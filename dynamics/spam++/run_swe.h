
rm *.png *.nc swe
make
mpirun -n $1 ./swe input.txt
python3 plot_swe.py
