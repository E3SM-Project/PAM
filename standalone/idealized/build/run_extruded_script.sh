rm *.png
./make_script.sh
mpirun.mpich -n 1 ./driver ../inputs/pamc_input_extruded_$2.yaml
python3 plot_extrudedmodel2D.py $1
