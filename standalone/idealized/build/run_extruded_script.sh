rm *.png *.nc driver
./make_script.sh
mpirun.mpich -n $3 ./driver ../inputs/pamc_input_extruded_$2.yaml
#python3 plot_extrudedmodel2D_parallel.py $1 ../inputs/pamc_input_extruded_$2.yaml $3
python3 plot_extrudedmodel2D.py $1
