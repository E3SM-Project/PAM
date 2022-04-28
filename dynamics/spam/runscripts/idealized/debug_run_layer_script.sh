rm *.png *.nc driver
./make_debug_script.sh
mpirun.mpich -n $3 ./driver ../inputs/pamc_input_layer_$2.yaml
