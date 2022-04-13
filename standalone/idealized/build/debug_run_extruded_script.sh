rm *.png
./make_debug_script.sh
mpirun.mpich -n 1 ./driver ../inputs/pamc_input_extruded_$2.yaml
