rm *.png *.nc driver
./make_script.sh
mpirun.mpich -n $3 ./driver ../inputs/pamc_input_layer_$2.yaml
python3 plot_layermodel2D.py $1
