rm *.png
./make_script.sh
mpirun.mpich -n 1 ./driver ../inputs/pamc_input_layer.yaml
python3 plot_layermodel2D.py swe
