rm *.png *.nc
cd ../build
rm driver
./../pam-c/make_script.sh $3
mpirun.mpich -n $3 ./driver ../inputs/pamc_input_layer_$2.yaml
cd ../pam-c
mv ../build/*.nc .
python3 plot_layermodel2D.py $1
