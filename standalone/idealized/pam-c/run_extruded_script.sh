rm *.png *.nc
cd ../build
rm driver
./../pam-c/make_script.sh $3
mpirun.mpich -n $3 ./driver ../inputs/pamc_input_extruded_$2.yaml
cd ../pam-c
mv ../build/*.nc .
python3 plot_extrudedmodel2D.py $1
#python3 plot_extrudedmodel2D_parallel.py $1 ../inputs/pamc_input_extruded_$2.yaml $3
