#clean up any existing files
rm *.png *.nc
cd ../build
rm driver

#build model
source ../../machines/linux_laptop_gnu_mpi_cpu.env
./cmakescript_pamc.sh PAM_SGS=none PAM_MICRO=none PAMC_MODEL=layermodel PAMC_HAMIL=tswe PAMC_THERMO=none PAMC_IO=serial
make -j 4

#linux_laptop_gnu_mpi_cpu_debug linux_laptop_gnu_mpi_cpu
#swe tswe

#run model
mpirun.mpich -n $1 ./driver ../inputs/pamc_idealized/pamc_input_layer_doublevortex.yaml
cd ../pam-c
mv ../build/*.nc .

#pamc_input_layer_bickleyjet
#pamc_input_layer_doublevortex


#plot model
python3 plot_layermodel2D.py tswe

#swe tswe
