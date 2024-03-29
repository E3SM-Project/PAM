#clean up any existing files
rm *.png *.nc
cd ../build
rm driver

#build model
source ../../machines/linux_laptop_gnu_mpi_cpu.env
./cmakescript_pamc.sh PAM_SGS=none PAM_MICRO=none PAMC_MODEL=extrudedmodel PAMC_HAMIL=an PAMC_THERMO=idealgaspottemp PAMC_IO=serial
make -j 4

#linux_laptop_gnu_mpi_cpu_debug linux_laptop_gnu_mpi_cpu
#p3 none
#shoc none
#ce mce_rho an man
#idealgaspottemp constkappavirpottemp


#run model
mpirun.mpich -n $1 ./driver ../inputs/pamc_idealized/pamc_input_extruded_densitycurrent.yaml
cd ../pam-c
mv ../build/*.nc .

#pamc_input_extruded_densitycurrent
#pamc_input_extruded_gravitywave
#pamc_input_extruded_largerisingbubble
#pamc_input_extruded_moistrisingbubble
#pamc_input_extruded_risingbubble
#pamc_input_extruded_twobubbles


#plot model
python3 plot_extrudedmodel2D.py an

#an ce man mce
