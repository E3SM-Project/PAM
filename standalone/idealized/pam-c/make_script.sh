cd ../build
source ../../machines/linux_laptop_gnu_mpi_cpu.env
./cmakescript_pamc.sh ../pam-c/set_pamc_cmakevars.sh
make -j $1
