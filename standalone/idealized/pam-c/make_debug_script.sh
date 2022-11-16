cd ../build
source ../../machines/linux_laptop_gnu_mpi_cpu_debug.env
./cmakescript.sh ../pam-c/set_pamc_cmakevars.sh
make -j $1
