cd ../build
source ../../machines/linux_laptop_debug.env
./cmakescript.sh ../pam-c/set_pamc_cmakevars.sh
make -j
