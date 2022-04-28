cd ../build
source ../../machines/linux_laptop.env
./cmakescript.sh ../pam-c/set_pamc_cmakevars.sh
make -j
