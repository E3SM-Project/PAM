# Workflow for building standalone

## Clone PAM

```bash
# cd into your home directory or a place you want to clone PAM
git clone git@github.com:E3SM-Project/PAM.git
cd PAM
git submodule update --init --recursive
```

## If you want to build with C++ scream physics

```bash
# cd into your home directory or a place you want to clone scream
git clone git@github.com:E3SM-Project/scream.git
cd scream
git submodule update --init --recursive
cd /path/to/PAM
cd standalone/machines
# edit linux_laptop_gnu_mpi_cpu.env to match your paths and options
#   Set PAM_SCREAM_USE_CXX to ON for C++ physics and to OFF to use Fortran CPU physics
source linux_laptop_gnu_mpi_cpu.env
#####################################
## Build P3 nd SHOC
#####################################
cd ../mmf_simplified/build_p3_shoc_cxx
./cmakescript.sh
make -j
#########################################
## Build P3 nd SHOC interfaces for PAM
#########################################
cd ../build_p3_shoc_cxx_interface
./cmakescript.sh
make -j
```

## Build PAM

```bash
cd /path/to/PAM
cd standalone/mmf_simplified/build
# Edit linus_laptop_gnu_mpi_cpu.env to match your paths and options
source ../../machines/linux_laptop_gnu_mpi_cpu.env
# Edit cmakescript to match your options
./cmakescript.sh
make -j
```

## Run PAM

```bash
cd /path/to/PAM
cd standalone/mmf_simplified/build
./driver ../inputs/input_euler3d.yaml
```


