
PNETCDF_PATH=/ccs/home/imn/parallel-netcdf-1.11.2_clang

CXX = mpic++ -std=c++14
LINK = mpic++
CXXFLAGS = -O3 -ffast-math -I${PNETCDF_PATH}/include -I/ccs/home/imn/hipcub/hipcub/include -I/ccs/home/imn/rocprim/include -D__USE_HIP__ # --amdgpu-target=gfx906 -D__HIP_ARCH_GFX906__=1
LDFLAGS = -L${PNETCDF_PATH}/lib -lpnetcdf

export OMPI_CXX=hipcc

