# Portable Atmosphere Model (PAM)

PAM is a portable atmospheric model written in C++ with performance portability in mind using a kernel launching approach. It works on CPUs, Nvidia GPUs, and AMD GPUs. The goal for PAM is to allow easily exchangeable dynamical core and physics options through a simple and clear interface layer that's common to each. An emphasis is placed on algorithms that give good hardware uitilization on accelerated architectures and MPI patterns that give good scalability.

git clone https://github.com/E3SM-Project/PAM.git
git submodule update --init
