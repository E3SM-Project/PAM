# SCAMPAM: SCAlable, Modular, and Portable Atmospheric Model

PAM is a portable atmospheric model written in C++ with performance portability in mind using a kernel launching approach. It works on CPUs, Nvidia GPUs, and AMD GPUs. The goal for PAM is to allow easily exchangeable dynamical core and physics options through a simple and clear interface layer that's common to each. An emphasis is placed on algorithms that give good hardware uitilization on accelerated architectures and MPI patterns that give good scalability.

## Directory Structure

* `coupler`: DataManager class and coupled state allocator function
* `dynamics`: Dynamical cores, currently "AWFL"
* `externals`: git submodule dependencies
* `include`: PAM include files for the whole project
* `physics`: physics modules, currently micro only
  * `micro`: The various microphysics options for PAM
* `standalone`: standalone driver for PAM
* `utils`: Various PAM utilities, mostly in python

## Running standalone PAM

1. Clone the repo
2. `cd SCAMPAM && git submodule update --init --recursive`
3. `cd standalone/build`
4. Edit `cmakescript.sh` to choose your dycore and microphysics option
5. `source [machine].sh`
6. `./cmakescript`
7. Use `SCAMPAM/utils/generate_vertical_levels.py` to generate a vertical coordinates file (see the file's documentation for help)
8. Edit `../inputs/input_euler3d.yaml` to whatever parameters you want
9. `./driver ../inputs/input_euler3d.yaml`

