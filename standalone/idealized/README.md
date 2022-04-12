# Instructions for running idealized standalone

* `cd PAM/standalone/idealized/build`
* Generate the vertical coordinates you want, e.g.,
```bash
python3 ../../../utils/generate_vertical_levels.py --function=equal \
                                                   --nlev=100       \
                                                   --ztop=20000     \
                                                   --output=vcoords_equal_100_20km.nc
```
* Setup the configuration you want (dycore, micro, sgs) in `cmakescript.sh`
* Configure the runtime options in `../input/input_[option].yaml`
* Source the correct machine file: `source ../../machines/[machine]_[options].sh`
* `./cmakescript.sh`
* `make -j`
* `./driver ../input/input_[option].sh`
