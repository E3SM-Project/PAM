


# HOW DO WE ADD RHO VS. RHO_D VARIANTS HERE?

#rm *.png *.nc *.build *.out layermodel2D
mkdir -p MLRB
for SIZE in 200 400
do
  for DIFFORDER in 2 4
  do
    for RECONTYPE in 'WENOFUNC' 'CFV'
    do
      for RECONORDER in 1 3 5 7 9
      do
        mkdir -p MLRB/MLRB${SIZE}-${RECONTYPE}${RECONORDER}-HODGE${DIFFORDER}
        cd MLRB/MLRB${SIZE}-${RECONTYPE}${RECONORDER}-HODGE${DIFFORDER}
        rm *.png *.nc *.build *.out
        cd ../..
        cp runsuites/MLRB/compile_consts.template build_consts.build
        sed -i "s/HODGEORDER/${DIFFORDER}/" build_consts.build
        sed -i "s/RECONORDER/${RECONORDER}/" build_consts.build
        sed -i "s/RECONTYPE/${RECONTYPE}/" build_consts.build
        cp runsuites/MLRB/model_consts.template model_consts.build
#FIX THIS TO MODIFY RHO VS. RHO_D
        sed -i "s/HODGEORDER/${DIFFORDER}/" model_consts.build
        #./build/cmake_linuxlaptop.sh mce2D build_consts.build model_consts.build >& build.out
        #make
        #mpirun.mpich -n 4 ./layermodel2D runsuites/MLRB/MLRB${SIZE}.input >& run.out
        mv *.png *.nc *.build *.out MLRB/MLRB${SIZE}-${RECONTYPE}${RECONORDER}-HODGE${DIFFORDER}
      done
    done  
  done
done
 
#python3 plot_MLRB.py
#mv *.png MLRB

