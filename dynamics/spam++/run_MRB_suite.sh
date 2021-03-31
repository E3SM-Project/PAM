

rm *.png *.nc *.build *.out layermodel2D
mkdir -p MRB
for SIZE in 100
do
  for DIFFORDER in 2
  do
    for RECONTYPE in 'WENOFUNC' 'CFV'
    do
      for RECONORDER in 1 3 5 7 9
      do
        mkdir -p MRB/MRB${SIZE}-${RECONTYPE}${RECONORDER}-HODGE${DIFFORDER}
        cd MRB/MRB${SIZE}-${RECONTYPE}${RECONORDER}-HODGE${DIFFORDER}
        rm *.png *.nc *.build *.out
        cd ../..
        cp runsuites/MRB/compile_consts.template build_consts.build
        sed -i "s/HODGEORDER/${DIFFORDER}/" build_consts.build
        sed -i "s/RECONORDER/${RECONORDER}/" build_consts.build
        sed -i "s/RECONTYPE/${RECONTYPE}/" build_consts.build
        ./build/cmake_linuxlaptop.sh mce2D build_consts.build &> build.out
        make &>> build.out
        mpirun.mpich -n 4 ./layermodel2D runsuites/MRB/MRB${SIZE}.input &> run.out
        mv *.png *.nc *.build *.out MRB/MRB${SIZE}-${RECONTYPE}${RECONORDER}-HODGE${DIFFORDER}
      done
    done  
  done
done
 
python3 plot_MRB.py
mv *.png MRB
