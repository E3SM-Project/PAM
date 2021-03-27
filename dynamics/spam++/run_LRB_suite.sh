

rm *.png *.nc *.build *.out layermodel2D
mkdir -p LRB
for SIZE in 200 400
do
  for DIFFORDER in 2 4
  do
    for RECONTYPE in 'WENOFUNC' 'CFV'
    do
      for RECONORDER in 1 3 5 7 9
      do
        mkdir -p LRB/LRB${SIZE}-${RECONTYPE}${RECONORDER}-HODGE${DIFFORDER}
        cd LRB/LRB${SIZE}-${RECONTYPE}${RECONORDER}-HODGE${DIFFORDER}
        rm *.png *.nc *.build *.out
        cd ../..
        cp runsuites/LRB/compile_consts.template build_consts.build
        sed -i "s/HODGEORDER/${DIFFORDER}/" build_consts.build
        sed -i "s/RECONORDER/${RECONORDER}/" build_consts.build
        sed -i "s/RECONTYPE/${RECONTYPE}/" build_consts.build
        ./build/cmake_linuxlaptop.sh ce2D build_consts.build &> build.out
        make &>> build.out
        mpirun.mpich -n 4 ./layermodel2D runsuites/LRB/LRB${SIZE}.input &> run.out
        mv *.png *.nc *.build *.out LRB/LRB${SIZE}-${RECONTYPE}${RECONORDER}-HODGE${DIFFORDER}
      done
    done  
  done
done
 
python3 plot_LRB.py
mv *.png LRB

