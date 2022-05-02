

rm *.png *.nc *.build *.out layermodel2D extrudedmodel2D
mkdir -p RB
for SIZE in 100
do
  for DIFFORDER in 2
  #for DIFFORDER in 2 4
  do
    #for RECONTYPE in 'WENOFUNC' 'CFV'
    for RECONTYPE in 'WENOFUNC'
    do
      for RECONORDER in 1 3 5 7 9
      do
        mkdir -p RB/RB${SIZE}-${RECONTYPE}${RECONORDER}-HODGE${DIFFORDER}
        cd RB/RB${SIZE}-${RECONTYPE}${RECONORDER}-HODGE${DIFFORDER}
        rm *.png *.nc *.build *.out
        cd ../..
        cp runsuites/RB/compile_consts.template build_consts.build
        sed -i "s/VERTHODGEORDER/${DIFFORDER}/" build_consts.build
        sed -i "s/VERTRECONORDER/${RECONORDER}/" build_consts.build
        sed -i "s/VERTRECONTYPE/${RECONTYPE}/" build_consts.build
        sed -i "s/HODGEORDER/${DIFFORDER}/" build_consts.build
        sed -i "s/RECONORDER/${RECONORDER}/" build_consts.build
        sed -i "s/RECONTYPE/${RECONTYPE}/" build_consts.build
        #./build/cmake_linuxlaptop.sh ce2D build_consts.build &> build.out
        ./build/cmake_linuxlaptop.sh ce2Dext build_consts.build &> build.out
        make &>> build.out
        #mpirun.mpich -n 4 ./layermodel2D runsuites/RB/RB${SIZE}.input &> run.out
        mpirun.mpich -n 4 ./extrudedmodel2D runsuites/RB/RB${SIZE}.input &> run.out
        mv *.png *.nc *.build *.out RB/RB${SIZE}-${RECONTYPE}${RECONORDER}-HODGE${DIFFORDER}
	cp runsuites/RB/RB${SIZE}.input RB/RB${SIZE}-${RECONTYPE}${RECONORDER}-HODGE${DIFFORDER}
      done
    done  
  done
done
 
python3 plot_RB.py
mv *.png RB
