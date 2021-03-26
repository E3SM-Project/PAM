# for some list of nx, ny, initcond, ndims, compile parameters, etc.





recon_types = ['WENOFUNC','CFV']
recon_orders = [1,3,5,9]
diff_orders = [2,4]
sizes = [100,200]

mkdir -p RB

for size in sizes:
for recon_type in recon_types:
    for recon_order in recon_orders:
        for diff_order in diff_orders:

mkdir -p RB/RB$SIZE-$RECONTYPE$RECONORDER-HODGE$DIFFORDER
cd RB/RB$SIZE-$RECONTYPE$RECONORDER-HODGE$DIFFORDER
rm *.png *.nc
cd ../..
./run_layermodel2D.sh 4 ce
cp runsuites/RB/compile_consts.template build_consts.build
MODIFY IT!

mv *.png *.nc *.build RB/RB$SIZE-$RECONTYPE$RECONORDER-HODGE$DIFFORDER

***************

python3 plot_RB.py
