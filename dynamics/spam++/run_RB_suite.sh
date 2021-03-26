# for some list of nx, ny, initcond, ndims, compile parameters, etc.





recon_types = ['WENOFUNC','CFV']
recon_orders = [1,3,5,9]
diff_orders = [2,4]
sizes = [100,]

mkdir -p RB

DSdict = {}
energy_dict = {}
pv_dict = {}
mass_dict = {}
for recon_type in recon_types:
    for recon_order in recon_orders:
        for diff_order in diff_orders:
            for size in sizes:
*********************
IDEALLY HERE THERE IS A COMPILE TIME OPTIONS FILE THAT THE BUILD SCRIPT READS FROM AND USES TO SET STUFF

THEN THIS SCRIPT CAN JUST POINT TO APPROPRIATE VERSION OF THAT FILE, OR AUTOMATICALLY BUILD IT FROM A TEMPLATE
SET RECONTYPE, RECON ORDER, DIFF ORDER
MOVE RELEVANT INPUT FILE IN HERE- OR MAKE THIS PART OF RUN LAYER MODEL IE RUN LAYER MODEL TAKES INPUT FILE AS A ARGUMENT?
mkdir -p RB/RB$SIZE-$RECONTYPE$RECONORDER-HODGE$DIFFORDER
cd RB/RB$SIZE-$RECONTYPE$RECONORDER-HODGE$DIFFORDER
rm *.png *.nc
cd ../..
./run_layermodel2D.sh 4 ce 0 0
mv *.png *.nc RB/RB$SIZE-$RECONTYPE$RECONORDER-HODGE$DIFFORDER

***************

python3 plot_RB.py
