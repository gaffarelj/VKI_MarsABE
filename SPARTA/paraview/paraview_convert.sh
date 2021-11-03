#!/bin/sh
cd surf
rm -rf *
pvpython ../../tools/surf2paraview.py ../../setup/data/data.CS_0021 CS_0021 
cd ../grid

rm -rf CS_0021_85km_0 
rm -rf CS_0021_85km_0.pvd 
echo 'Converting result grid of CS_0021 at 85km (refinement 0) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_0021_85km_0 vals_CS_0021_85km_0 -r ../../setup/results_sparta/CS_0021/gridvals_85km_0.*.dat 
pvpython ../../tools/grid2paraview_original.py def/grid.CS_0021_85km_0 Kn_CS_0021_85km_0 -r ../../setup/results_sparta/CS_0021/gridKn_85km_0.*.dat 

rm -rf CS_0021_85km_1 
rm -rf CS_0021_85km_1.pvd 
echo 'Converting result grid of CS_0021 at 85km (refinement 1) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_0021_85km_1 vals_CS_0021_85km_1 -r ../../setup/results_sparta/CS_0021/gridvals_85km_1.*.dat 
pvpython ../../tools/grid2paraview_original.py def/grid.CS_0021_85km_1 Kn_CS_0021_85km_1 -r ../../setup/results_sparta/CS_0021/gridKn_85km_1.*.dat 

rm -rf CS_0021_85km_2 
rm -rf CS_0021_85km_2.pvd 
echo 'Converting result grid of CS_0021 at 85km (refinement 2) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_0021_85km_2 vals_CS_0021_85km_2 -r ../../setup/results_sparta/CS_0021/gridvals_85km_2.*.dat 
pvpython ../../tools/grid2paraview_original.py def/grid.CS_0021_85km_2 Kn_CS_0021_85km_2 -r ../../setup/results_sparta/CS_0021/gridKn_85km_2.*.dat 
