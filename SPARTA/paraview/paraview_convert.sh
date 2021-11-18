#!/bin/sh
cd surf
rm -rf *
pvpython ../../tools/surf2paraview.py ../../setup/data/data.CS_1021 CS_1021 
pvpython ../../tools/surf2paraview.py ../../setup/data/data.CS_2021 CS_2021 
pvpython ../../tools/surf2paraview.py ../../setup/data/data.CS_2120 CS_2120 
pvpython ../../tools/surf2paraview.py ../../setup/data/data.CS_3021 CS_3021 
cd ../grid

rm -rf vals_CS_1021_85km_0 
rm -rf vals_CS_1021_85km_0.pvd 
echo 'Converting results of CS_1021 at 85km (refinement 0) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_1021_85km_0 vals_CS_1021_85km_0 -r ../../setup/results_sparta/CS_1021/vals_85km_0.*.dat 

rm -rf vals_CS_1021_85km_1 
rm -rf vals_CS_1021_85km_1.pvd 
echo 'Converting results of CS_1021 at 85km (refinement 1) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_1021_85km_1 vals_CS_1021_85km_1 -r ../../setup/results_sparta/CS_1021/vals_85km_1.*.dat 

rm -rf vals_CS_1021_85km_2 
rm -rf vals_CS_1021_85km_2.pvd 
echo 'Converting results of CS_1021 at 85km (refinement 2) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_1021_85km_2 vals_CS_1021_85km_2 -r ../../setup/results_sparta/CS_1021/vals_85km_2.*.dat 

rm -rf vals_CS_1021_85km_3 
rm -rf vals_CS_1021_85km_3.pvd 
echo 'Converting results of CS_1021 at 85km (refinement 3) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_1021_85km_3 vals_CS_1021_85km_3 -r ../../setup/results_sparta/CS_1021/vals_85km_3.*.dat 

rm -rf vals_CS_1021_115km_0 
rm -rf vals_CS_1021_115km_0.pvd 
echo 'Converting results of CS_1021 at 115km (refinement 0) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_1021_115km_0 vals_CS_1021_115km_0 -r ../../setup/results_sparta/CS_1021/vals_115km_0.*.dat 

rm -rf vals_CS_1021_115km_1 
rm -rf vals_CS_1021_115km_1.pvd 
echo 'Converting results of CS_1021 at 115km (refinement 1) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_1021_115km_1 vals_CS_1021_115km_1 -r ../../setup/results_sparta/CS_1021/vals_115km_1.*.dat 

rm -rf vals_CS_1021_115km_2 
rm -rf vals_CS_1021_115km_2.pvd 
echo 'Converting results of CS_1021 at 115km (refinement 2) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_1021_115km_2 vals_CS_1021_115km_2 -r ../../setup/results_sparta/CS_1021/vals_115km_2.*.dat 

rm -rf vals_CS_1021_115km_3 
rm -rf vals_CS_1021_115km_3.pvd 
echo 'Converting results of CS_1021 at 115km (refinement 3) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_1021_115km_3 vals_CS_1021_115km_3 -r ../../setup/results_sparta/CS_1021/vals_115km_3.*.dat 

rm -rf vals_CS_1021_150km_0 
rm -rf vals_CS_1021_150km_0.pvd 
echo 'Converting results of CS_1021 at 150km (refinement 0) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_1021_150km_0 vals_CS_1021_150km_0 -r ../../setup/results_sparta/CS_1021/vals_150km_0.*.dat 

rm -rf vals_CS_1021_150km_1 
rm -rf vals_CS_1021_150km_1.pvd 
echo 'Converting results of CS_1021 at 150km (refinement 1) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_1021_150km_1 vals_CS_1021_150km_1 -r ../../setup/results_sparta/CS_1021/vals_150km_1.*.dat 

rm -rf vals_CS_1021_150km_2 
rm -rf vals_CS_1021_150km_2.pvd 
echo 'Converting results of CS_1021 at 150km (refinement 2) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_1021_150km_2 vals_CS_1021_150km_2 -r ../../setup/results_sparta/CS_1021/vals_150km_2.*.dat 

rm -rf vals_CS_1021_150km_3 
rm -rf vals_CS_1021_150km_3.pvd 
echo 'Converting results of CS_1021 at 150km (refinement 3) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_1021_150km_3 vals_CS_1021_150km_3 -r ../../setup/results_sparta/CS_1021/vals_150km_3.*.dat 

rm -rf vals_CS_2021_85km_0 
rm -rf vals_CS_2021_85km_0.pvd 
echo 'Converting results of CS_2021 at 85km (refinement 0) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2021_85km_0 vals_CS_2021_85km_0 -r ../../setup/results_sparta/CS_2021/vals_85km_0.*.dat 

rm -rf vals_CS_2021_85km_1 
rm -rf vals_CS_2021_85km_1.pvd 
echo 'Converting results of CS_2021 at 85km (refinement 1) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2021_85km_1 vals_CS_2021_85km_1 -r ../../setup/results_sparta/CS_2021/vals_85km_1.*.dat 

rm -rf vals_CS_2021_85km_2 
rm -rf vals_CS_2021_85km_2.pvd 
echo 'Converting results of CS_2021 at 85km (refinement 2) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2021_85km_2 vals_CS_2021_85km_2 -r ../../setup/results_sparta/CS_2021/vals_85km_2.*.dat 

rm -rf vals_CS_2021_85km_3 
rm -rf vals_CS_2021_85km_3.pvd 
echo 'Converting results of CS_2021 at 85km (refinement 3) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2021_85km_3 vals_CS_2021_85km_3 -r ../../setup/results_sparta/CS_2021/vals_85km_3.*.dat 

rm -rf vals_CS_2021_115km_0 
rm -rf vals_CS_2021_115km_0.pvd 
echo 'Converting results of CS_2021 at 115km (refinement 0) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2021_115km_0 vals_CS_2021_115km_0 -r ../../setup/results_sparta/CS_2021/vals_115km_0.*.dat 

rm -rf vals_CS_2021_115km_1 
rm -rf vals_CS_2021_115km_1.pvd 
echo 'Converting results of CS_2021 at 115km (refinement 1) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2021_115km_1 vals_CS_2021_115km_1 -r ../../setup/results_sparta/CS_2021/vals_115km_1.*.dat 

rm -rf vals_CS_2021_115km_2 
rm -rf vals_CS_2021_115km_2.pvd 
echo 'Converting results of CS_2021 at 115km (refinement 2) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2021_115km_2 vals_CS_2021_115km_2 -r ../../setup/results_sparta/CS_2021/vals_115km_2.*.dat 

rm -rf vals_CS_2021_115km_3 
rm -rf vals_CS_2021_115km_3.pvd 
echo 'Converting results of CS_2021 at 115km (refinement 3) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2021_115km_3 vals_CS_2021_115km_3 -r ../../setup/results_sparta/CS_2021/vals_115km_3.*.dat 

rm -rf vals_CS_2021_150km_0 
rm -rf vals_CS_2021_150km_0.pvd 
echo 'Converting results of CS_2021 at 150km (refinement 0) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2021_150km_0 vals_CS_2021_150km_0 -r ../../setup/results_sparta/CS_2021/vals_150km_0.*.dat 

rm -rf vals_CS_2021_150km_1 
rm -rf vals_CS_2021_150km_1.pvd 
echo 'Converting results of CS_2021 at 150km (refinement 1) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2021_150km_1 vals_CS_2021_150km_1 -r ../../setup/results_sparta/CS_2021/vals_150km_1.*.dat 

rm -rf vals_CS_2021_150km_2 
rm -rf vals_CS_2021_150km_2.pvd 
echo 'Converting results of CS_2021 at 150km (refinement 2) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2021_150km_2 vals_CS_2021_150km_2 -r ../../setup/results_sparta/CS_2021/vals_150km_2.*.dat 

rm -rf vals_CS_2021_150km_3 
rm -rf vals_CS_2021_150km_3.pvd 
echo 'Converting results of CS_2021 at 150km (refinement 3) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2021_150km_3 vals_CS_2021_150km_3 -r ../../setup/results_sparta/CS_2021/vals_150km_3.*.dat 

rm -rf vals_CS_2120_85km_0 
rm -rf vals_CS_2120_85km_0.pvd 
echo 'Converting results of CS_2120 at 85km (refinement 0) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2120_85km_0 vals_CS_2120_85km_0 -r ../../setup/results_sparta/CS_2120/vals_85km_0.*.dat 

rm -rf vals_CS_2120_85km_1 
rm -rf vals_CS_2120_85km_1.pvd 
echo 'Converting results of CS_2120 at 85km (refinement 1) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2120_85km_1 vals_CS_2120_85km_1 -r ../../setup/results_sparta/CS_2120/vals_85km_1.*.dat 

rm -rf vals_CS_2120_85km_2 
rm -rf vals_CS_2120_85km_2.pvd 
echo 'Converting results of CS_2120 at 85km (refinement 2) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2120_85km_2 vals_CS_2120_85km_2 -r ../../setup/results_sparta/CS_2120/vals_85km_2.*.dat 

rm -rf vals_CS_2120_85km_3 
rm -rf vals_CS_2120_85km_3.pvd 
echo 'Converting results of CS_2120 at 85km (refinement 3) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2120_85km_3 vals_CS_2120_85km_3 -r ../../setup/results_sparta/CS_2120/vals_85km_3.*.dat 

rm -rf vals_CS_2120_115km_0 
rm -rf vals_CS_2120_115km_0.pvd 
echo 'Converting results of CS_2120 at 115km (refinement 0) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2120_115km_0 vals_CS_2120_115km_0 -r ../../setup/results_sparta/CS_2120/vals_115km_0.*.dat 

rm -rf vals_CS_2120_115km_1 
rm -rf vals_CS_2120_115km_1.pvd 
echo 'Converting results of CS_2120 at 115km (refinement 1) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2120_115km_1 vals_CS_2120_115km_1 -r ../../setup/results_sparta/CS_2120/vals_115km_1.*.dat 

rm -rf vals_CS_2120_115km_2 
rm -rf vals_CS_2120_115km_2.pvd 
echo 'Converting results of CS_2120 at 115km (refinement 2) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2120_115km_2 vals_CS_2120_115km_2 -r ../../setup/results_sparta/CS_2120/vals_115km_2.*.dat 

rm -rf vals_CS_2120_115km_3 
rm -rf vals_CS_2120_115km_3.pvd 
echo 'Converting results of CS_2120 at 115km (refinement 3) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2120_115km_3 vals_CS_2120_115km_3 -r ../../setup/results_sparta/CS_2120/vals_115km_3.*.dat 

rm -rf vals_CS_2120_150km_0 
rm -rf vals_CS_2120_150km_0.pvd 
echo 'Converting results of CS_2120 at 150km (refinement 0) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2120_150km_0 vals_CS_2120_150km_0 -r ../../setup/results_sparta/CS_2120/vals_150km_0.*.dat 

rm -rf vals_CS_2120_150km_1 
rm -rf vals_CS_2120_150km_1.pvd 
echo 'Converting results of CS_2120 at 150km (refinement 1) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2120_150km_1 vals_CS_2120_150km_1 -r ../../setup/results_sparta/CS_2120/vals_150km_1.*.dat 

rm -rf vals_CS_2120_150km_2 
rm -rf vals_CS_2120_150km_2.pvd 
echo 'Converting results of CS_2120 at 150km (refinement 2) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2120_150km_2 vals_CS_2120_150km_2 -r ../../setup/results_sparta/CS_2120/vals_150km_2.*.dat 

rm -rf vals_CS_2120_150km_3 
rm -rf vals_CS_2120_150km_3.pvd 
echo 'Converting results of CS_2120 at 150km (refinement 3) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_2120_150km_3 vals_CS_2120_150km_3 -r ../../setup/results_sparta/CS_2120/vals_150km_3.*.dat 

rm -rf vals_CS_3021_85km_0 
rm -rf vals_CS_3021_85km_0.pvd 
echo 'Converting results of CS_3021 at 85km (refinement 0) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_3021_85km_0 vals_CS_3021_85km_0 -r ../../setup/results_sparta/CS_3021/vals_85km_0.*.dat 

rm -rf vals_CS_3021_85km_1 
rm -rf vals_CS_3021_85km_1.pvd 
echo 'Converting results of CS_3021 at 85km (refinement 1) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_3021_85km_1 vals_CS_3021_85km_1 -r ../../setup/results_sparta/CS_3021/vals_85km_1.*.dat 

rm -rf vals_CS_3021_85km_2 
rm -rf vals_CS_3021_85km_2.pvd 
echo 'Converting results of CS_3021 at 85km (refinement 2) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_3021_85km_2 vals_CS_3021_85km_2 -r ../../setup/results_sparta/CS_3021/vals_85km_2.*.dat 

rm -rf vals_CS_3021_85km_3 
rm -rf vals_CS_3021_85km_3.pvd 
echo 'Converting results of CS_3021 at 85km (refinement 3) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_3021_85km_3 vals_CS_3021_85km_3 -r ../../setup/results_sparta/CS_3021/vals_85km_3.*.dat 

rm -rf vals_CS_3021_115km_0 
rm -rf vals_CS_3021_115km_0.pvd 
echo 'Converting results of CS_3021 at 115km (refinement 0) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_3021_115km_0 vals_CS_3021_115km_0 -r ../../setup/results_sparta/CS_3021/vals_115km_0.*.dat 

rm -rf vals_CS_3021_115km_1 
rm -rf vals_CS_3021_115km_1.pvd 
echo 'Converting results of CS_3021 at 115km (refinement 1) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_3021_115km_1 vals_CS_3021_115km_1 -r ../../setup/results_sparta/CS_3021/vals_115km_1.*.dat 

rm -rf vals_CS_3021_115km_2 
rm -rf vals_CS_3021_115km_2.pvd 
echo 'Converting results of CS_3021 at 115km (refinement 2) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_3021_115km_2 vals_CS_3021_115km_2 -r ../../setup/results_sparta/CS_3021/vals_115km_2.*.dat 

rm -rf vals_CS_3021_115km_3 
rm -rf vals_CS_3021_115km_3.pvd 
echo 'Converting results of CS_3021 at 115km (refinement 3) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_3021_115km_3 vals_CS_3021_115km_3 -r ../../setup/results_sparta/CS_3021/vals_115km_3.*.dat 

rm -rf vals_CS_3021_150km_0 
rm -rf vals_CS_3021_150km_0.pvd 
echo 'Converting results of CS_3021 at 150km (refinement 0) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_3021_150km_0 vals_CS_3021_150km_0 -r ../../setup/results_sparta/CS_3021/vals_150km_0.*.dat 

rm -rf vals_CS_3021_150km_1 
rm -rf vals_CS_3021_150km_1.pvd 
echo 'Converting results of CS_3021 at 150km (refinement 1) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_3021_150km_1 vals_CS_3021_150km_1 -r ../../setup/results_sparta/CS_3021/vals_150km_1.*.dat 

rm -rf vals_CS_3021_150km_2 
rm -rf vals_CS_3021_150km_2.pvd 
echo 'Converting results of CS_3021 at 150km (refinement 2) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_3021_150km_2 vals_CS_3021_150km_2 -r ../../setup/results_sparta/CS_3021/vals_150km_2.*.dat 

rm -rf vals_CS_3021_150km_3 
rm -rf vals_CS_3021_150km_3.pvd 
echo 'Converting results of CS_3021 at 150km (refinement 3) to ParaView...'
pvpython ../../tools/grid2paraview_original.py def/grid.CS_3021_150km_3 vals_CS_3021_150km_3 -r ../../setup/results_sparta/CS_3021/vals_150km_3.*.dat 
