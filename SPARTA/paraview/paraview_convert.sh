#!/bin/sh
cd surf
rm -rf *
pvpython ../../tools/surf2paraview.py ../../setup/data/data.CS_0020 CS_0020 
pvpython ../../tools/surf2paraview.py ../../setup/data/data.CS_0021 CS_0021 
pvpython ../../tools/surf2paraview.py ../../setup/data/data.CS_1020 CS_1020 
pvpython ../../tools/surf2paraview.py ../../setup/data/data.CS_1021 CS_1021 
pvpython ../../tools/surf2paraview.py ../../setup/data/data.CS_2020 CS_2020 
pvpython ../../tools/surf2paraview.py ../../setup/data/data.CS_2021 CS_2021 
pvpython ../../tools/surf2paraview.py ../../setup/data/data.CS_2120 CS_2120 
pvpython ../../tools/surf2paraview.py ../../setup/data/data.CS_3020 CS_3020 
pvpython ../../tools/surf2paraview.py ../../setup/data/data.CS_3021 CS_3021 
cd ../grid
rm -rf CS_0020_85km 
rm -rf CS_0020_85km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_0020_85km CS_0020_85km -r ../../setup/results_sparta/CS_0020/npart_85km.*.gz 
rm -rf CS_0020_115km 
rm -rf CS_0020_115km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_0020_115km CS_0020_115km -r ../../setup/results_sparta/CS_0020/npart_115km.*.gz 
rm -rf CS_0020_150km 
rm -rf CS_0020_150km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_0020_150km CS_0020_150km -r ../../setup/results_sparta/CS_0020/npart_150km.*.gz 
rm -rf CS_0021_85km 
rm -rf CS_0021_85km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_0021_85km CS_0021_85km -r ../../setup/results_sparta/CS_0021/npart_85km.*.gz 
rm -rf CS_0021_115km 
rm -rf CS_0021_115km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_0021_115km CS_0021_115km -r ../../setup/results_sparta/CS_0021/npart_115km.*.gz 
rm -rf CS_0021_150km 
rm -rf CS_0021_150km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_0021_150km CS_0021_150km -r ../../setup/results_sparta/CS_0021/npart_150km.*.gz 
rm -rf CS_1020_85km 
rm -rf CS_1020_85km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_1020_85km CS_1020_85km -r ../../setup/results_sparta/CS_1020/npart_85km.*.gz 
rm -rf CS_1020_115km 
rm -rf CS_1020_115km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_1020_115km CS_1020_115km -r ../../setup/results_sparta/CS_1020/npart_115km.*.gz 
rm -rf CS_1020_150km 
rm -rf CS_1020_150km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_1020_150km CS_1020_150km -r ../../setup/results_sparta/CS_1020/npart_150km.*.gz 
rm -rf CS_1021_85km 
rm -rf CS_1021_85km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_1021_85km CS_1021_85km -r ../../setup/results_sparta/CS_1021/npart_85km.*.gz 
rm -rf CS_1021_115km 
rm -rf CS_1021_115km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_1021_115km CS_1021_115km -r ../../setup/results_sparta/CS_1021/npart_115km.*.gz 
rm -rf CS_1021_150km 
rm -rf CS_1021_150km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_1021_150km CS_1021_150km -r ../../setup/results_sparta/CS_1021/npart_150km.*.gz 
rm -rf CS_2020_85km 
rm -rf CS_2020_85km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_2020_85km CS_2020_85km -r ../../setup/results_sparta/CS_2020/npart_85km.*.gz 
rm -rf CS_2020_115km 
rm -rf CS_2020_115km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_2020_115km CS_2020_115km -r ../../setup/results_sparta/CS_2020/npart_115km.*.gz 
rm -rf CS_2020_150km 
rm -rf CS_2020_150km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_2020_150km CS_2020_150km -r ../../setup/results_sparta/CS_2020/npart_150km.*.gz 
rm -rf CS_2021_85km 
rm -rf CS_2021_85km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_2021_85km CS_2021_85km -r ../../setup/results_sparta/CS_2021/npart_85km.*.gz 
rm -rf CS_2021_115km 
rm -rf CS_2021_115km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_2021_115km CS_2021_115km -r ../../setup/results_sparta/CS_2021/npart_115km.*.gz 
rm -rf CS_2021_150km 
rm -rf CS_2021_150km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_2021_150km CS_2021_150km -r ../../setup/results_sparta/CS_2021/npart_150km.*.gz 
rm -rf CS_2120_85km 
rm -rf CS_2120_85km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_2120_85km CS_2120_85km -r ../../setup/results_sparta/CS_2120/npart_85km.*.gz 
rm -rf CS_2120_115km 
rm -rf CS_2120_115km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_2120_115km CS_2120_115km -r ../../setup/results_sparta/CS_2120/npart_115km.*.gz 
rm -rf CS_2120_150km 
rm -rf CS_2120_150km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_2120_150km CS_2120_150km -r ../../setup/results_sparta/CS_2120/npart_150km.*.gz 
rm -rf CS_3020_85km 
rm -rf CS_3020_85km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_3020_85km CS_3020_85km -r ../../setup/results_sparta/CS_3020/npart_85km.*.gz 
rm -rf CS_3020_115km 
rm -rf CS_3020_115km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_3020_115km CS_3020_115km -r ../../setup/results_sparta/CS_3020/npart_115km.*.gz 
rm -rf CS_3020_150km 
rm -rf CS_3020_150km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_3020_150km CS_3020_150km -r ../../setup/results_sparta/CS_3020/npart_150km.*.gz 
rm -rf CS_3021_85km 
rm -rf CS_3021_85km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_3021_85km CS_3021_85km -r ../../setup/results_sparta/CS_3021/npart_85km.*.gz 
rm -rf CS_3021_115km 
rm -rf CS_3021_115km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_3021_115km CS_3021_115km -r ../../setup/results_sparta/CS_3021/npart_115km.*.gz 
rm -rf CS_3021_150km 
rm -rf CS_3021_150km.pvd 
pvpython ../../tools/grid2paraview.py grid.CS_3021_150km CS_3021_150km -r ../../setup/results_sparta/CS_3021/npart_150km.*.gz 
