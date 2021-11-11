#!/bin/sh
module load openmpi
mpirun -np 16 spa_ < in.CS_0021_85km | tee ../results_sparta/CS_0021/stats_85km.dat
mpirun -np 16 spa_ < in.CS_0021_115km | tee ../results_sparta/CS_0021/stats_115km.dat
mpirun -np 16 spa_ < in.CS_0021_150km | tee ../results_sparta/CS_0021/stats_150km.dat
mpirun -np 16 spa_ < in.CS_1021_115km | tee ../results_sparta/CS_1021/stats_115km.dat
mpirun -np 16 spa_ < in.CS_1021_150km | tee ../results_sparta/CS_1021/stats_150km.dat
mpirun -np 16 spa_ < in.CS_2021_115km | tee ../results_sparta/CS_2021/stats_115km.dat
mpirun -np 16 spa_ < in.CS_2021_150km | tee ../results_sparta/CS_2021/stats_150km.dat
mpirun -np 16 spa_ < in.CS_2120_85km | tee ../results_sparta/CS_2120/stats_85km.dat
mpirun -np 16 spa_ < in.CS_2120_115km | tee ../results_sparta/CS_2120/stats_115km.dat
mpirun -np 16 spa_ < in.CS_2120_150km | tee ../results_sparta/CS_2120/stats_150km.dat
mpirun -np 16 spa_ < in.CS_3021_85km | tee ../results_sparta/CS_3021/stats_85km.dat
mpirun -np 16 spa_ < in.CS_3021_115km | tee ../results_sparta/CS_3021/stats_115km.dat
mpirun -np 16 spa_ < in.CS_3021_150km | tee ../results_sparta/CS_3021/stats_150km.dat
