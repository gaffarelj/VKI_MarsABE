#!/bin/sh
module load openmpi
mpirun -np 16 spa_ < in.CS_0021_85km | tee ../results_sparta/CS_0021/stats_85km.dat
