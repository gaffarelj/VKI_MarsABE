import sys
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
import numpy as np
from MCD.parallel_mcd import parallel_mcd as PMCD
from netCDF4 import Dataset

mcd = PMCD()

h = 85e3
n = 0

for Ls in np.arange(0, 360.1, 30):                              # solar longitude [deg], between 0 and 360 deg, by 30 deg
    for loctime in np.arange(0, 24.1, 2):                       # local time [hr], between 0 and 24 Martian hrs, by 2 hrs
        for h in np.arange(85e3, 150e3, 3):                     # altitude [m], between 85km and 150km, by 3m
            for lat in np.arange(-90, 90.1, 3.75):              # latitude [deg], between -90 and 90, by 3.75 deg
                for lon in np.arange(-180, 180.1, 5.625):       # longitude [deg], between -180 and 180, by 5.625 deg
                    n += 1
print(n) # n = 9,937,352,880 (10 billion)