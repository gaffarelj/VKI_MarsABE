import numpy as np
from numpy.core.shape_base import vstack
from scipy.interpolate import interpn
import sys
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from tools import time_conversions as TC


file_path = "/mnt/c/Mars-GRAM 2010/Release1.0_Nov10/txtFiles/tpdmsy11.txt"

class GRAM_atmo:

    def __init__(self):
        # Load the atmospheric data file. The columns are as follows:
        # Ls [deg]: solar longitude
        # Hgt [km]: altitude
        # Lat [def]: latitude
        # TA0K [K]: temperature
        # Pa0Nm2 [Pa]: pressure
        # DA0kgm3 [kg/m3]: density
        self.atmo_data = np.loadtxt(file_path, skiprows=1, usecols=(0, 1, 2, 3, 8, 13))

        # Setup list with all possible values for each parameter (Ls, h, lat)
        def possible_values(column):
            # Make a uniformely distributed range from the column data
            diff = np.array(column[1:] - column[:-1])
            # Find the step between the values
            dx = min(np.fabs(diff[diff != 0]))
            return np.arange(min(column), max(column)+dx/1e3, dx) # use dx/1e3 to make sure that the last value is also used
        self.possible_Ls = possible_values(self.atmo_data[:,0])
        self.possible_h = possible_values(self.atmo_data[:,1])
        self.possible_lat = possible_values(self.atmo_data[:,2])


        self.table_params = (self.possible_Ls, self.possible_h, self.possible_lat)
        self.table_values = np.empty((len(self.possible_Ls), len(self.possible_h), len(self.possible_lat), 3))
        i_t = 0
        for i, Ls in enumerate(self.possible_Ls):
            table_Ls = self.atmo_data[self.atmo_data[:,0]==Ls]
            for j, h in enumerate(self.possible_h):
                table_h = table_Ls[table_Ls[:,1]==h]
                for k, lat in enumerate(self.possible_lat):
                    self.table_values[i, j, k] = (table_h[k,3:])
                    i_t += 1
                    
    def get_density(self, h, _lon, lat, time, time_is_JD=True, JD_Tudat=True):
        if time_is_JD:
            Ls, Ds = TC.JD_to_Ls(time, JD_Tudat=JD_Tudat)
        else:
            Ls, Ds = time
        T, P, rho = interpn(self.table_params, self.table_values, [Ls, h/1e3, lat], bounds_error=True, fill_value=None)[0]
        print(h/1e3, _lon, lat, Ls)
        print(rho), input()
        if lat != 0:
            input()
        return rho