import sys
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
import numpy as np
import time
from MCD.parallel_mcd import parallel_mcd as PMCD
from netCDF4 import Dataset

# Load the Mars Climate Database
mcd = PMCD()

# List all of the variables
Ls_list = np.arange(0, 330.1, 30)           # solar longitude [deg], between 0 and 360 deg, by 30 deg
localtime_list = np.arange(0, 22.1, 2)      # local time [hr], between 0 and 24 Martian hrs, by 2 hrs
h_list = np.arange(85000, 150001.1, 3)      # altitude [m], between 85km and 150km, by 3m
lat_list = np.arange(-90, 88, 3.75)         # latitude [deg], between -90 and 90, by 3.75 deg
lon_list = np.arange(-180, 175, 5.625)      # longitude [deg], between -180 and 180, by 5.625 deg

# Compute the total number of elements to save
N = len(Ls_list)*len(localtime_list)*len(h_list)*len(lat_list)*len(lon_list)
n = 1
print("Total elements to save: %i (%.1fB)" % (N, N/1e9))
times = []

for i_L, Ls in enumerate(Ls_list):
    # Create the nc file    
    nc_data = Dataset("MCD/nc/Ls_%03d.nc" % Ls, "w", format="NETCDF4")

    # Create the dimensions
    Ls_dim          = nc_data.createDimension("Ls", 12)
    localtime_dim   = nc_data.createDimension("localtime", 12)
    latitude_dim    = nc_data.createDimension("latitude", 48)
    longitude_dim   = nc_data.createDimension("longitude", 64)
    altitude_dim    = nc_data.createDimension("altitude", 21668)

    # Create the variables
    Ls_s        = nc_data.createVariable("Ls","f4",("Ls",), least_significant_digit=4,zlib=True)
    Ls_s.units = "deg"
    Ls_s[:] = Ls_list
    localtime_s = nc_data.createVariable("localtime","f4",("localtime",), least_significant_digit=2,zlib=True)
    localtime_s.units = "hr"
    localtime_s[:] = localtime_list
    latitude_s  = nc_data.createVariable("latitude","f4",("latitude",), least_significant_digit=3,zlib=True)
    latitude_s.units = "deg"
    latitude_s[:] = lat_list
    longitude_s = nc_data.createVariable("longitude","f4",("longitude",), least_significant_digit=3,zlib=True)
    longitude_s.units = "deg"
    longitude_s[:] = lon_list
    altitude_s  = nc_data.createVariable("altitude","f4",("altitude",), least_significant_digit=6,zlib=True)
    altitude_s.units = "m"
    altitude_s[:] = h_list

    density     = nc_data.createVariable("density","f4",("Ls","localtime","latitude","longitude","altitude"), least_significant_digit=5,zlib=True)
    density.units = "kg/m3"
    temperature = nc_data.createVariable("temperature","f4",("Ls","localtime","latitude","longitude","altitude"), least_significant_digit=5,zlib=True)
    temperature.units = "K"
    pressure    = nc_data.createVariable("pressure","f4",("Ls","localtime","latitude","longitude","altitude"), least_significant_digit=5,zlib=True)
    pressure.units = "Pa"
    wind_v      = nc_data.createVariable("wind_v","f4",("Ls","localtime","latitude","longitude","altitude"), least_significant_digit=5,zlib=True)
    wind_v.units = "m/s"
    wind_N      = nc_data.createVariable("wind_N","f4",("Ls","localtime","latitude","longitude","altitude"), least_significant_digit=5,zlib=True)
    wind_N.units = "m/s"
    wind_E      = nc_data.createVariable("wind_E","f4",("Ls","localtime","latitude","longitude","altitude"), least_significant_digit=5,zlib=True)
    wind_E.units = "m/s"
    fCO2        = nc_data.createVariable("fCO2","f4",("Ls","localtime","latitude","longitude","altitude"), least_significant_digit=5,zlib=True)
    fCO2.units = "mol/mol"
    fN2         = nc_data.createVariable("fN2","f4",("Ls","localtime","latitude","longitude","altitude"), least_significant_digit=5,zlib=True)
    fN2.units = "mol/mol"
    fAr         = nc_data.createVariable("fAr","f4",("Ls","localtime","latitude","longitude","altitude"), least_significant_digit=5,zlib=True)
    fAr.units = "mol/mol"
    fCO         = nc_data.createVariable("fCO","f4",("Ls","localtime","latitude","longitude","altitude"), least_significant_digit=5,zlib=True)
    fCO.units = "mol/mol"
    fO          = nc_data.createVariable("fO","f4",("Ls","localtime","latitude","longitude","altitude"), least_significant_digit=5,zlib=True)
    fO.units = "mol/mol"
    fO2         = nc_data.createVariable("fO2","f4",("Ls","localtime","latitude","longitude","altitude"), least_significant_digit=5,zlib=True)
    fO2.units = "mol/mol"

    # Loop trough all variables
    for i_l, localtime in enumerate(localtime_list):
        for i_h, h in enumerate(h_list):
            for i_lat, lat in enumerate(lat_list):
                for i_lon, lon in enumerate(lon_list):
                    # Show progress
                    if n % 1e2 == 0:
                        if len(times) < 50:
                            times.append(time.time())
                            remaining = 0
                        else:
                            dts = np.array(times)[1:] - np.array(times)[:-1]
                            avg_dt = np.mean(dts)
                            times.pop(0), times.append(time.time())
                            remaining = avg_dt * (N - n)
                        print("Extracting data point %09d/%i. Estimated %.2f hrs remaining." % (n, N, remaining/3600), end="\r")
                    n += 1
                    # Call the MCD
                    mcd.call(Ls=Ls, localtime=localtime, lat=lat, lon=lon, h=h)
                    # Save data
                    density[i_L, i_l, i_lat, i_lon, i_h]        = mcd.dens
                    temperature[i_L, i_l, i_lat, i_lon, i_h]    = mcd.temp
                    pressure[i_L, i_l, i_lat, i_lon, i_h]       = mcd.pres
                    wind_v[i_L, i_l, i_lat, i_lon, i_h]         = mcd.vertwind
                    wind_N[i_L, i_l, i_lat, i_lon, i_h]         = mcd.merwind
                    wind_E[i_L, i_l, i_lat, i_lon, i_h]         = mcd.zonwind
                    fCO2[i_L, i_l, i_lat, i_lon, i_h]           = mcd.species_frac[0]
                    fN2[i_L, i_l, i_lat, i_lon, i_h]            = mcd.species_frac[1]
                    fAr[i_L, i_l, i_lat, i_lon, i_h]            = mcd.species_frac[2]
                    fCO[i_L, i_l, i_lat, i_lon, i_h]            = mcd.species_frac[3]
                    fO[i_L, i_l, i_lat, i_lon, i_h]             = mcd.species_frac[4]
                    fO2[i_L, i_l, i_lat, i_lon, i_h]            = mcd.species_frac[5]
                

    nc_data.close()
    print("/n Done with Ls = %i deg" % Ls)