import numpy as np
import sys
sys.path.insert(0,"\\".join(sys.path[0].split("\\")[:-2]))
from tools import time_conversions as TC

class parallel_mcd:

    def __init__(self, default_inputs=True, load_on_init=True, load_parallel=True):
        # The following two arrays contain the starting time in Ls and Julian date at which the given module should be used
        self.limiting_Ls = np.arange(0, 330.01, 30)
        self.call_mcd_list = []
        self.n_species = (56, 65)
        self.Mars_R = 3389.5e3
        self.load_parallel = load_parallel
        if default_inputs:
            self.default_inputs()
        if load_parallel:
            self.load_mcd(charge_files=load_on_init)
        else:
            from MCD.fmcd import call_mcd as call_mcd
            self.call_mcd_list = [call_mcd]

    def default_inputs(self):
        # Load default inputs for the MCD.
        self.zkey = 1                    # xz is thus the radial distance from the centre of Mars
        self.xz = self.Mars_R+200e3      # [m], distance from the centre of Mars (200km altitude as default)
        self.xlon = 137.4                # [deg], East longitude
        self.xlat = -4.6                 # [deg], Latitude
        self.hireskey = 0                # use the lower resolution grid of 5.625x3.75 deg
        self.localtime = 12              # time of the day in hours
        self.dset = "/mnt/c/MCD/data/"   # path to the MCD data
        self.scena = 1                   # Climatology Scenario, solar EUV average conditions; use 2 for minimum solar conditions, and 3 for maximum
        self.perturkey = 0               # do not additional perturbations
        self.seedin = 1                  # (not) used for the random number generator
        self.gwlength = 0                # (not) used for gravity waves length
        # select which external variables to return
        self.extvarkeys = np.zeros(100)
        self.extvarkeys[self.n_species[0]:self.n_species[1]+1] = 1
        self.datekey = 1                 # 0 = Julian, 1 = Ls
        self.xdate = 0

    def load_mcd(self, charge_files=True):
        # Import the same Fortran interface compiled with different names in Python.
        # This is not an elegant solution, but a very efficient one.
        # This way, the required file can be loaded only once per module,
        # and data from the MCD can be obtained by switching modules depending on the time
        from MCD.fmcd_1 import call_mcd as call_mcd_1
        from MCD.fmcd_2 import call_mcd as call_mcd_2
        from MCD.fmcd_3 import call_mcd as call_mcd_3
        from MCD.fmcd_4 import call_mcd as call_mcd_4
        from MCD.fmcd_5 import call_mcd as call_mcd_5
        from MCD.fmcd_6 import call_mcd as call_mcd_6
        from MCD.fmcd_7 import call_mcd as call_mcd_7
        from MCD.fmcd_8 import call_mcd as call_mcd_8
        from MCD.fmcd_9 import call_mcd as call_mcd_9
        from MCD.fmcd_10 import call_mcd as call_mcd_10
        from MCD.fmcd_11 import call_mcd as call_mcd_11
        from MCD.fmcd_12 import call_mcd as call_mcd_12
        from MCD.fmcd_13 import call_mcd as call_mcd_13
        
        # Load the MCD corresponding to each Ls range (solar longitude)
        for i_module, xdate in enumerate(self.limiting_Ls):
            call_mcd = eval("call_mcd_%i" % (i_module+1))
            self.call_mcd_list.append(call_mcd)
            # Load the MCD data
            if charge_files:
                # Load all of the MCD files for one Martian year, in different modules (to keep them loaded)
                print("The MCD is loading for a full Martian year, please wait up to a few minutes...")
                call_mcd(self.zkey,self.xz,self.xlon,self.xlat,self.hireskey,self.datekey,xdate,\
                    self.localtime,self.dset,self.scena,self.perturkey,self.seedin,self.gwlength,self.extvarkeys)

    def call(self, Ls=None, localtime=None, lat=None, lon=None, h=None, print_results=False):
        # Call the MCD
        if Ls is not None: self.xdate = Ls
        if localtime is not None: self.localtime = localtime
        if lat is not None: self.xlat = lat
        if lon is not None: self.xlon = lon
        if h is not None: self.xz = self.Mars_R+h
        # Find what module should be used for the given solar longitude
        if self.load_parallel:
            i_module = np.where(self.limiting_Ls - self.xdate <= 15)[0][-1]
        else:
            i_module = 0
        call_mcd = self.call_mcd_list[i_module]
        self.pres,self.dens,self.temp,self.zonwind,self.merwind,self.meanvar,self.extvars,self.seedout,self.ier = \
            call_mcd(self.zkey,self.xz,self.xlon,self.xlat,self.hireskey,self.datekey,self.xdate,\
            self.localtime,self.dset,self.scena,self.perturkey,self.seedin,self.gwlength,self.extvarkeys)
        # Extract the volumetric ratio of each species in the atmosphere
        self.species_name = ["CO2", "N2", "Ar", "CO", "O", "O2", "O3", "H", "H2"] # note: He also accessible if required later, at index 77
        self.species_frac = []
        for n in range(*self.n_species):
            self.species_frac.append(self.extvars[n])
        self.species_dict = dict(zip(self.species_name, self.species_frac))

        # Print the results
        if print_results:
            print("Call MCD at Ls=%.2f deg and %.2f hrs, Lat=%.2f deg, Lon=%.2f deg, and h=%.2f km" % (self.xdate, self.localtime, self.xlat, self.xlon, (self.xz-self.Mars_R)/1e3))
            print("Density = %.5e [kg/m3]" % self.dens)
            print("Wind = %.5f E / %.5f N [m/s]" % (self.zonwind, self.merwind))
            print("Species [mol/mol] = %s" % self.species_dict)

    def density(self, h, lat, lon, time, update_inputs=True, time_is_JD=True, JD_Tudat=True):
        self.xz = self.Mars_R + h
        self.xlat = lat
        self.xlon = lon
        # TODO: check the date conversion
        if time_is_JD:
            Ls, Ds = TC.JD_to_Ls(time, JD_Tudat=JD_Tudat)
        else:
            Ls, Ds = time
        self.localtime = Ds % 1 * 24
        self.xdate = Ls
        self.call()
        return self.dens