import sys
sys.path = [p for p in sys.path if p != ""]
while sys.path[0].split("/")[-1] != "VKI_MarsABE":
    sys.path.insert(0,"/".join(sys.path[0].split("/")[:-1]))
from tudatpy.kernel.math import interpolators
import numpy as np


class satellite:
    def __init__(self, name, mass, Cd, prop_mass=0, S_ref=0.01, Cd_h=[85e3, 115e3, 150e3], \
        SA_areas=[0, 0, 0], SA_frac=0.7042, SA_eff=0.29, EPS_eff=0.89, S_t=0, comp_ratio=1, verbose=False):
        """
        Satellite class, containing all relevant information for all analysed satellites configurations
        Inputs:
         * name (str): name of the satellite configuration
         * mass (float OR [float]*2): float OR list containing the constant satellite mass OR respectively the wet and dry mass of the satellite
         * Cd (float OR [float]*N): float OR list containing the drag coefficient of the satellite fixed OR at some altitudes
         * S_ref (float): reference surface area of the satellite
         * Cd_h ([float]*N): list at which the drag coefficients have been given
         * SA_areas ([float]*3): list containing the solar panel surface area of the satellite in each plane
         * SA_frac (float): fraction of the satellite areas that is actually the solar array
         * SA_eff (float): efficiency of the solar array
         * EPS_eff (float): efficiency of the Electrical Power System
         * S_t (float): end of ABE inlet area (throat)
         * comp_ratio (float OR [float]*N): compression ratio between the free stream density and the throat density; either fixed value or depending on the altitude from CD_h
        """
        # Save the satellite properties
        self.name = name
        self.area_x = SA_areas[0]*SA_frac
        self.area_y = SA_areas[1]*SA_frac
        self.area_z = SA_areas[2]*SA_frac
        self.wet_mass = mass
        self.dry_mass = mass - prop_mass
        self.Cd_list = Cd
        self.Cd_h = Cd_h
        self.S_ref = S_ref
        self.SA_eff = SA_eff
        self.EPS_eff = EPS_eff
        self.h_warning = [verbose, verbose]
        self.S_t = S_t
        self.comp_ratio_list = comp_ratio
        
        if type(self.Cd_list) == list:
            # Create a cubic spline interpolator, capped at the boundaries
            interpolator_settings = interpolators.cubic_spline_interpolation(boundary_interpolation=interpolators.throw_exception_at_boundary)
            drag_dict = dict(zip(self.Cd_h, self.Cd_list))
            # Setup the drag interpolator
            self.drag_interpolator = interpolators.create_one_dimensional_scalar_interpolator(drag_dict, interpolator_settings)

        if type(self.comp_ratio_list) == list:
            # Create a cubic spline interpolator, capped at the boundaries
            interpolator_settings = interpolators.cubic_spline_interpolation(boundary_interpolation=interpolators.throw_exception_at_boundary)
            comp_ratio_dict = dict(zip(self.Cd_h, self.comp_ratio_list))
            # Setup the drag interpolator
            self.comp_ratio_interpolator = interpolators.create_one_dimensional_scalar_interpolator(comp_ratio_dict, interpolator_settings)

    def __str__(self):
        # Return the satellite represented as a string
        if type(self.Cd_list) == list:
            cd_l = (min(self.Cd_list), max(self.Cd_list))
            cd_str = "Cd ranging in (%.2f; %.2f)" % cd_l
        else:
            cd_l = (self.Cd_list, self.Cd_list)
            cd_str = "constant Cd of %.2f" % self.Cd_list
        b_c = [self.dry_mass/(cd_l[1]*self.S_ref), self.wet_mass/(cd_l[0]*self.S_ref)]
        return "Satellite '%s' with wet mass of %.2f kg (dry of %.2f kg), reference area of %.4f m2, %s, ballistic coefficient ranging in (%.2f; %.2f)" % \
            (self.name, self.wet_mass, self.dry_mass, self.S_ref, cd_str, *b_c)

    def __repr__(self):
        return self.name

    def get_cd(self, h):
        if type(self.Cd_list) != list:
            return self.Cd_list
        # Return the interpolated drag coefficient at the given altitude h [m]
        if (h > max(self.Cd_h) and self.h_warning[0]) or (h < min(self.Cd_h) and self.h_warning[1]):
            print("Warning: the specified altitude of %i km is outside the measured range of (%i; %i) km." % (h/1e3, min(self.Cd_h)/1e3, max(self.Cd_h)/1e3))
            # Suppress further warnings
            if h > max(self.Cd_h):
                self.h_warning[0] = False
            else:
                self.h_warning[1] = False
        # If the interpolated Cd is out of the boundaries, take the average
        try:
            Cd = self.drag_interpolator.interpolate(h)
        except RuntimeError:
            Cd = np.mean(self.Cd_list)
        return Cd

    def get_comp_ratio(self, h):
        if type(self.comp_ratio_list) != list:
            return self.comp_ratio_list
        # Return the interpolated drag coefficient at the given altitude h [m]
        # If the interpolated Cd is out of the boundaries, take the average
        try:
            comp_ratio = self.comp_ratio_interpolator.interpolate(h)
        except RuntimeError:
            comp_ratio = np.mean(self.comp_ratio_list)
        return comp_ratio

# Define atmosphere-breathing satellite properties
satellites = {
    "CS_0021": satellite("CS_0021", 2.93225, [1.277012, 1.803986, 1.666586], S_ref=0.0104, S_t=5e-4, \
        SA_areas=[0.031058, 0.083343, 0.083343], comp_ratio=[37.82965, 50.49326, 60.37834]),
    "CS_1021": satellite("CS_1021", 3.44025, [1.469716, 1.897386, 1.761229], S_ref=0.0108, S_t=5e-4, \
        SA_areas=[0.031058, 0.143343, 0.083343], comp_ratio=[35.10988, 50.17104, 60.29815]),
    "CS_2021": satellite("CS_2021", 3.94825, [1.721637, 2.006168, 1.872265], S_ref=0.0112, S_t=5e-4, \
        SA_areas=[0.031058, 0.203343, 0.083343], comp_ratio=[35.27233, 49.77842, 59.51885]),
    "CS_2120": satellite("CS_2120", 3.94825, [5.776334, 3.991432, 4.195638], S_ref=0.041858285, S_t=5e-4, \
        SA_areas=[0, 0.282426, 0.042426], comp_ratio=[36.17590, 50.10605, 58.88298]),
    "CS_3021": satellite("CS_3021", 4.45625, [1.959728, 2.141222, 2.007698], S_ref=0.042258285, S_t=5e-4, \
        SA_areas=[0.031058, 0.263343, 0.083343], comp_ratio=[42.28136, 58.14033, 64.84591])
}

# Define the same satellites, but with mass difference to account for the lack of atmosphere-breathing inlet, and the addition of a Xenon propellant tank
satellites_with_tank = satellites.copy()
for name, sat in satellites_with_tank.items():
    sat.dry_mass -= 0.22
    sat.wet_mass += 0.43

# Plot satellite the CD vs altitude
if False:
    from matplotlib import pyplot as plt
    plt.rcParams.update({'font.size': 13, 'figure.figsize': (10.5, 7), 'savefig.format': 'pdf'})
    for s_name, sat in satellites.items():
        plt.plot(np.array(sat.Cd_h)/1e3, sat.Cd_list, label=s_name, linestyle="dotted", marker="o")
    plt.xlabel("Altitude [km]"), plt.ylabel("Drag coefficient [-]")
    plt.xticks([85, 115, 150])
    plt.legend(), plt.grid(), plt.tight_layout()
    plt.savefig(sys.path[0]+"/figures/sat_Cd_h.pdf")