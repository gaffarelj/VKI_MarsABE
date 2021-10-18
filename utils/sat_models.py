from tudatpy.kernel.math import interpolators
import numpy as np


class satellite:
    def __init__(self, name, mass, Cd, S_ref=0.01, Cd_h=[85e3, 115e3, 150e3], \
        SA_areas=[0, 0, 0], SA_frac=0.7042, SA_eff=0.29, EPS_eff=0.89, verbose=False):
        """
        Satellite class, containing all relevant information for all analysed satellites configurations
        Inputs:
         * name (str): name of the satellite configuration
         * mass (float OR [float]*2): float OR list containing the constant satellite mass OR respectively the wet and dry mass of the satellite
         * Cd (float [float]*N): float OR list containing the drag coefficient of the satellite fixed OR at some altitudes
         * S_ref (float): reference surface area of the satellite
         * Cd_h ([float]*N): list at which the drag coefficients have been given
         * SA_areas ([float]*3): list containing the solar panel surface area of the satellite in each plane
         * SA_frac (float): fraction of the satellite areas that is actually the solar array
         * SA_eff (float): efficiency of the solar array
         * EPS_eff (float): efficiency of the Electrical Power System
        """
        # Save the satellite properties
        self.name = name
        self.area_x = SA_areas[0]*SA_frac
        self.area_y = SA_areas[1]*SA_frac
        self.area_z = SA_areas[2]*SA_frac
        if type(mass) != list:
            mass = [mass, mass]
        self.wet_mass = mass[0]
        self.dry_mass = mass[1]
        self.Cd_list = Cd
        self.Cd_h = Cd_h
        self.S_ref = S_ref
        self.SA_eff = SA_eff
        self.EPS_eff = EPS_eff
        self.h_warning = [verbose, verbose]
        
        if type(self.Cd_list) == list:
            # Create a cubic spline interpolator, capped at the boundaries
            interpolator_settings = interpolators.cubic_spline_interpolation(boundary_interpolation=interpolators.use_default_value)
            drag_dict = dict(zip(self.Cd_h, self.Cd_list))
            # Setup the drag interpolator
            self.drag_interpolator = interpolators.create_one_dimensional_interpolator(drag_dict, interpolator_settings)

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
            # Supress further warnings
            if h > max(self.Cd_h):
                self.h_warning[0] = False
            else:
                self.h_warning[1] = False
        Cd = self.drag_interpolator.interpolate(h)
        # If the interpolated Cd is out of the boundaries, take the average
        if Cd == 0:
            Cd = np.mean(self.Cd_list)
        return Cd

satellites = {
    "CS_0020": satellite("CS_0020", [4.29425, 3.64425], [2.83794, 2.74090, 2.80022], SA_areas=[0, 0.042426, 0.042426]),
    "CS_0021": satellite("CS_0021", [4.80225, 4.15225], [4.13116, 3.49154, 3.49964], SA_areas=[0.031058, 0.083343, 0.083343]),
    "CS_1020": satellite("CS_1020", [4.80225, 4.15225], [6.53680, 8.47435, 8.21728], SA_areas=[0, 0.102426, 0.042426]),
    "CS_1021": satellite("CS_1021", [5.31025, 4.66025], [5.33072, 4.29959, 4.24327], SA_areas=[0.031058, 0.143343, 0.083343]),
    "CS_2020": satellite("CS_2020", [5.31025, 4.66025], [7.37655, 9.09419, 8.75681], SA_areas=[0, 0.162426, 0.042426]),
    "CS_2021": satellite("CS_2021", [5.81825, 5.16825], [6.97704, 5.73235, 5.50560], SA_areas=[0.031058, 0.203343, 0.083343]),
    "CS_2120": satellite("CS_2120", [5.81825, 5.16825], [6.43496, 5.09859, 4.97959], SA_areas=[0, 0.282426, 0.042426]),
    "CS_3020": satellite("CS_3020", [6.32625, 5.67625], [8.61881, 9.91262, 9.48819], SA_areas=[0, 0.222426, 0.042426]),
    "CS_3021": satellite("CS_3021", [6.32625, 5.67625], [9.77352, 10.80229, 10.21743], SA_areas=[0.031058, 0.263343, 0.083343])
}