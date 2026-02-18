# Python3 script DensToIsoDens tool
# Purpose : convert water density to isotopic density
# Author : R. Guasch
# Date : 18/09/2024, updated 19/11/2025 : use IAPWS to get density at T/P conditions

from iapws import IAPWS97

def THMPT(Pressure, Temperature):
    """
    Parameters :
    ----------
        - Pressure: Pressure in Pascal
        - Temperature: Temperature in Kelvin

    Returns:
    ----------
        - Density: Density in Kg/m^3
        - Enthalpy: Enthalpy in J/Kg
        - Conductivity (k): Thermal conductivity in W/m/K
        - Viscosity (mu): Dynamic viscosity in Pa.s
        - Cp: Specific heat capacity in J/Kg/K
    """

    # Create an instance of the IAPWS97 class for water
    water = IAPWS97(P=Pressure/1e6, T=Temperature)

    # Extract properties from the water object
    Density = water.rho
    Enthalpy = water.h
    Conductivity = water.k
    Viscosity = water.mu
    Cp = water.cp

    return Density, Enthalpy*1e3, Conductivity, Viscosity, Cp

def THMTX(Temperature, Quality):
    """
    Get Density, enthalpy, conductivity, dynamic viscosity and specific heat at Temperature and Quality
    Parameters :
    ----------
        - Temperature: Temperature in Kelvin
        - Quality: Quality of the mixture (0 for liquid, 1 for vapor)
    Returns:
    ----------
        - Density (rho): Density in Kg/m^3 at 
        - Enthalpy (h): Enthalpy in J/Kg
        - Conductivity (k): Thermal conductivity in W/m/K
        - Viscosity (mu): Dynamic viscosity in Pa.s
        - Specific heat capacity (Cp) in J/Kg/K

    """
    # Create an instance of the IAPWS97 class for water
    print(f"THMTX : Temperature: {Temperature} K, Quality: {Quality}")
    water = IAPWS97(T=Temperature, x=Quality)

    # Extract properties from the water object
    Density = water.rho
    Enthalpy = water.h
    Conductivity = water.k
    Viscosity = water.mu
    Cp = water.cp
    return Density, Enthalpy*1e3, Conductivity, Viscosity, Cp


def DensToIsoDens(coolant_density):
    """
    density : density of the material in g/cm3
    isotopic_mass : mass of the isotope in g/mol
    """
    # Calculation of moderator data
    # AVOGADRO's number
    A = 6.022094E-1 # Normalizing by 10E-24 to obtain isotopic density in # / b*cm

    M_H2O = 15.9994 + 2.0*1.00794
    N_H2O = coolant_density*A/M_H2O
    N_O = N_H2O
    N_H = 2.0*N_H2O 
    return N_H, N_O

def densForPT(P , T):
    """Calculate the density for a given pressure and temperature.

    Args:
        P (float): The pressure in Pa.
        T (float): The temperature in K.

    Returns:
        float: The density in kg/m^3.
    """

    density, _, _, _, _ = THMPT(P, T)
    N_H, N_O = DensToIsoDens(density/1e3)  # Convert kg/m^3 to g/cm^3
    return N_H, N_O


def densForVoidFraction(T, void_fraction, N_H_H2O, N_O16_H2O):
    """
    Calculate the density for a given void fraction.
    Args:
        T (float): The temperature in K.
        void_fraction (float): The void fraction (between 0 and 1).
        N_H_H2O (float): The isotopic density of hydrogen in water (atoms/b.cm).
        N_O16_H2O (float): The isotopic density of oxygen in water (atoms/b.cm).
    Returns:
        tuple: The isotopic densities of hydrogen and oxygen (atoms/b.cm).
    """
    rho_l, _, _, _, _ = THMTX(T, 0.0)
    rho_vap, _, _, _, _ = THMTX(T, 1.0)
    r = rho_vap / rho_l
    fact = (1-void_fraction + void_fraction*r)
    rho_mixture = rho_l * (1.0 - void_fraction) + rho_vap * void_fraction
    N_H1, N_O1 = DensToIsoDens(rho_mixture/1e3)

    N_H = N_H_H2O * fact
    N_O = N_O16_H2O * fact
    return N_H, N_O, N_H1, N_O1




if __name__ == "__main__":
    # Example usage
    P = 7.2e6  # Pressure in Pa
    T = 559.0  # Temperature in K
    N_H, N_O = densForPT(P, T)
    print(f"For P={P} Pa and T={T} K:")
    print(f"N_H = {N_H:E} atoms/b.cm")
    print(f"N_O = {N_O:E} atoms/b.cm")
    ###
    N_H_H20 = 4.94546E-02  # Hydrogen in water
    N_O16_H2O = 2.47298E-02  # Oxygen in water
    void_fractions_to_test = [0.0, 0.4, 0.6, 0.8]

    for void_fraction in void_fractions_to_test:
        N_H_void, N_O_void, N_H1_void, N_O1_void = densForVoidFraction(T, void_fraction, N_H_H20, N_O16_H2O)
        print(f"For void_fraction={void_fraction}:")
        print(f"N_H = {N_H_void:E} atoms/b.cm")
        print(f"N_O = {N_O_void:E} atoms/b.cm")
        print(f"N_H (from mixture) = {N_H1_void:E} atoms/b.cm")
        print(f"N_O (from mixture) = {N_O1_void:E} atoms/b.cm")
    
    #for dens in [0.409331512518504]:
#    print(f"N_H = {DensToIsoDens(dens)[0]:E}, N_O = {DensToIsoDens(dens)[1]:E}")


