# Python3 script DensToIsoDens tool
# Purpose : convert water density to isotopic density
# Author : R. Guasch
# Date : 18/09/2024

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


#for dens in [0.7398, 0.6897]:
#    print(f"N_H = {DensToIsoDens(dens)[0]:E}, N_O = {DensToIsoDens(dens)[1]:E}")


