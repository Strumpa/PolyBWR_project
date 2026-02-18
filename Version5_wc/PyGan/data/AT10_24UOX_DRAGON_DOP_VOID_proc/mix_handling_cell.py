#####################################################################
#                                                                   #
# Description : PyGan script for BWR simulation with DRAGON5        #
# Author      : R.Guasch                                            #
# Date        : 2025/06/26                                          #
# Purpose     : MIX definition for pin cell, as a function of :     #
#                               - TFuel,                            #
#                               - TCool,                            #
#                               - Dcool,                            #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def dens_to_iso(coolant_density):
    """
    Parameters
    ----------
    coolant_density : float
       Density of water in kg/m^3

    Returns
    -------
    N_H : float
        Isotopic density of Hydrgen in water in atom / b-cm
    N_O : float
        Isotopic density of Oxygen in water in atom / b-cm
    """
    # Convert density to g/cm^3
    coolant_density = coolant_density/1000
    # Calculation of moderator data
    # AVOGADRO's number
    A = 6.022094E-1 # Normalizing by 10E-24 to obtain isotopic density in # / b*cm

    M_H2O = 15.9994 + 2.0*1.00794
    N_H2O = coolant_density*A/M_H2O
    N_O = N_H2O
    N_H = 2.0*N_H2O 
    return float(f"{N_H:4E}"), float(f"{N_O:4E}")


def LIB_24UOX(Library,ssh_option,TFuel,TCool,Dcool):
    """
    Parameters :
    -------------
    Library : str
        Name of the cross sectional libary to be used for DRAGON5 calculations ie draglib name given in .access symbolic link
    ssh_option : str
        Name of the selected self-shielding method, only PT and RSE are supported here
    TFuel : 
        Effective fuel temperature (K)
    TCool :
        Coolant mixture temperature (K)
    DCool : 
        Coolant mixture density (kg/m^3)

    Returns :
    -----------
    pyLIB : lcm object 
        Microscopic library to be self-shielded
    N_H : float
        Number density of Hydrogen in cooling water, in atom/b-cm
    N_O : float
        Number density of Oxygen in cooling water, in atom/b-cm

    """

    # Compute N_H and N_O from DCool
    N_H, N_O = dens_to_iso(Dcool)
    # Lifo
    myLifo=lifo.new()
    myLifo.pushEmpty("LIBRARY", "LCM")
    myLifo.push(Library)
    myLifo.push(ssh_option)
    myLifo.push(TFuel)
    myLifo.push(TCool)
    myLifo.push(N_H)
    myLifo.push(N_O)
    myLifo.lib()

    # Execution
    mixBWR = cle2000.new('MIX_C',myLifo,1)
    mixBWR.exec()

    # Recover
    myLifo.lib()
    pyLIB = myLifo.node("LIBRARY")

    # Clear stack before next execution
    while myLifo.getMax() > 0:
      myLifo.pop()

    return pyLIB, N_H, N_O
