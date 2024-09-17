# Python3 -- PyGan script for generating FEM "full core" geometry
# Purpose : perform DONJON geometry generation and meshing for a single fuel assembly made of a single fuel rod.
# The geometry is generated with GEO module, meshed by TRIVAT: and material indexation is performed by MATEX:.
# Author : R. Guasch
# Date : 16/09/2024
import lifo
import cle2000

def INIT_DONJON(I_z, ene_groups):
    """
    pitch : real number representing the pitch of the fuel assembly / pincell
    height : real number representing the height of the active core / pincell
    I_z : integer representing the number of axial nodes in the fuel rod (Finite volumes for THM solution) 
    ene_groups : integer representing the number of energy groups for the neutron transport calculation   

    Tf_ini : real number initial fuel temperature
    Tc_ini : real number initial coolant temperature
    DensC_ini : real number initial coolant density
    power : real number initial reactor power
    """

    # Lifo
    myLifo=lifo.new()
    myLifo.push(I_z)
    myLifo.push(ene_groups)

    myLifo.lib()
    calcINIT = cle2000.new('INIT_DONJON',myLifo,1)
    calcINIT.exec()

    #recover
    #myLifo.lib()
    #pyTRACKS = myLifo.node(TRACK)
    #pyMATEX = myLifo.node(MaterialIndexation)
    
    return #pyTRACKS, pyMATEX