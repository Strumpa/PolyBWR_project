# Version: 2021.09.22
# Author : Clément HUET, Raphaël Guasch
from multiPhysics_proc.THM import Version5_THM_prototype
from multiPhysics_proc.THM import plotting
from iapws import IAPWS97
import numpy as np

## User choice:
solveConduction = True
zPlotting = [0.8]

########## Thermal hydraulics parameters ##########
## Fluid parameters
# T_inlet, T_outlet = 270, 287 Celcius
# Nominal coolant flow rate = 1530 kg/s
# Nominal operating pressure = 7.2 MPa (abs)
hInlet = 
pOutlet =  
qFlow = 

## Geometric parameters
canalType = "square"
waterRadius = 1.295e-2 # m ATRIUM10 pincell pitch
fuelRadius = 0.4435e-2 # m : fuel rod radius
gapRadius = 0.4520e-2 # m : expansion gap radius : "void" between fuel and clad - equivalent to inner clad radius
cladRadius = 0.5140e-2 # m : clad external radius
height = 3.8 # m : height : active core height in BWRX-300 SMR

## Meshing parameters:
If = 8
I1 = 3
Iz1 = 10 # number of control volumes in the axial direction

## Thermalhydraulics correlation
voidFractionCorrel = "HEM1"
frfaccorel = "base"
P2Pcorel = "base"

############ Nuclear Parameters ###########
## Fission parameters
# specific power = 38.6 W/g
qFiss = 

## Material parameters
kFuel = 4.18 # W/m.K, TECHNICAL REPORTS SERIES No. 59 : Thermal Conductivity of Uranium Dioxide, IAEA, VIENNA, 1966
Hgap = 10000 
kClad = 21.5 # W/m.K, Thermal Conductivity of Zircaloy-2 (as used in BWRX-300) according to https://www.matweb.com/search/datasheet.aspx?MatGUID=eb1dad5ce1ad4a1f9e92f86d5b44740d
# k_Zircaloy-4 = 21.6 W/m.K too so check for ATRIUM-10 clad material but should have the same thermal conductivity
########## Algorithm parameters ###########
nIter = 1000
tol = 1e-4
underRelaxationFactor = 0.5

########## Fields of the problem ##########
TeffFuel = []
Twater = []
rho = []
Qfiss = []

## Fonction used for the convergence algorithm
def underRelaxation(Field, OldField, underRelaxationFactor):
    return underRelaxationFactor*Field + (1-underRelaxationFactor)*OldField

def convergence(Field, OldField, tol):
    return np.abs(Field-OldField) < tol

## Initial thermal hydraulic resolution
case1 = Version5_THM_prototype("Testing THM Prototype", canalType, waterRadius, fuelRadius, gapRadius, cladRadius,
                            height, hInlet, pOutlet, qFlow, qFiss, kFuel, Hgap, kClad, Iz1, If, I1, zPlotting, 
                            solveConduction, dt=0, t_tot=0, frfaccorel, P2Pcorel, voidFractionCorrel)

## MultiPhysics resolution
for i in range(nIter):
    
    ################## Nuclear part ##################

    ############# Thermalhydraulic part ##############
    case1 = Version5_THM_prototype("Testing THM Prototype", canalType, waterRadius, fuelRadius, gapRadius, cladRadius,
                            height, hInlet, pOutlet, qFlow, qFiss, kFuel, Hgap, kClad, Iz1, If, I1, zPlotting, 
                            solveConduction, dt=0, t_tot=0, frfaccorel, P2Pcorel, voidFractionCorrel)    ##### qFiss to be updated

    TeffTEMP, TwaterTEMP, rhoTEMP = case1.get_nuclear_parameters()
    TeffFuel.append(TeffTEMP)
    Twater.append(TwaterTEMP)
    rho.append(rhoTEMP)

    ############## Under relaxation #################
    TeffFuel[-1] = underRelaxation(TeffFuel[-1], TeffFuel[-2], underRelaxationFactor)
    Twater[-1] = underRelaxation(Twater[-1], Twater[-2], underRelaxationFactor)
    rho[-1] = underRelaxation(rho[-1], rho[-2], underRelaxationFactor)

    ############## Convergence test #################
    if convergence(TeffFuel[-1], TeffFuel[-2], tol) and convergence(Twater[-1], Twater[-2], tol) and convergence(rho[-1], rho[-2], tol):
        print("Convergence reached after ", i, " iterations")
        break

    if i == nIter-1:
        print("Convergence not reached after ", i, " iterations")
    


    
