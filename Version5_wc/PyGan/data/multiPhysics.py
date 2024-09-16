# Version: 2021.09.22
# Author : Clément HUET, Raphaêl Guasch
from multiPhysics_proc.THM import Version5_THM_prototype
from multiPhysics_proc.THM import plotting
from iapws import IAPWS97
import numpy as np

## User choice:
solveConduction = True
zPlotting = [0.8]

########## Thermal hydraulics parameters ##########
## Fluid parameters
hInlet = 
pOutlet = 
qFlow = 

## Geometric parameters
canalType = "circular"
waterRadius = 
fuelRadius = 
gapRadius = 
cladRadius = 
height = 

## Meshing parameters:
If = 8
I1 = 3
Iz1 = 10

## Thermalhydraulics correlation
voidFractionCorrel = "HEM1"
frfaccorel = "base"
P2Pcorel = "base"

############ Nuclear Parameters ###########
## Fission parameters
qFiss = 

## Material parameters
kFuel = 
Hgap = 
kClad = 

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
    


    
