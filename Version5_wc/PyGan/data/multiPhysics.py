## PyGan script to couple DONJON procedures with THM solution
# Date : 20/09/2024
# Author : Clément HUET, Raphaël GUASCH
# Purpose : test and validate neutronics and thermalhydraulics coupling on a single BWR pincell
from THM_main import Version5_THM_prototype as THM_prototype
from THM_main import plotting
from iapws import IAPWS97
import numpy as np
import os, shutil
import lifo
import lcm
import cle2000
from assertS import *


########## Begin functions declaration ##########
## Fonction used for the convergence algorithm
def underRelaxation(Field, OldField, underRelaxationFactor):
    return underRelaxationFactor*Field + (1-underRelaxationFactor)*OldField

def convergence(Field, OldField, tol):
    return np.abs(Field-OldField) < tol

def guessAxialPowerShape(amplitude, Iz, height, Fuel_volume):
    """
    Amplitude : float : amplitude of the axial power shape (W)
    Iz : int : number of control volumes in the axial direction
    height : float : height of the fuel rod (m)
    Fuel_volume : float : volume of the fuel rod (m3)
    return : np.array : axial power shape with a sine shape units (W/m3)
                        --> corresponds to the power density in each control volume
    """
    z_space = np.linspace(0, height, Iz)
    return (amplitude/(Fuel_volume))*np.sin(z_space*np.pi/height)



######## End functions declaration ##########


########## User input ##########

solveConduction = True
zPlotting = [0.8]

########## Thermal hydraulics parameters ##########
## Geometric parameters
canalType = "square"
waterRadius = 1.295e-2 # m ATRIUM10 pincell pitch
fuelRadius = 0.4435e-2 # m : fuel rod radius
gapRadius = 0.4520e-2 # m : expansion gap radius : "void" between fuel and clad - equivalent to inner clad radius
cladRadius = 0.5140e-2 # m : clad external radius
height = 3.8 # m : height : active core height in BWRX-300 SMR


## Fluid parameters
tInlet = 270 + 273.15 # K
# T_inlet, T_outlet = 270, 287 Celcius
pOutlet =  7.2e6 # Pa
pressureDrop = 186737 #Pa/m
falsePInlet = pOutlet - height * pressureDrop
rhoInlet = IAPWS97(T = tInlet, P = falsePInlet*10**(-6)).rho #kg/m3
flowArea = waterRadius ** 2 - np.pi * cladRadius ** 2

# Nominal coolant flow rate = 1530 kg/s
# Nominal operating pressure = 7.2 MPa (abs)
#hInlet =  # to fill

#qFlow =  # to fill
massFlowRate = 1530  / (200*91)  # kg/s

## Additional parameters needed for the calculation
solveConduction = True
volumic_mass_U = 19000 # kg/m3
Fuel_volume = np.pi*fuelRadius**2*height # m3
Fuel_mass = Fuel_volume*volumic_mass_U # kg

## Meshing parameters:
If = 8
I1 = 3
Iz1 = 20 # number of control volumes in the axial direction

## Thermalhydraulics correlation
voidFractionCorrel = "HEM1"
frfaccorel = "base"
P2Pcorel = "base"
numericalMethod = "FVM"

############ Nuclear Parameters ###########
## Fission parameters
# specific power = 38.6 W/g
specificPower = 38.6 # W/g
PFiss = specificPower*Fuel_mass*1000 # W

#qFiss = PFiss/Fuel_volume # W/m3

## Material parameters
kFuel = 4.18 # W/m.K, TECHNICAL REPORTS SERIES No. 59 : Thermal Conductivity of Uranium Dioxide, IAEA, VIENNA, 1966
Hgap = 10000 
kClad = 21.5 # W/m.K, Thermal Conductivity of Zircaloy-2 (as used in BWRX-300) according to https://www.matweb.com/search/datasheet.aspx?MatGUID=eb1dad5ce1ad4a1f9e92f86d5b44740d
# k_Zircaloy-4 = 21.6 W/m.K too so check for ATRIUM-10 clad material but should have the same thermal conductivity
########## Algorithm parameters ###########
nIter = 1000
tol_TH = 1e-4 # Convergence criterion for the TH solution
tol_Keff = 1e-5 # Convergence criterion for the Keff value
underRelaxationFactor = 0.5

########## Fields of the TH problem ##########
TeffFuel = []
Twater = []
rho = []
Qfiss = []



##### Begin Calculation scheme for coupled neutronics and thermalhydraulics solution to the BWR pincell problem.

# 1.) Guess the axial power shape : used to initialize the TH solution
qFiss = guessAxialPowerShape(PFiss, Iz1, height, Fuel_volume)
print(f"qFiss = {qFiss}")


# 2.) TH solution for initial guess of power shape : sine shape (if devide by Iz : doesnt work and the initial TH solution isn't in IAPWS domain)
# This can probably be fixed by giving a different initial guess for the total fission power.
## 2.1) Initial thermal hydraulic resolution
THComponent = THM_prototype("BWR Pincell equivalent canal", canalType, waterRadius, fuelRadius, gapRadius, cladRadius, 
                            height, tInlet, pOutlet, massFlowRate, qFiss, kFuel, Hgap, kClad, Iz1, If, I1, zPlotting, 
                            solveConduction, dt = 0, t_tot = 0, frfaccorel = frfaccorel, P2Pcorel = P2Pcorel, voidFractionCorrel = voidFractionCorrel, 
                            numericalMethod = numericalMethod)

TeffTEMP, TwaterTEMP, rhoTEMP = THComponent.get_nuclear_parameters() # rename this function to be clear about what it does 
#                                                           ---> Fuel and coolant temperature + coolant density aren't nuclear parameters per-se
#
TeffFuel.append(TeffTEMP)
Twater.append(TwaterTEMP)
rho.append(rhoTEMP)

# Overwriting the initial Teff field since THM_conduction is broken ?
TeffTEMP = [750.0 for i in range(Iz1)]

print(f"After initialization : TeffFuel = {TeffFuel}, Twater = {Twater}, rho = {rho}")

# 2.2) Create LCM object to store the TH data to be used in the neutronics solution
THData = lcm.new('LCM','THData')
THData['TFuelList']    = np.array(TeffTEMP, dtype='f')
THData['TCoolList'] = np.array(TwaterTEMP, dtype='f')
THData['DCoolList'] = np.array(rhoTEMP/1000, dtype='f')
THData.close() # close without erasing

## 3.) Initializing Neutronics solution
# 3.1) construct the Lifo stack for IniDONJON
ipLifo1=lifo.new()
ipLifo1.pushEmpty("Fmap", "LCM") # Fuel Map
ipLifo1.pushEmpty("Matex", "LCM") # Material Indexation
ipLifo1.pushEmpty("Cpo", "LCM") # Compo
ipLifo1.pushEmpty("Track", "LCM") # Tracking data for FEM
ipLifo1.push(THData) # Thermal Hydraulic data for initialization

# 3.2) call IniDONJON Cle-2000 procedure
IniDONJON = cle2000.new('IniDONJON',ipLifo1,1)
IniDONJON.exec()
print("IniDONJON execution completed")

# recover the output LCM objects
Fmap = ipLifo1.node("Fmap")
Matex = ipLifo1.node("Matex")
Cpo = ipLifo1.node("Cpo")
Track = ipLifo1.node("Track")
RecoverTHData = ipLifo1.node("THData") # Recover the TH data after the initialization is this needed ? 
#       --> THData is just used to initialize the Fuel Map object but not modified so can use initial THData object

# State vector and number of parameters can be used to modify the Fmap object directly ?
stateVector = Fmap["STATE-VECTOR"]
mylength = stateVector[0]*stateVector[1]
npar = stateVector[7]

print("Recovered stateVector: ", stateVector)
print("Number of parameters: ", npar)

# empty the Lifo stack for IniDONJON
while ipLifo1.getMax() > 0:
  ipLifo1.pop();

# 4.) Iterative scheme for coupled neutronics and thermalhydraulics solution
# Attempting to run multiPhysics scheme.
# Everytime, solve neutron transport, update axial power shape, solve TH, update TH data, and repeat.
powi = 0.1722 # Reference at 0.1722 MW from AT10_24UOX test ?

# check powi == PFiss
print(f"powi = {powi} MW and PFiss = {PFiss} W")

# Create Lifo stack for Neutronics solution
ipLifo2 = lifo.new()
Neutronics = cle2000.new('Neutronics',ipLifo2,1)
Keffs = [] # list to store the Keff values at each iteration and check convergence on it

## MultiPhysics resolution

iter = 0


conv = False
while not conv:
    iter+=1
    # 4.1) Neutronics solution for TH parameters obtained with initial TH solution

    ################## Neutronics part ##################
    # fill the Lifo stack for Neutronics solution
    print(f"iter = {iter}")
    ipLifo2.push(Fmap);
    ipLifo2.push(Matex);
    if iter == 1: # At the first iteration, create empty LCM objects to host the Flux and Power fields
        print("in iter = 1")
        Flux = ipLifo2.pushEmpty("Flux", "LCM")
        Power = ipLifo2.pushEmpty("Power", "LCM")
    else:
        ipLifo2.push(Flux)
        ipLifo2.push(RecoveredPower)
        print("Flux and Power at iteration > 1")

    ipLifo2.push(Cpo) # Push COMPO object
    ipLifo2.push(Track) # Push Tracking object for FEM solution : obtained from IniDONJON
    if iter == 1:
        ipLifo2.push(THData) # Push Thermal Hydraulic data 
    else:
        ipLifo2.push(RecoveredTHData)
    ipLifo2.push(iter)
    ipLifo2.push(powi) 

    # 4.2) Call Neutronics component :
    print("call Neutronics procedure at iter=", iter)
    Neutronics.exec()
    print("Neutronics.c2m execution completed")
    Flux = ipLifo2.node("Flux") # Recover the Flux field
    RecoveredPower = ipLifo2.node("Power") # Recover the Power field --> used to update the axial power shape
    Keff = Flux["K-EFFECTIVE"][0] # Recover the Keff value
    Keffs.append(Keff)
    print(f"At iter {iter} : Keff = {Keff}")
    PowerDistribution = RecoveredPower["POWER-DISTR"]
    print(f"Power distribution : {PowerDistribution} kW")
    #print(f"Uniform TH data used for initialization : {THData["THData"]}")
    qFiss = PowerDistribution*1000 # Updating the axial power shape, converting to W from kW
    
    #Power.val()
    #Power.close() # close without erasing
    

    ############# Thermalhydraulic part ##############
    THMComponent = THM_prototype("BWR Pincell equivalent canal", canalType, waterRadius, fuelRadius, gapRadius, cladRadius, 
                            height, tInlet, pOutlet, massFlowRate, qFiss, kFuel, Hgap, kClad, Iz1, If, I1, zPlotting, 
                            solveConduction, dt = 0, t_tot = 0, frfaccorel = frfaccorel, P2Pcorel = P2Pcorel, voidFractionCorrel = voidFractionCorrel, 
                            numericalMethod = numericalMethod)    ##### qFiss to be updated

    TeffTEMP, TwaterTEMP, rhoTEMP = THComponent.get_nuclear_parameters()
    TeffTEMP = np.array([750.0 for i in range(Iz1)])
    print(f"THM resolution at iter={iter} : TeffFuel = {TeffTEMP}, Twater = {TwaterTEMP}, rho = {rhoTEMP}")
    TeffFuel.append(TeffTEMP)
    Twater.append(TwaterTEMP)
    rho.append(rhoTEMP)
    print(f"THM resolution at iter={iter} : TeffFuel = {TeffFuel}, Twater = {Twater}, rho = {rho}")

    RecoveredTHData = lcm.new('LCM','THData')
    RecoveredTHData['TFuelList']    = np.array(TeffTEMP, dtype='f')
    RecoveredTHData['TCoolList'] = np.array(TwaterTEMP, dtype='f')
    RecoveredTHData['DCoolList'] = np.array(rhoTEMP/1000, dtype='f')
    RecoveredTHData.close() # close without erasing

    ############## Under relaxation of TH fields #################
    TeffFuel[-1] = underRelaxation(TeffFuel[-1], TeffFuel[-2], underRelaxationFactor)
    Twater[-1] = underRelaxation(Twater[-1], Twater[-2], underRelaxationFactor)
    rho[-1] = underRelaxation(rho[-1], rho[-2], underRelaxationFactor)

    ############## Convergence test #################
    #print(f"Convergence test at iter={iter} : TeffFuel = {TeffFuel}, Twater = {Twater}, rho = {rho}")
    #if convergence(TeffFuel[-1], TeffFuel[-2], tol) and convergence(Twater[-1], Twater[-2], tol) and convergence(rho[-1], rho[-2], tol):
    #    print("Convergence reached after ", iter, " iterations")
    #    break

    # empty the ipLifo2 Lifo stack
    while ipLifo2.getMax() > 0:
        print("Clearing ipLifo2 stack") 
        ipLifo2.pop();
    
    if iter>1 and convergence(Keffs[-1], Keffs[-2], tol_Keff):
        conv = True

    if iter == nIter-1:
        print("Convergence not reached after ", iter, " iterations")



    """

    """






    
