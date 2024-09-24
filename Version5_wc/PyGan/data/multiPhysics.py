## PyGan script to couple DONJON procedures with THM solution
# Date : 20/09/2024
# Author : Clément HUET, Raphaël GUASCH
# Purpose : test and validate neutronics and thermalhydraulics coupling on a single BWR pincell
from THM_main import Version5_THM_prototype as THM_prototype
from THM_main import plotting
from iapws import IAPWS97
import matplotlib.pyplot as plt
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
    print(f"Field = {Field}")
    print(f"OldField = {OldField}")
    print(f"Field - OldField = {Field - OldField}")
    conv=False
    for i in range(len(Field)):
        if np.abs(Field[i] - OldField[i]) > tol:
            conv = False
            break
        else:
            conv = True
    return conv

def guessAxialPowerShape(Ptot, Iz, height, radius):
    """
    Ptot : float : total power released (W)
    Iz : int : number of control volumes in the axial direction
    height : float : height of the fuel rod (m)
    radius : float : radius of the fuel rod (m)
    return : np.array : axial power shape with a sine shape units (W/m3)
                        --> corresponds to the power density in each control volume 
                        !! Issue with IAPWS tables when dividing by Iz
    """
    volume = np.pi * radius**2 * height
    
    # Heights of each control volume (equally spaced along the tube height)
    heights = np.linspace(0, height, Iz + 1)
    
    # Define the power profile as a sine function of height
    power_profile = lambda h: np.sin(np.pi * h / height)
    
    # Compute the volumic power for each control volume
    volumic_powers = []
    total_integral = 0
    
    for i in range(Iz):
        # Midpoint of the control volume
        h_mid = 0.5 * (heights[i] + heights[i + 1])
        print(f"Height = {h_mid}")
        
        # Power density at this control volume
        power_density = power_profile(h_mid)
        print(f"Power density = {power_density}")
        
        # Volume of this control volume
        dz = (heights[i + 1] - heights[i])
        
        # Store the volumic power (W/m^3)
        volumic_powers.append(power_density)
        
        # Update total integral for normalization
        total_integral += power_density * dz

    print(f"Total_integral = {total_integral}")
    
    # Normalize the volumetric powers so the total power matches Ptot
    volumic_powers = np.array(volumic_powers) * Ptot /(total_integral*np.pi*radius**2)/Iz
    print(f"Volumic powers = {volumic_powers}")
    total_power = np.sum(volumic_powers) * volume
    print(f"Total power = {total_power}")
    
    return volumic_powers   

def compute_difference_fields(field):
    """
    Compute the difference between the last two fields of the list field
    field : list : list of the fields to compare
    """
    print(field)
    diff = np.abs(field[-1] - field[-2])
    print(f"Difference between the last two fields = {diff}")
    return diff

def compute_residuals(field):
    """
    Compute the residuals of the field
    field : list : list of the fields to compute the residuals
    """
    residuals = np.abs(field[-1] - field[-2])/field[-2]
    print(f"Residuals of the field = {residuals}")
    return residuals

def quickPlot(x, y, title, xlabel, ylabel, saveName, path, SAVE_DIR):
    fig,ax = plt.subplots()
    if len(y) == len(x):
        ax.scatter(x, y)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        os.chdir(SAVE_DIR)
        fig.savefig(saveName)
        os.chdir(path)
    else:
        for i in range(len(y)):
            data = y[i]
            ax.plot(x, data, label=f"iteration {i}")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend()
        ax.set_title(title)
        os.chdir(SAVE_DIR)
        fig.savefig(saveName)
        os.chdir(path)
    return
######## End functions declaration ##########


########## User input ##########

solveConduction = True
zPlotting = [] #If empty, no plotting of the axial distribution of the fields, otherwise, list of the axial positions where the fields are plotted

########## Thermal hydraulics parameters ##########
## Geometric parameters
canalType = "square"
waterRadius = 1.295e-2 # m ATRIUM10 pincell pitch
fuelRadius = 0.4435e-2 # m : fuel rod radius
gapRadius = 0.4520e-2 # m : expansion gap radius : "void" between fuel and clad - equivalent to inner clad radius
cladRadius = 0.5140e-2 # m : clad external radius
height = 3.8 # m : height : active core height in BWRX-300 SMR


## Fluid parameters

# T_inlet, T_outlet = 270, 287 Celcius
tInlet = 270 + 273.15 # K
#Nominal operating pressure = 7.2 MPa (abs)
pOutlet =  7.2e6 # Pa 
# Nominal coolant flow rate = 1530 kg/s
massFlowRate = 1530  / (200*91)  # kg/s

## Additional parameters needed for the calculation
solveConduction = True
volumic_mass_U = 19000 # kg/m3
Fuel_volume = np.pi*fuelRadius**2*height # m3
Fuel_mass = Fuel_volume*volumic_mass_U # kg
print(f"Fuel mass = {Fuel_mass} kg")

## Meshing parameters:
If = 8
I1 = 3
Iz1 = 20 # number of control volumes in the axial direction for now only 20 control volumes are supported for DONJON solution

## Thermalhydraulics correlation
voidFractionCorrel = "HEM1"
frfaccorel = "base"
P2Pcorel = "base"
numericalMethod = "FVM"

######## Creation of results directory ##########
path=os.getcwd()
a=os.path.exists(f"multiPhysics_PyGan_24UOX_cell")
if a==False:
	os.mkdir(f"multiPhysics_PyGan_24UOX_cell")
print(path)

SAVE_DIR = f"multiPhysics_PyGan_24UOX_cell/{numericalMethod}/{voidFractionCorrel}/"

a=os.path.exists(SAVE_DIR)
if a==False:
	os.makedirs(SAVE_DIR)
########## End creation of results directory ##########
    
########### Print options ###########
FULL_PRINT = True # Print all the fields at each iteration, create figures in tmp/rundir/ ?

############ Nuclear Parameters ###########
## Fission parameters
# specific power = 38.6 W/g
specificPower = 38.60 # W/g, multiplied by 5 to have a more realistic value and create boiling
PFiss = specificPower*Fuel_mass*1000 # W
print("PFiss = ", PFiss)

compo_name = "_COMPO_24UOX"

#qFiss = PFiss/Fuel_volume # W/m3

## Material parameters
kFuel = 4.18 # W/m.K, TECHNICAL REPORTS SERIES No. 59 : Thermal Conductivity of Uranium Dioxide, IAEA, VIENNA, 1966
Hgap = 10000 
kClad = 21.5 # W/m.K, Thermal Conductivity of Zircaloy-2 (as used in BWRX-300) according to https://www.matweb.com/search/datasheet.aspx?MatGUID=eb1dad5ce1ad4a1f9e92f86d5b44740d
# k_Zircaloy-4 = 21.6 W/m.K too so check for ATRIUM-10 clad material but should have the same thermal conductivity
########## Algorithm parameters ###########
nIter = 1000
tol_TH = 1e-4 # Convergence criterion for the TH solution
tol_POW = 1e-3 # Convergence criterion for the power axial distribution
tol_Keff = 1e-5 # Convergence criterion for the Keff value
TH_underRelaxationFactor = 0.5
Pow_underRelaxationFactor = 0.5 # Under relaxation factor for the power axial distribution to be tested

########## Fields of the TH problem ##########
TeffFuel = []
Twater = []
rho = []
Qfiss = []

Residuals_TeffFuel = []
Residuals_Twater = []
Residual_rho = []
Residual_Qfiss = []

####### Power shape to store the axial power shape at each iteration
Keffs = [] # list to store the Keff values at each iteration and check convergence on it
Power_Distrib = []
Residuals_Power_Distrib = []

##### Begin Calculation scheme for coupled neutronics and thermalhydraulics solution to the BWR pincell problem.

# 1.) Guess the axial power shape : used to initialize the TH solution
qFiss_init = guessAxialPowerShape(PFiss, Iz1, height, fuelRadius)
print("$$$ - multiPhysics.py : BEGIN INITIALIZATION- $$$")

Qfiss.append(qFiss_init)


# 2.) TH solution for initial guess of power shape : sine shape (if devide by Iz : doesnt work and the initial TH solution isn't in IAPWS domain)
# This can probably be fixed by giving a different initial guess for the total fission power.
## 2.1) Initial thermal hydraulic resolution
THComponent = THM_prototype("BWR Pincell equivalent canal", canalType, waterRadius, fuelRadius, gapRadius, cladRadius, 
                            height, tInlet, pOutlet, massFlowRate, qFiss_init, kFuel, Hgap, kClad, Iz1, If, I1, zPlotting, 
                            solveConduction, dt = 0, t_tot = 0, frfaccorel = frfaccorel, P2Pcorel = P2Pcorel, voidFractionCorrel = voidFractionCorrel, 
                            numericalMethod = numericalMethod)

TeffTEMP, TwaterTEMP, rhoTEMP = THComponent.get_TH_parameters() # renamed this function to be clear about what it does 
#                                                           ---> Fuel and coolant temperature + coolant density aren't nuclear parameters per-se
#
# Overwriting the initial Teff field : Conduction works again but current COMPO doesnt have low enough TFuel
#TeffTEMP = np.array([750.0 for i in range(Iz1)]) --> This should be fixed!
TeffFuel.append(TeffTEMP)
Twater.append(TwaterTEMP)
rho.append(rhoTEMP)



print(f"$$ - FIRST THM resolution iter = 0 : TeffFuel = {TeffTEMP}, Twater = {TwaterTEMP}, rho = {rhoTEMP}")
print(f"After initialization : TeffFuel = {TeffFuel}, Twater = {Twater}, rho = {rho}")

# 2.2) Create LCM object to store the TH data to be used in the neutronics solution
print("Creating initial THData object")
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
ipLifo1.push(compo_name) # Compo name

# 3.2) call IniDONJON Cle-2000 procedure
IniDONJON = cle2000.new('IniDONJON',ipLifo1,1)
IniDONJON.exec()
print("IniDONJON execution completed")

# recover the output LCM objects
Fmap = ipLifo1.node("Fmap")
Matex = ipLifo1.node("Matex")
Cpo = ipLifo1.node("Cpo")
Track = ipLifo1.node("Track")
InitTHData = ipLifo1.node("THData") # Recover the TH data after the initialization is this needed ? 
#       --> THData is just used to initialize the Fuel Map object but not modified so can use initial THData object
print(f"Initial THData object TFuel = : {InitTHData['TFuelList']}, TCool = {InitTHData['TCoolList']}, DCool = {InitTHData['DCoolList']}")
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
powi = PFiss/1e6 # Reference at 0.1722 MW from AT10_24UOX test ?

# check powi == PFiss
print(f"powi = {powi} MW and PFiss = {PFiss} W")

# Create Lifo stack for Neutronics solution
ipLifo2 = lifo.new()
Neutronics = cle2000.new('Neutronics',ipLifo2,1)


## Multi-Physics resolution
iter = 0

conv = False
while not conv:
    iter+=1
    # 4.1) Neutronics solution for TH parameters obtained with initial TH solution

    ################## Neutronics part ##################
    # fill the Lifo stack for Neutronics solution
    print(f"$$ - BEGIN iter = {iter}")
    ipLifo2.push(Fmap);
    ipLifo2.push(Matex);
    if iter == 1: # At the first iteration, create empty LCM objects to host the Flux and Power fields
        print("in iter = 1")
        Flux = ipLifo2.pushEmpty("Flux", "LCM")
    else:
        ipLifo2.push(Flux)
        print("Flux and Power at iteration > 1")
    Power = ipLifo2.pushEmpty("Power", "LCM")
    ipLifo2.push(Cpo) # Push COMPO object
    ipLifo2.push(Track) # Push Tracking object for FEM solution : obtained from IniDONJON
    if iter%2 == 1:
        if iter == 1:
            ipLifo2.push(InitTHData)
            print("Pushing THData object at iter = 1")
        else:
            ipLifo2.push(UpdatedTHData2)
            print("Pushing THData object 2")
    elif iter%2 == 0:
        print("Pushing THData object 1")
        ipLifo2.push(UpdatedTHData1)
    #ipLifo2.push(UpdatedTHData) # Push the TH data obtained from the TH solution
    ipLifo2.push(iter)
    ipLifo2.push(powi) 

    # 4.2) Call Neutronics component :
    print("call Neutronics procedure at iter=", iter)
    Neutronics.exec()
    print("Neutronics.c2m execution completed")
    Flux = ipLifo2.node("Flux") # Recover the Flux field
    RecoveredPower = ipLifo2.node("Power") # Recover the Power field (LCM object)
    # 4.2.1) Recover the Keff value
    Keff = Flux["K-EFFECTIVE"][0] 
    Keffs.append(Keff)
    print(f"At iter {iter} : Keff = {Keff}")
    # 4.2.2) Recover the Power axial distribution obtained from neutronics calculation --> used to update the axial power shape
    PowerDistribution = RecoveredPower["POWER-DISTR"] 
    Power_Distrib.append(PowerDistribution) # Store the axial power shape at each iteration
    print(f"Power distribution : {PowerDistribution} kW")
    if iter > 1:
        Residuals_Power_Distrib.append(compute_difference_fields(Power_Distrib))
        if FULL_PRINT:
            quickPlot(range(Iz1), Power_Distrib, "Power distribution convergence", "axial position (ctrl vol)", "Power (kW)", "Power_distribution_convergence.png", path, SAVE_DIR)
    # 4.3) Update the axial power shape to be used in the TH solution
    qFiss = PowerDistribution*1000 # Updating the axial power shape, converting to W from kW
    Qfiss.append(qFiss)

    diff_source = compute_difference_fields(Qfiss)
    print(f"updating source: difference = {diff_source}")

    print(f"Updated axial power shape : {qFiss}")
    print(f"Initial axial power shape : {qFiss_init}")
    # 4.4) Erase the THData object to store the updated TH data
    if iter == 1:
        InitTHData.erase()
    if iter > 1 and iter % 2 == 0:
        print("Erasing THData object 1")
        UpdatedTHData1.erase() # erase the THData object
    elif iter > 1 and iter % 2 == 1:
        print("Erasing THData object 2")
        UpdatedTHData2.erase()

    # 5.) TH procedure for updated power shape
    ############# Thermalhydraulic part ##############
    # 5.1) TH resolution with updated power shape :
    THMComponent = THM_prototype("BWR Pincell equivalent canal", canalType, waterRadius, fuelRadius, gapRadius, cladRadius, 
                            height, tInlet, pOutlet, massFlowRate, qFiss, kFuel, Hgap, kClad, Iz1, If, I1, zPlotting, 
                            solveConduction, dt = 0, t_tot = 0, frfaccorel = frfaccorel, P2Pcorel = P2Pcorel, voidFractionCorrel = voidFractionCorrel, 
                            numericalMethod = numericalMethod)    ##### qFiss to be updated

    # 5.2) recover the new TH data
    TeffTEMP, TwaterTEMP, rhoTEMP = THComponent.get_TH_parameters()
    #TeffTEMP = np.array([750.0 for i in range(Iz1)])
    #TeffTEMP = [750.0 for i in range(Iz1)]
    print(f"THM resolution at iter={iter} : TeffFuel = {TeffTEMP}, Twater = {TwaterTEMP}, rho = {rhoTEMP}")
    TeffFuel.append(TeffTEMP)
    Residuals_TeffFuel.append(compute_residuals(TeffFuel))
    Twater.append(TwaterTEMP)
    Residuals_Twater.append(compute_residuals(Twater))
    rho.append(rhoTEMP)
    Residual_rho.append(compute_residuals(rho))
    print(f"THM resolution at iter={iter} : TeffFuel = {TeffFuel}, Twater = {Twater}, rho = {rho}")

    if FULL_PRINT:
        quickPlot(range(Iz1), TeffFuel, "Fuel temperature convergence", "axial position (ctrl vol)", "TFuel (K)", "TFuel_convergence.png", path, SAVE_DIR)
        quickPlot(range(Iz1), Twater, "Coolant temperature convergence", "axial position (ctrl vol)", "TCool (K)", "TCool_convergence.png", path, SAVE_DIR)
        quickPlot(range(Iz1), rho, "Coolant density convergence", "axial position (ctrl vol)", "DCool (kg/m3)", "DCool_convergence.png", path, SAVE_DIR)
        quickPlot(range(Iz1), Qfiss, "Axial power shape", "axial position (ctrl vol)", "Power (W)", "Power_shape.png", path, SAVE_DIR)

    
    diff = compute_difference_fields(rho)
    print(f"updating rho: difference = {diff}")
    # 5.3) Update the TH data to be used in the next neutronics solution
    
    if iter%2 == 1:
        UpdatedTHData1 = lcm.new('LCM','THData') # Create new LCM object to store the updated TH data
        print("Creating updated THData object 1")
        print(UpdatedTHData1)
        UpdatedTHData1['TFuelList']    = np.array(TeffTEMP, dtype='f')
        UpdatedTHData1['TCoolList'] = np.array(TwaterTEMP, dtype='f')
        UpdatedTHData1['DCoolList'] = np.array(rhoTEMP/1000, dtype='f') # Convert density from kg/m3 to g/cm3
        #UpdatedTHData.val() # validate the LCM object
        UpdatedTHData1.close() # close without erasing
    elif iter%2 == 0:
        UpdatedTHData2 = lcm.new('LCM','THData')
        print("Creating updated THData object 2")
        print(UpdatedTHData2)
        UpdatedTHData2['TFuelList']    = np.array(TeffTEMP, dtype='f')
        UpdatedTHData2['TCoolList'] = np.array(TwaterTEMP, dtype='f')
        UpdatedTHData2['DCoolList'] = np.array(rhoTEMP/1000, dtype='f') # Convert density from kg/m3 to g/cm3
        UpdatedTHData2.close() # close without erasing

    print("Updated THData object created")
    #THData_List.append(UpdatedTHData) # Store the TH data in a list to be used in the next neutronics resolution
    # 5.4) Under relaxation of TH fields for the next iteration : check if it is needed ?
    TeffFuel[-1] = underRelaxation(TeffFuel[-1], TeffFuel[-2], TH_underRelaxationFactor)
    Twater[-1] = underRelaxation(Twater[-1], Twater[-2], TH_underRelaxationFactor)
    rho[-1] = underRelaxation(rho[-1], rho[-2], TH_underRelaxationFactor)

    # Under relaxation of the Power distribution?
    if iter > 1:

        Power_Distrib[-1] = underRelaxation(Power_Distrib[-1], Power_Distrib[-2], Pow_underRelaxationFactor)
        if FULL_PRINT:
            quickPlot(range(Iz1), Power_Distrib, "Power distribution convergence", "axial position (ctrl vol)", "Power (kW)", "Power_distribution_convergence_relaxed.png", path, SAVE_DIR)


    # 6.) Empty the ipLifo2 Lifo stack to prepare for the next iteration
    while ipLifo2.getMax() > 0:
        print("Clearing ipLifo2 stack at iter = ", iter) 
        ipLifo2.pop();
    
    print(f"$$ - END iter = {iter}")
    # 7.) Check for convergence on Keff value : add more converegence criteria for other fields ?
    if iter>1 and np.abs(Keffs[-1]-Keffs[-2])< tol_Keff:
        print("Convergence on Keff reached after ", iter, " iterations")
        if iter>1 and convergence(rho[-1], rho[-2], tol_TH) and convergence(TeffFuel[-1], TeffFuel[-2], tol_TH) and convergence(Twater[-1], Twater[-2], tol_TH):
            if iter>1 and convergence(Power_Distrib[-1], Power_Distrib[-2], tol_POW):
                conv = True
                print("Convergence on TH fields and Power reached after ", iter, " iterations")
    if iter == nIter-1:
        print("Convergence not reached after ", iter, " iterations")


    print("$$$ - multiPhysics.py : END ITERATIONS - $$$")


# 8.) Plot the results

# Creation of results directory
path=os.getcwd()
a=os.path.exists(f"multiPhysics_PyGan_24UOX_cell")
if a==False:
	os.mkdir(f"multiPhysics_PyGan_24UOX_cell")
print(path)

SAVE_DIR = f"multiPhysics_PyGan_24UOX_cell/{numericalMethod}/{voidFractionCorrel}/"

a=os.path.exists(SAVE_DIR)
if a==False:
	os.makedirs(SAVE_DIR)

print(SAVE_DIR)
print(f"Keffs = {Keffs}")
# Plot the results
quickPlot(range(len(Keffs)), Keffs, "Keff convergence", "iteration", "Keff", "Keff_convergence.png", path, SAVE_DIR)

quickPlot(range(Iz1), TeffFuel, "Fuel temperature convergence", "axial position (ctrl vol)", "TFuel (K)", "TFuel_convergence.png", path, SAVE_DIR)
quickPlot(range(Iz1), Twater, "Coolant temperature convergence", "axial position (ctrl vol)", "TCool (K)", "TCool_convergence.png", path, SAVE_DIR)
quickPlot(range(Iz1), rho, "Coolant density convergence", "axial position (ctrl vol)", "DCool (kg/m3)", "DCool_convergence.png", path, SAVE_DIR)
quickPlot(range(Iz1), Power_Distrib, "Power distribution convergence", "axial position (ctrl vol)", "Power (kW)", "Power_distribution_convergence.png", path, SAVE_DIR)
quickPlot(range(Iz1), Qfiss, "Axial power shape", "axial position (ctrl vol)", "Power (W)", "Power_shape.png", path, SAVE_DIR)
    
quickPlot(range(len(Residuals_TeffFuel)), Residuals_TeffFuel, "Fuel temperature residuals", "iteration", "Residuals", "TFuel_residuals.png", path, SAVE_DIR)
quickPlot(range(len(Residuals_Twater)), Residuals_Twater, "Coolant temperature residuals", "iteration", "Residuals", "TCool_residuals.png", path, SAVE_DIR)
quickPlot(range(len(Residual_rho)), Residual_rho, "Coolant density residuals", "iteration", "Residuals", "DCool_residuals.png", path, SAVE_DIR)
quickPlot(range(len(Residuals_Power_Distrib)), Residuals_Power_Distrib, "Power distribution residuals", "iteration", "Residuals", "Power_distribution_residuals.png", path, SAVE_DIR)
quickPlot(range(len(Residual_Qfiss)), Residual_Qfiss, "Axial power shape residuals", "iteration", "Residuals", "Power_shape_residuals.png", path, SAVE_DIR)