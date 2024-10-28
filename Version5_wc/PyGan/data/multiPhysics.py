## PyGan script to couple DONJON procedures with THM solution
# Date : 20/09/2024
# Author : Clément HUET, Raphaël GUASCH
# Purpose : test and validate neutronics and thermalhydraulics coupling on a single BWR pincell
from THM_main import Version5_THM_prototype as THM_prototype
from THM_main import plotting
from iapws import IAPWS97
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simps
import os, shutil
import lifo
import lcm
import cle2000
from assertS import *
import time


########## Begin helper functions declaration ##########
## Fonction used for the convergence algorithm
def underRelaxation(Field, OldField, underRelaxationFactor):
    print(f"Field = {Field}")
    print(f"OldField = {OldField}")
    underRelaxed_Field = underRelaxationFactor*Field + (1-underRelaxationFactor)*OldField
    print(f"Under relaxed Field = {underRelaxed_Field}")
    return underRelaxationFactor*Field + (1-underRelaxationFactor)*OldField

def convergence_TH(Field, OldField, tol):
    conv=False
    print("Checking TH convergence")
    rms_error = np.sqrt(np.sum((Field - OldField)**2)/len(Field))
    print(f"TH rms error = {rms_error}")
    for i in range(len(Field)):
        if np.abs(Field[i] - OldField[i]) > tol:
            conv = False
            break
        else:
            conv = True
    return conv

def convergence_POW(Field, OldField, tol):
    """
    Tolerance in % of the nodal power, to allow for an easier comparison between different total powers.
    """
    conv=False
    print("Checking POW convergence")
    rms_error = np.sqrt(np.sum((Field - OldField)**2)/len(Field))
    print(f"POW rms error = {rms_error}")
    for i in range(len(Field)):
        if np.abs(Field[i] - OldField[i])/OldField[i] > tol:
            conv = False
            break
        else:
            conv = True
    return conv

def convergence_pcm(keffNew, keffOld, tol):
    print(f"keffNew = {keffNew}")
    print(f"keffOld = {keffOld}")
    print(f"keffNew - keffOld = {keffNew - keffOld}")
    conv=False
    print("Checking convergence on keff")
    print(f"Error on keff = {np.abs(keffNew - keffOld)*1e5} pcm")
    print(f"tol = {tol} pcm")   
    if np.abs(keffNew - keffOld)*1e5 > tol:
        conv = False
    else:
        conv = True
    return conv

def guessAxialPowerShape(Ptot, Iz, height, radius): # Old version of the function, introduced errrors when scaling up to more control volumes but seemed to work for 20 control volumes
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

def guess_power_density_sine(Ptot, h, r, num_control_volumes): # New version, supposed to be more accurate but introduced alot of instability in the TH solution
    """
    Compute the power density in a cylinder with a sine-shaped profile along the axial direction.

    Parameters:
    Ptot (float): Total power released in the cylinder [W].
    h (float): Total height of the cylinder [m].
    r (float): Radius of the cylinder [m].
    num_control_volumes (int): Number of control volumes along the z-axis.

    Returns:
    z_values (numpy array): Array of z positions (center of each control volume) along the height of the cylinder.
    power_density (numpy array): Array of power densities corresponding to each control volume.
    q0 (float): The peak power density constant.
    """
    
    # Cross-sectional area of the cylinder
    A_cross_section = np.pi * r**2

    # Function for the power density profile along the axial direction
    def q(z, q0):
        return q0 * np.sin(np.pi * z / h)

    # Calculate q0 to ensure total power matches Ptot
    q0 = Ptot / (A_cross_section * (2 * h / np.pi))  # derived from the integral

    # Compute the boundaries and midpoints of each control volume
    z_boundaries = np.linspace(0, h, num_control_volumes + 1)
    z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes

    # Compute power density at each midpoint
    power_density = q(z_values, q0)

    return z_values, power_density

def guess_power_density_cosine(Ptot, h, r, num_control_volumes):
    """
    Compute the power density in a cylinder with a cosine-shaped profile along the axial direction.

    Parameters:
    Ptot (float): Total power released in the cylinder [W].
    h (float): Total height of the cylinder [m].
    r (float): Radius of the cylinder [m].
    num_control_volumes (int): Number of control volumes along the z-axis.

    Returns:
    z_values (numpy array): Array of z positions (center of each control volume) along the height of the cylinder.
    power_density (numpy array): Array of power densities corresponding to each control volume.
    """
    
    # Cross-sectional area of the cylinder
    A_cross_section = np.pi * r**2

    # Function for the power density profile along the axial direction (cosine shape)
    def q(z, q0):
        return q0 * np.cos(np.pi * z / (2 * h))

    # Calculate q0 to ensure total power matches Ptot
    #q0 = Ptot / (A_cross_section * (h / 2))  # derived from the integral of cosine
    q0 = Ptot*np.pi / (A_cross_section * (h * 2))  # derived from the integral of cosine

    # Compute the boundaries and midpoints of each control volume
    z_boundaries = np.linspace(0, h, num_control_volumes + 1)
    z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes

    # Compute power density at each midpoint
    power_density = q(z_values, q0)

    return z_values, power_density

def compute_power_densities(integrated_powers, r, Iz, height):
    """
    Convert a list of integrated powers at each axial node into power densities.

    Parameters:
    integrated_powers (numpy array): List of integrated powers at each node [W].
    r (float): Radius of the cylinder [m].
    Iz (int) : number of axial-subdivisions.

    Returns:
    power_densities (numpy array): Power densities at each axial node [W/m^3].
    """
    print(len(integrated_powers))
    print(f"integrated_powers = {integrated_powers} W")

    # Cross-sectional area of the cylinder (constant for all nodes)
    A_cross_section = np.pi * r**2
    z_values = np.linspace(0, height, Iz+1)
    print(f"z_values = {z_values}")
    # Calculate the heights (Δz) of each control volume
    dz = np.diff(z_values)  # Heights between adjacent z positions

    # Compute the volumes of each control volume
    volumes = A_cross_section * dz  # Volume of each cylinder segment

    # Calculate the power densities for each node
    power_densities = integrated_powers / volumes

    return power_densities

def integrate_power_density(z_values, power_density, r):
    """
    Numerically integrate the power density profile to compute total power and ensure it's normalized to Ptot.

    Parameters:
    z_values (numpy array): Array of z positions (center of each control volume).
    power_density (numpy array): Array of power densities corresponding to each z position.
    r (float): Radius of the cylinder [m].

    Returns:
    integrated_power (float): Numerically integrated total power from the power density profile.
    """
    
    # Cross-sectional area of the cylinder
    A_cross_section = np.pi * r**2

    # Integrate using the trapezoidal rule to get the total power
    #integrated_power = np.trapz(power_density * A_cross_section, z_values)
    
    # Integrate using Simpson's rule to get the total power
    integrated_power = simps(power_density * A_cross_section, z_values)

    return integrated_power

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
    field : list : list of the fields to compute the residuals in %
    """
    residuals = (field[-1] - field[-2])*100/field[-2]
    print(f"Residuals of the field = {residuals} %")
    return residuals

def quickPlot(x, y, title, xlabel, ylabel, saveName, path, SAVE_DIR):
    fig,ax = plt.subplots()
    if len(y) == len(x) and isinstance(y[0], np.float32):
        #print("Plotting a scatter plot")
        ax.scatter(x, y)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        os.chdir(SAVE_DIR)
        fig.savefig(saveName)
        os.chdir(path)
    else:
        for i in range(len(y)):
            if i%10==0:
                data = y[i]
                ax.plot(x, data, '2-',linewidth=1, label=f"iteration {i}")
        ax.plot(x, y[-1], '2-',linewidth=1, label=f"iteration {len(y)}")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend()
        ax.set_title(title)
        os.chdir(SAVE_DIR)
        fig.savefig(saveName)
        os.chdir(path)
    #fig.close()
    plt.close(fig)
    return
######## End helper functions declaration ##########


##
########## User input ##########
########## Algorithm parameters ###########
start_time = time.time()
solveConduction = True
nIter = 1000 #1000
tol_TH_Tf = 1e-3 # Convergence criterion for TFuel, in K
tol_TH_Tc = 1e-3 # Convergence criterion for TCool, in K
tol_TH_Dc = 1e-3 # Convergence criterion for DCool, in kg/m3
tol_POW = 0.01 # Convergence criterion for the power axial distribution, in % of nodal power
tol_Keff =  0.1 # Convergence criterion for the Keff value in pcm

relax_Pow = False # Under relaxation of the Power distribution for the next iteration
Pow_underRelaxationFactor = 0.5 # 0.8, 0.5, 0.2, Under relaxation factor for the power axial distribution to be tested, 0.1 used in Serpent/OpenFoam coupling

relax_TH = False # Under relaxation of the TH fields for the next iteration
TH_underRelaxationFactor = 0.9 # 0.8, 0.5, 0.2, Under relaxation factor for the TH fields to be tested

########## Mesh parameters ###########
zPlotting = [] #If empty, no plotting of the axial distribution of the fields, otherwise, list of the axial positions where the fields are plotted
## Meshing parameters:
If = 8
I1 = 3
# Sensitivity to the meshing parameters
Iz1 = 70 # number of control volumes in the axial direction, added 70 for comparison with GeN-Foam
# Iz1 = 10, 20, 40, 50, 70, 80 and 160 are supported for the DONJON solution


power_scaling_factor = 1 # 1, 2, 4, 8 # Scaling factor for the power axial distribution

########## Choice of Thermalhydraulics correlation ##########
voidFractionCorrel = 'EPRIvoidModel' # 'modBestion', 'HEM1', 'GEramp', 'EPRIvoidModel'
frfaccorel = "Churchill" # 'base', 'blasius', 'Churchill', 'Churchill_notOK' ?
P2Pcorel = "HEM1" # 'base', 'HEM1', 'HEM2', 'MNmodel'
numericalMethod = "BiCG" # "FVM": Solves the system using matrix inversion with preconditioning.
                        # "GaussSiedel" : Applies the Gauss-Seidel iterative solver.
                        # "BiCG" : Uses the BiConjugate Gradient method for solving non-symmetric or indefinite matrices.
                        # "BiCGStab" : Applies the BiCGStab (BiConjugate Gradient Stabilized) method to ensure faster and more stable convergence.

########## Thermal hydraulics parameters ##########
## Geometric parameters
canalType = "square" # "square", "cylindrical"
pitch = 1.295e-2 # m : ATRIUM10 pincell pitch
fuelRadius = 0.4435e-2 # m : fuel rod radius
gapRadius = 0.4520e-2 # m : expansion gap radius : "void" between fuel and clad - equivalent to inner clad radius
cladRadius = 0.5140e-2 # m : clad external radius
height = 1.555 # m : height : 3.8 m : active core height in BWRX-300 SMR, 1.555 m : for GeNFoam comparison.


## Fluid parameters

# T_inlet, T_outlet = 270, 287 Celcius
#tInlet = 270 + 273.15 # K, for BWRX-300 SMR core, try lowering the inlet temperature to set boiling point back and reduce the void fraction increase in the first few cm
tInlet = 270 + 273.15 # K, for BWRX-300 SMR core
#Nominal operating pressure = 7.2 MPa (abs)
pOutlet =  7.2e6 # Pa 
# Nominal coolant flow rate = 1530 kg/s
massFlowRate = 1530  / (200*91)  # kg/s

## Material parameters
kFuel = 4.18 # W/m.K, TECHNICAL REPORTS SERIES No. 59 : Thermal Conductivity of Uranium Dioxide, IAEA, VIENNA, 1966
Hgap = 10000 
#Hgap = 9000
kClad = 21.5 # W/m.K, Thermal Conductivity of Zircaloy-2 (as used in BWRX-300) according to https://www.matweb.com/search/datasheet.aspx?MatGUID=eb1dad5ce1ad4a1f9e92f86d5b44740d
# k_Zircaloy-4 = 21.6 W/m.K too so check for ATRIUM-10 clad material but should have the same thermal conductivity


    
########### Print options ###########
FULL_PRINT = True # Print all the fields at each iteration, create figures in tmp/rundir/ ?

############ Nuclear Parameters ###########
# Number of fuel rods and assemblies for a small modular Boiling Water Reactor core

n_rods = 91 # This value is for the ATRIUM-10 fuel assembly (10*10-9), in a GNF2 fuel assembly (BWRX-300) there are 92 fuel rods per assembly
n_assmblies = 240 # This value is for the BWRX-300 SMR core, ref : "Status Report – BWRX-300 (GE Hitachi and Hitachi GE Nuclear Energy)"
full_core_power = 870e6 # W, full core thermal power of the BWRX-300 SMR core, ref : "Status Report – BWRX-300 (GE Hitachi and Hitachi GE Nuclear Energy)"

volumic_mass_U = 10970 # kg/m3 'Thermophysical Properties of MOX and UO2 Fuels Including the Effects of Irradiation' - ORNL/TM-2000/351, Popov, Carbajo, Ivanov, Yoder. November 2000.

## Fuel rod scale parameters :
Fuel_volume = np.pi*fuelRadius**2*height # m3

#Fuel_rod_volume = np.pi*cladRadius**2*height # m3
Fuel_mass = Fuel_volume*volumic_mass_U*1000 # g
print(f"Fuel mass = {Fuel_mass} g")

print(f"$$ - BEGIN Iz1 = {Iz1}, power_scaling_factor = {power_scaling_factor}")
Bundle_volume = Fuel_volume / Iz1 # m3, Bundle <=> 1 axial slice of the fuel channel
fuel_rod_power = full_core_power/(n_rods*n_assmblies) # W
#specificPower = fuel_rod_power/Fuel_mass # W/g 
print(f"Fuel rod power before scaling = {fuel_rod_power} W")
#print(f"Specific power = {specificPower} W/g")


compo_name = "_COMPO_24UOX" # Name of the COMPO object to be used in the neutronics solution


PFiss = fuel_rod_power/power_scaling_factor # W
print(f"PFiss = {PFiss} = fuel_rod_power (scaled) W")

########## Fields of the TH problem ##########
TeffFuel = []
Twater = []
rho = []
voidFraction = []
Volumic_Powers = []
Relaxed_Volumic_Powers = []


Residuals_TeffFuel = []
Residuals_Twater = []
Residuals_rho = []
Residuals_voidFraction = []
Residuals_Volumic_Powers = []
Residuals_Relaxed_Volumic_Powers = []

####### Power shape to store the axial power shape at each iteration
Keffs = [] # list to store the Keff values at each iteration and check convergence on it
Power_Distribs = []
Relaxed_Power_Distribs = []

Residuals_Power_Distribs = []
Residuals_Relaxed_Power_Distribs = []


######## Creation of results directory ##########
path=os.getcwd()
a=os.path.exists(f"multiPhysics_PyGan_24UOX_cell")
if a==False:
    os.mkdir(f"multiPhysics_PyGan_24UOX_cell")
print(path)

SAVE_DIR = f"multiPhysics_PyGan_24UOX_cell/{numericalMethod}/{voidFractionCorrel}_{frfaccorel}_{P2Pcorel}/"

a=os.path.exists(SAVE_DIR)
if a==False:
    os.makedirs(SAVE_DIR)
########## End creation of results directory ##########

##### Begin Calculation scheme for coupled neutronics and thermalhydraulics solution to the BWR pincell problem.

# 1.) Guess the axial power shape : used to initialize the TH solution
z_mesh, qFiss_init = guess_power_density_sine(PFiss, height, fuelRadius, Iz1)
#z_mesh, qFiss_init = guess_power_density_cosine(PFiss, height, fuelRadius, Iz1)
print(f"qFiss_init = {qFiss_init}")

# Check the total power
total_power = integrate_power_density(z_mesh, qFiss_init, fuelRadius)
print(f"Total power = {total_power} W")
print(f"Error on integrated (total) power = {total_power - PFiss} W <-> {100*(total_power - PFiss)/PFiss} %")
print("$$$ - multiPhysics.py : BEGIN INITIALIZATION- $$$")

Volumic_Powers.append(qFiss_init)

initial_volumic_power = np.sum(qFiss_init)
print(f"Initial volumic power = {initial_volumic_power} W")


# 2.) TH solution for initial guess of power shape : sine shape (if devide by Iz : doesnt work and the initial TH solution isn't in IAPWS domain)
# This can probably be fixed by giving a different initial guess for the total fission power.
## 2.1) Initial thermal hydraulic resolution
THComponentIni = THM_prototype("Initialization of BWR Pincell equivalent canal", canalType, pitch, fuelRadius, gapRadius, cladRadius, 
                            height, tInlet, pOutlet, massFlowRate, qFiss_init, kFuel, Hgap, kClad, Iz1, If, I1, zPlotting, 
                            solveConduction, dt = 0, t_tot = 0, frfaccorel = frfaccorel, P2Pcorel = P2Pcorel, voidFractionCorrel = voidFractionCorrel, 
                            numericalMethod = numericalMethod)

TeffIni, TwaterIni, rhoIni, voidFracIni = THComponentIni.get_TH_parameters() # renamed this function to be clear about what it does 
#                                                           ---> Fuel and coolant temperature + coolant density aren't nuclear parameters per-se
#

THM_ini_time = (time.time() - start_time)
print(f"THM_ini_time = {THM_ini_time} s = total elapsed time = {THM_ini_time} s")


TeffFuel.append(TeffIni)
Twater.append(TwaterIni)
rho.append(rhoIni)
voidFraction.append(voidFracIni)



print(f"$$ - FIRST THM resolution iter = 0 : TeffFuel = {TeffIni}, Twater = {TwaterIni}, rho = {rhoIni}")
print(f"After initialization : TeffFuel = {TeffFuel}, Twater = {Twater}, rho = {rho}")

# 2.2) Create LCM object to store the TH data to be used in the neutronics solution
print("Creating initial THData object")
THData = lcm.new('LCM','THData')
THData['TFuelList']    = np.array(TeffIni, dtype='f')
THData['TCoolList'] = np.array(TwaterIni, dtype='f')
THData['DCoolList'] = np.array(rhoIni/1000, dtype='f')
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
ipLifo1.push(int(Iz1)) # Number of axial subdivisions
ipLifo1.push(Fuel_mass) # Mass of the fuel
ipLifo1.push(height) # Height of the fuel rod
ipLifo1.push(pitch) # pitch of the fuel assembly/cell

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

current_time0 = time.time()
DONJON_ini_time = (current_time0 - THM_ini_time)
total_elapsed_time = (current_time0 - start_time)
print(f"DONJON initialization time = {DONJON_ini_time} s")
print(f"Total elapsed time = {total_elapsed_time} s")
## Multi-Physics resolution
iter = 0


conv = False
while not conv:
    current_time1 = time.time()
    convKeff = False
    convTH = False
    convPow = False
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
    if iter == 1:
        ipLifo2.push(InitTHData)
        print("Pushing THData object at iter = 1")
    else:
        ipLifo2.push(UpdatedTHData)
        print("Pushing THData object 2")
    #ipLifo2.push(UpdatedTHData) # Push the TH data obtained from the TH solution
    ipLifo2.push(iter)
    ipLifo2.push(powi)
    ipLifo2.push(int(Iz1)) 

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
    FLUX = RecoveredPower["FLUX"]
    current_time2 = time.time()
    DONJON_iter_time = (current_time2 - current_time1)
    total_elapsed_time = (current_time2 - start_time)
    print(f"DONJON_iter_time = {DONJON_iter_time} s = total elapsed time = {total_elapsed_time} s")
    print(f"Power distribution : {PowerDistribution} W")
    qFiss = compute_power_densities(PowerDistribution, fuelRadius, Iz1, height) # Updating the axial power shape, converting to W from kW, and dividing by the fuel volume to get W/m3
    Volumic_Powers.append(qFiss)
    Power_Distribs.append(PowerDistribution) # Store the axial power shape at each iteration
    # 4.3) Update the axial power shape to be used in the TH solution
    print(f"Initial axial power shape : {qFiss_init}")
    print(f"Power distribution : {PowerDistribution}")
    print(f"Updated axial power/vol shape : {qFiss}")

    # 4.3.1) Under relaxation of the Power distribution for the next iteration
    if iter == 1 and relax_Pow:
        PowerDistribution_relaxed = underRelaxation(PowerDistribution, qFiss_init*Bundle_volume, Pow_underRelaxationFactor)
        # Under relaxation of the Power distribution?
        qFiss = compute_power_densities(PowerDistribution_relaxed, fuelRadius, Iz1, height)
    elif iter > 1 and relax_Pow:
            PowerDistribution_relaxed = underRelaxation(PowerDistribution, Relaxed_Power_Distribs[-1], Pow_underRelaxationFactor)
            qFiss = compute_power_densities(PowerDistribution_relaxed, fuelRadius, Iz1, height)
    if relax_Pow:
        Relaxed_Power_Distribs.append(PowerDistribution)
        Relaxed_Volumic_Powers.append(qFiss)


    if iter == 1:
        SAVE_DATA_DIR = f"multiPhysics_PyGan_24UOX_cell/{numericalMethod}/{voidFractionCorrel}_{frfaccorel}_{P2Pcorel}/mesh{Iz1}_{power_scaling_factor}/Figures/"
        a=os.path.exists(SAVE_DATA_DIR)
        if a==False:
            os.makedirs(SAVE_DATA_DIR)
        os.chdir(SAVE_DATA_DIR)
        heights = np.linspace(0, height, Iz1 + 1)
        #mid_heights = z_mesh
        fig, ax = plt.subplots()
        ax.plot(z_mesh, qFiss, '2-',linewidth=1, label="Initial DONJON Volumic Power Distribution")
        ax.plot(z_mesh, qFiss_init, '2-', linewidth=1, label="Guess Axial Volumic Power Distribution")
        ax.set_xlabel("height (m)")
        ax.set_ylabel("Volumic Power (W)")
        ax.legend()
        ax.set_title("Volumic Power Distributions at iter = 1")
        fig.savefig(f"Power_distributions_iter{iter}.png")
        os.chdir(path)

    if iter > 1:
        Residuals_Power_Distribs.append(compute_residuals(Power_Distribs))
        Residuals_Volumic_Powers.append(compute_residuals(Volumic_Powers))
        if relax_Pow:
            Residuals_Relaxed_Power_Distribs.append(compute_residuals(Relaxed_Power_Distribs))
            Residuals_Relaxed_Volumic_Powers.append(compute_residuals(Relaxed_Volumic_Powers))
    print(f"Initial axial power shape : {qFiss_init}")
    # 4.4) Erase the THData object to store the updated TH data
    if iter == 1:
        InitTHData.erase()
    else:
        print("Erasing THData object")
        UpdatedTHData.erase() # erase the THData object

    # 5.) TH procedure for updated power shape
    ############# Thermalhydraulic part ##############

    updated_qFiss = qFiss

    current_time3 = time.time()
    time_spent_in_plotting_res_under_relax = (current_time3 - current_time2)
    print(f"Time spent in plotting residuals under relaxation = {time_spent_in_plotting_res_under_relax} s")
    total_elapsed_time = (current_time3 - start_time)
    print(f"Total elapsed time = {total_elapsed_time} s")
    # 5.1) TH resolution with updated power shape :
    print(f"At iter {iter} : THM resolution with updated power shape")
    print(f"Updated axial power shape : {updated_qFiss}")
    THMComponent = THM_prototype("BWR Pincell equivalent canal", canalType, pitch, fuelRadius, gapRadius, cladRadius, 
                            height, tInlet, pOutlet, massFlowRate, updated_qFiss, kFuel, Hgap, kClad, Iz1, If, I1, zPlotting, 
                            solveConduction, dt = 0, t_tot = 0, frfaccorel = frfaccorel, P2Pcorel = P2Pcorel, voidFractionCorrel = voidFractionCorrel, 
                            numericalMethod = numericalMethod)    ##### qFiss updated

    # 5.2) recover the new TH data
    TeffTEMP, TwaterTEMP, rhoTEMP, voidFracTEMP = THMComponent.get_TH_parameters()
    current_time4 = time.time()
    time_spent_in_THM = (current_time4 - current_time3)
    total_elapsed_time = (current_time4 - start_time)
    print(f"Time spent in THM iter {iter} = {time_spent_in_THM} s")
    print(f"Total elapsed time = {total_elapsed_time} s")
    TeffFuel.append(TeffTEMP)
    Residuals_TeffFuel.append(compute_residuals(TeffFuel))
    Twater.append(TwaterTEMP)
    Residuals_Twater.append(compute_residuals(Twater))
    rho.append(rhoTEMP)
    Residuals_rho.append(compute_residuals(rho))
    voidFraction.append(voidFracTEMP)
    Residuals_voidFraction.append(compute_residuals(voidFraction))

    print(f"THM resolution at iter={iter} : All lists TeffFuel = {TeffFuel}, Twater = {Twater}, rho = {rho}")

    if FULL_PRINT:
        print(z_mesh)
        print(TeffFuel)
        quickPlot(z_mesh, TeffFuel, "Fuel temperature convergence", "height (m)", "TFuel (K)", "TFuel_convergence.png", path, SAVE_DIR)
        quickPlot(z_mesh, Twater, "Coolant temperature convergence", "height (m)", "TCool (K)", "TCool_convergence.png", path, SAVE_DIR)
        quickPlot(z_mesh, rho, "Coolant density convergence", "height (m)", "DCool (kg/m3)", "DCool_convergence.png", path, SAVE_DIR)
        quickPlot(z_mesh, Volumic_Powers, "Axial power shape", "height (m)", "Power (W)", "Power_shape.png", path, SAVE_DIR)

    
    # 5.3) Update the TH data to be used in the next neutronics solution
    
    UpdatedTHData = lcm.new('LCM','THData') # Create new LCM object to store the updated TH data
    UpdatedTHData['TFuelList']    = np.array(TeffTEMP, dtype='f')
    UpdatedTHData['TCoolList'] = np.array(TwaterTEMP, dtype='f')
    UpdatedTHData['DCoolList'] = np.array(rhoTEMP/1000, dtype='f') # Convert density from kg/m3 to g/cm3
    UpdatedTHData.close() # close without erasing

    # 5.4) Under relaxation of TH fields for the next iteration : check if it is needed ?
    if relax_TH:
        TeffFuel[-1] = underRelaxation(TeffFuel[-1], TeffFuel[-2], TH_underRelaxationFactor)
        Twater[-1] = underRelaxation(Twater[-1], Twater[-2], TH_underRelaxationFactor)
        rho[-1] = underRelaxation(rho[-1], rho[-2], TH_underRelaxationFactor)


    # 6.) Empty the ipLifo2 Lifo stack to prepare for the next iteration
    while ipLifo2.getMax() > 0: 
        ipLifo2.pop()
    
    print(f"$$ - END iter = {iter}")
    # 7.) Check for convergence on Keff value : add more converegence criteria for other fields ?
    if iter>1 and convergence_pcm(Keffs[-1], Keffs[-2], tol_Keff):
        print("Convergence on Keff reached after ", iter, " iterations")
        convKeff = True
    if iter>1 and convergence_TH(rho[-1], rho[-2], tol_TH_Dc) and convergence_TH(TeffFuel[-1], TeffFuel[-2], tol_TH_Tf) and convergence_TH(Twater[-1], Twater[-2], tol_TH_Tc):
        print("Convergence on TH fields reached after ", iter, " iterations")
        convTH = True
    if iter>1 and convergence_POW(Power_Distribs[-1], Power_Distribs[-2], tol_POW):
        convPOW = True
        print("Convergence on Power reached after ", iter, " iterations")
    if convKeff and convTH and convPOW:
        conv = True
        total_iter = iter
        print("Convergence on Keff, TH fields and Power reached after ", iter, " iterations")
        print("Convergence criteria enforced : tol_Keff = ", tol_Keff, " tol_TH_Tf = ", tol_TH_Tf, " tol_TH_Tc = ", tol_TH_Tc, " tol_TH_Dc = ", tol_TH_Dc, " tol_POW = ", tol_POW)
    if iter == nIter-1:
        total_iter = iter
        print("Convergence not reached after ", iter, " iterations")
        break

    current_time5 = time.time()
    time_spent_in_res_and_conv = (current_time5 - current_time4)
    total_elapsed_time = (current_time5 - start_time)
    print(f"Time spent in residuals and convergence check = {time_spent_in_res_and_conv} s")
    print(f"Total elapsed time = {total_elapsed_time} s")
    total_iteration_time = (current_time5 - current_time1)
    print(f"Total iteration time = {total_iteration_time} s for iteration {iter}")

print("$$$ - multiPhysics.py : END ITERATIONS - $$$")

end_iters_time = time.time()
time_spent_in_multiPhysics = (end_iters_time - start_time)
print(f"Time spent in multiPhysics iterations = {time_spent_in_multiPhysics} s, for a total of {total_iter} iterations, convergence reached {conv} (T/F)")

# 8.) Check the results
print("$$$ - multiPhysics.py : BEGIN RESULTS - $$$")

print(f"Checking Fields lengths : Keff : {len(Keffs)}, TeffFuel : {len(TeffFuel)}, Twater : {len(Twater)}, rho : {len(rho)}, Power_Distrib : {len(Power_Distribs)}")
print(f"Checking residuals lengths : Residuals_TeffFuel : {len(Residuals_TeffFuel)}, Residuals_Twater : {len(Residuals_Twater)}, Residual_rho : {len(Residuals_rho)}, Residuals_Power_Distrib : {len(Residuals_Power_Distribs)}")


print("$$$ - multiPhysics.py : END RESULTS - $$$")


# 9.) Plot the results

if relax_Pow:
    relaxPOW_id = f"relaxedPOW_{Pow_underRelaxationFactor}"
else:
    relaxPOW_id = "non_relaxedPOW"
if relax_TH:
    relaxTH_id = f"relaxedTH_{TH_underRelaxationFactor}"
else:
    relaxTH_id = "non_relaxedTH"
# Creation of results directory
path=os.getcwd()
a=os.path.exists(f"multiPhysics_PyGan_24UOX_cell")
if a==False:
    os.mkdir(f"multiPhysics_PyGan_24UOX_cell")
print(path)

SAVE_FIG_DIR = f"multiPhysics_PyGan_24UOX_cell/{numericalMethod}/{voidFractionCorrel}_{frfaccorel}_{P2Pcorel}/mesh{Iz1}_{power_scaling_factor}/Figures"
SAVE_DATA_DIR = f"multiPhysics_PyGan_24UOX_cell/{numericalMethod}/{voidFractionCorrel}_{frfaccorel}_{P2Pcorel}/mesh{Iz1}_{power_scaling_factor}/Data"

a=os.path.exists(SAVE_FIG_DIR)
if a==False:
    os.makedirs(SAVE_FIG_DIR)

print(SAVE_DIR)
print(f"Keffs = {Keffs}")
# Plot the results
quickPlot(range(len(Keffs)), Keffs, "Keff convergence", "iteration", "Keff", "Keff_convergence.png", path, SAVE_FIG_DIR)

quickPlot(z_mesh, TeffFuel, "Fuel temperature convergence", "height (m)", "TFuel (K)", f"TFuel_convergence_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
quickPlot(z_mesh, Twater, "Coolant temperature convergence", "height (m)", "TCool (K)", f"TCool_convergence_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
quickPlot(z_mesh, rho, "Coolant density convergence", "height (m)", "DCool (kg/m3)", f"DCool_convergence_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
quickPlot(z_mesh, voidFraction, "Void fraction convergence", "height (m)", "Void fraction", f"Void_fraction_convergence_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
quickPlot(z_mesh, Power_Distribs, "Power distribution convergence", "height (m)", "Power (kW)", f"Power_distribution_convergence_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
quickPlot(z_mesh, Volumic_Powers, "Qfiss convergence", "height (m)", "Power (W/m^3)", f"Qfiss_convergence_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
    
quickPlot(z_mesh, Residuals_TeffFuel, "Fuel temperature residuals", "height (m)", "Residuals", f"TFuel_residuals_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
quickPlot(z_mesh, Residuals_Twater, "Coolant temperature residuals", "height (m)", "Residuals", f"TCool_residuals_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
quickPlot(z_mesh, Residuals_rho, "Coolant density residuals", "height (m)", "Residuals", f"DCool_residuals_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
quickPlot(z_mesh, Residuals_voidFraction, "Void fraction residuals", "height (m)", "Residuals", f"Void_fraction_residuals_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
quickPlot(z_mesh, Residuals_Power_Distribs, "Power distribution residuals", "height (m)", "Residuals", f"Power_distribution_residuals_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
#print(f"Residuals Power distribution : {Residuals_Power_Distribs}")
quickPlot(z_mesh, Residuals_Volumic_Powers, "Axial power shape residuals", "height (m)", "Residuals", f"Qfiss_residuals_{relaxPOW_id}_{relaxTH_id}.png", path, SAVE_FIG_DIR)
print("$$$ - multiPhysics.py : END of PLOTTING - $$$")


# 10.) Save the results for exportation to Serpent/OpenFoam/GeN-Foam
# Save the results in a file

a=os.path.exists(SAVE_DATA_DIR)
if a==False:
    os.makedirs(SAVE_DATA_DIR)
os.chdir(SAVE_DATA_DIR)
case = compo_name.split("_")[2]
TeffFuel = np.array(TeffFuel)
Twater = np.array(Twater)
rho = np.array(rho)
Qfiss = np.array(Volumic_Powers)
Power_Distrib = np.array(Power_Distribs)
Residuals_TeffFuel = np.array(Residuals_TeffFuel)
Residuals_Twater = np.array(Residuals_Twater)
Residual_rho = np.array(Residuals_rho)
Residual_Qfiss = np.array(Residuals_Volumic_Powers)
Residuals_Power_Distrib = np.array(Residuals_Power_Distribs)



np.savetxt(f"TeffFuel_{case}_mesh{Iz1}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", TeffFuel[-1])
np.savetxt(f"Keffs_{case}_mesh{Iz1}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Keffs)
np.savetxt(f"Twater_{case}_mesh{Iz1}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Twater[-1])
np.savetxt(f"rho_{case}_mesh{Iz1}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", rho[-1])
np.savetxt(f"voidFraction_{case}_mesh{Iz1}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", voidFraction[-1])
np.savetxt(f"Qfiss_{case}_mesh{Iz1}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Qfiss[-1])
np.savetxt(f"Power_Distrib_{case}_mesh{Iz1}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Power_Distrib[-1])

np.savetxt(f"Residuals_TeffFuel_{case}_mesh{Iz1}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Residuals_TeffFuel[-1])
np.savetxt(f"Residuals_Twater_{case}_mesh{Iz1}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Residuals_Twater[-1])
np.savetxt(f"Residuals_rho_{case}_mesh{Iz1}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Residual_rho[-1])
np.savetxt(f"Residuals_voidFraction_{case}_mesh{Iz1}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Residuals_voidFraction[-1])
np.savetxt(f"Residuals_Qfiss_{case}_mesh{Iz1}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Residual_Qfiss[-1])
np.savetxt(f"Residuals_Power_Distrib_{case}_mesh{Iz1}_{numericalMethod}_{voidFractionCorrel}_{relaxPOW_id}_{relaxTH_id}.txt", Residuals_Power_Distrib[-1])

# 11.) Save power axial distribution and flux axial distribution for the last iteration
os.chdir(path)
shutil.copyfile("Flux01.res", f"{SAVE_DATA_DIR}/Flux01.res")
shutil.copyfile("Flux02.res", f"{SAVE_DATA_DIR}/Flux02.res")
shutil.copyfile("Pdistr.res", f"{SAVE_DATA_DIR}/Pdistr.res")
shutil.copyfile("multiPhysics.result", f"{SAVE_DATA_DIR}/multiPhysics_{relaxPOW_id}_{relaxTH_id}.result")

#
print("$$$ - multiPhysics.py : END of EXPORTS - $$$")
last_time = time.time()
time_for_exports = (last_time - end_iters_time)
print(f"Time spent in exporting = {time_for_exports} s")
total_time = (last_time - start_time)
print(f"Total time spent in multiPhysics.py = {total_time} s")
print(relax_Pow, relax_TH)
print(relaxPOW_id, relaxTH_id)
print("$$$ - multiPhysics.py : END OF SCRIPT - $$$")