import numpy as np

def get_Initial_CosineProfile(PFiss, Iz):
    """
    PFiss : real number initial power of the fuel rod
    Iz : integer representing the number of axial nodes in the fuel rod (Finite volumes for THM solution) 
    """
    profile = np.zeros(Iz)
    for i in range(Iz):
        profile[i] = PFiss*np.cos(np.pi*i/Iz)
    return profile

def getTHMParams():
    """
    Tfuel : real number fuel temperature
    TCool : real number coolant temperature
    DensCool : real number coolant density
    """
    Tfuel = 750.0
    TCool = 559.0
    DensCool = 0.7803

    return [Tfuel, TCool, DensCool]

def getGeoParams():
    """
    pitch : real number fuel assembly pitch
    height : real number fuel rod height
    """
    pitch = 1.295
    height = 380.0

    return [pitch, height]

def getFuelProps():
    """
    volumic_mass_UOX : real number fuel volumic mass
    Fuel_volume : real number fuel volume
    Fuel_mass : real number fuel mass
    specificPower : real number specific power
    PFiss : real number fission power
    FuelProps : numpy array with two real numbers : PFiss, Fuel_mass
    """
    volumic_mass_UOX = 10970.0 # kg/m3
    fuelRadius = 0.4435e-2
    height = 380.0
    Fuel_volume = np.pi*fuelRadius**2*height # m3
    Fuel_mass = Fuel_volume*volumic_mass_UOX # kg
    specificPower = 38.6 # W/g
    PFiss = specificPower*Fuel_mass*1000.0 # W

    return [PFiss, Fuel_mass]

def getInitialParams(Iz):
    """
    FuelProps : array with two real numbers : PFiss, Fuel_mass
    axPow : array with I_z real numbers representing the initial power profile of the fuel rod
    THParams :  array with three real numbers : Tfuel, TCool, DensCool
    GeoParams : array with two real numbers : pitch, height
    """

    FuelProps = getFuelProps()
    axPow = get_Initial_CosineProfile(FuelProps[0], Iz)
    THParams = getTHMParams()
    GeoParams = getGeoParams()

    return FuelProps, axPow, THParams, GeoParams
