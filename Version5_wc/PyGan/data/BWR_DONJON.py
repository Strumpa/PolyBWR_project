# PyGan scrip to test DONJON procedures used in multiPhysics solution
# Author: R. Guasch
# Date 19/09/2024

import numpy as np
import os, shutil
import lifo
import lcm
import cle2000
import matplotlib

from INIT_DONJON import *

# User choice:
Tfuel = 750.0
TCool = 559.0
DensCool = 0.7803

fuelRadius = 0.4435e-2
pitch = 1.295
height = 380.0
Iz = 10 # or 20
ene_groups = 2

## Additional parameters needed for the calculation
volumic_mass_UOX = 10970 # kg/m3
Fuel_volume = np.pi*fuelRadius**2*height # m3
Fuel_mass = Fuel_volume*volumic_mass_UOX # kg
specificPower = 38.6 # W/g
PFiss = specificPower*Fuel_mass*1000 # W
print('PFiss = ', PFiss)
print('Fuel_mass = ', Fuel_mass)



