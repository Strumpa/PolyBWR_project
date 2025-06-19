# Python3 script to create BWR assembly geometry and test out the SYBILT prototype
# Author: R. Guasch
# Date : 04/05/2025
# Dependencies : GEOM.py

from GEOM import GEO
from CARCEL import CARCEL
from CAR2D import CAR2D
import numpy as np


### Assembly scale data :
assembly_pitch = 15.24 # Assembly pitch
Water_gap =  0.75 # water gap thickness

Box_out = 13.74 # Outer side box outer 
Box_in = 13.4  # Inner side box outer
Box_thi = (Box_out-Box_in)/2  # Thickness

Chan_out = 3.5 # Channel box outer side
Chan_in = 3.34 # Channel box inner side
Chan_thi =  (Chan_out - Chan_in)/2 # Channel box thickness

### Fuel pin scale data :
cell_pitch = 1.295
clad_outer_rad =  0.5140
clad_inner_rad = 0.4520
fuel_rad = 0.4435
UOX_radii = np.array([0.0, 0.313602, 0.396678, 0.43227, fuel_rad, clad_inner_rad, clad_outer_rad])
Gd_radii = np.array([0.0, 0.19834, 0.28049, 0.34353, 0.39668, 0.43227, fuel_rad, clad_inner_rad, clad_outer_rad])

color_map = {"C1": 'khaki', 
             "C2": 'steelblue', 
             "C3": 'hotpink', 
             "C4": 'orange', 
             "C5": 'green', 
             "C6": 'plum', 
             "C7": 'turquoise', 
             "C8": 'crimson'}
# Create a 2D cartesian geometry
BWR_Lattice = CAR2D("BWR_Lattice", level=1, nx=3, ny=3, nz=1, 
                    meshx=np.array([0.0, 4*cell_pitch, 7*cell_pitch, 10*cell_pitch]), 
                    meshy=np.array([0.0, 4*cell_pitch, 7*cell_pitch, 10*cell_pitch]), 
                    meshz=np.array([0.0, 1.0]))
# Create sub geometries
# Create the channel box
BWR_CHANNEL_BOX = CAR2D("BWR_CHANNEL_BOX", level=2, nx=1, ny=1, nz=1,
                    meshx=np.array([0.0, 3*cell_pitch]),
                    meshy=np.array([0.0, 3*cell_pitch]),
                    meshz=np.array([0.0, 1.0]))

# South West
# C01 in AT10_ASSBLY DRAGON5 case
BWR_4x4_SW = CAR2D("BWR_4x4_SW", level=2, nx=4, ny=4, nz=1,
                    meshx=np.array([0.0, cell_pitch, 2*cell_pitch, 3*cell_pitch, 4*cell_pitch]), 
                    meshy=np.array([0.0, cell_pitch, 2*cell_pitch, 3*cell_pitch, 4*cell_pitch]), 
                    meshz=np.array([0.0, 1.0]))
                        # Left to right
CELL_MAP_C01 = np.array([["C11", "C21", "C31", "C51"], # Bottom
                         ["C21", "C41", "C71", "C61"],
                         ["C31", "C71", "C61", "C61"],
                         ["C51", "C61", "C61", "C61"]]) # Top
# South Middle
# C02 in AT10_ASSBLY DRAGON5 case
BWR_3x4_SM = CAR2D("BWR_3x4_SM", level=2, nx=3, ny=4, nz=1,
                    meshx=np.array([0.0, cell_pitch, 2*cell_pitch, 3*cell_pitch]),
                    meshy=np.array([0.0, cell_pitch, 2*cell_pitch, 3*cell_pitch, 4*cell_pitch]),
                    meshz=np.array([0.0, 1.0]))
                        # Left to right
CELL_MAP_C02 = np.array([["C62", "C52", "C42"], # Bottom
                         ["C72", "C62", "C62"],
                         ["C62", "C72", "C62"],
                         ["C62", "C62", "C72"]]) # Top
# South East
# C03 in AT10_ASSBLY DRAGON5 case
BWR_3x4_SE = CAR2D("BWR_3x4_SE", level=2, nx=3, ny=4, nz=1,
                    meshx=np.array([0.0, cell_pitch, 2*cell_pitch, 3*cell_pitch]),
                    meshy=np.array([0.0, cell_pitch, 2*cell_pitch, 3*cell_pitch, 4*cell_pitch]),
                    meshz=np.array([0.0, 1.0]))
                        # Left to right
CELL_MAP_C03 = np.array([["C33", "C23", "C13"], # Bottom
                         ["C73", "C43", "C23"],
                         ["C63", "C73", "C33"],
                         ["C63", "C53", "C43"]]) # Top

# Middle West = Transpose of C02 = BWR_3x4_SM
BWR_4x3_MW = CAR2D("BWR_4x3_MW", level=2, nx=4, ny=3, nz=1,
                    meshx=np.array([0.0, cell_pitch, 2*cell_pitch, 3*cell_pitch, 4*cell_pitch]),
                    meshy=np.array([0.0, cell_pitch, 2*cell_pitch, 3*cell_pitch]),
                    meshz=np.array([0.0, 1.0]))
# Middle East
# C04 in AT10_ASSBLY DRAGON5 case 
BWR_3x3_ME = CAR2D("BWR_3x3_ME", level=2, nx=3, ny=3, nz=1,
                    meshx=np.array([0.0, cell_pitch, 2*cell_pitch, 3*cell_pitch]),
                    meshy=np.array([0.0, cell_pitch, 2*cell_pitch, 3*cell_pitch]),
                    meshz=np.array([0.0, 1.0]))
                        # Left to right
CELL_MAP_C04 = np.array([["C44", "C64", "C44"], # Bottom
                         ["C34", "C74", "C44"],
                         ["C44", "C44", "C44"]]) # Top

# North West
# Transpose of C03 = BWR_3x4_SE
BWR_4x3_NW = CAR2D("BWR_4x3_NW", level=2, nx=4, ny=3, nz=1,
                    meshx=np.array([0.0, cell_pitch, 2*cell_pitch, 3*cell_pitch, 4*cell_pitch]),
                    meshy=np.array([0.0, cell_pitch, 2*cell_pitch, 3*cell_pitch]),
                    meshz=np.array([0.0, 1.0]))
# North Middle
# Transpose of C04 = BWR_3x3_ME
BWR_3x3_NM = CAR2D("BWR_3x3_NM", level=2, nx=3, ny=3, nz=1,
                    meshx=np.array([0.0, cell_pitch, 2*cell_pitch, 3*cell_pitch]),
                    meshy=np.array([0.0, cell_pitch, 2*cell_pitch, 3*cell_pitch]),
                    meshz=np.array([0.0, 1.0]))
# North East
# C05 in AT10_ASSBLY DRAGON5 case
BWR_3x3_NE = CAR2D("BWR_3x4_NE", level=2, nx=3, ny=3, nz=1,
                    meshx=np.array([0.0, cell_pitch, 2*cell_pitch, 3*cell_pitch]),
                    meshy=np.array([0.0, cell_pitch, 2*cell_pitch, 3*cell_pitch]),
                    meshz=np.array([0.0, 1.0]))
                # Left to right
CELL_MAP_C05 = np.array([["C45", "C85", "C35"],# Bottom
                         ["C85", "C45", "C25"],
                         ["C35", "C25", "C15"]])# Top
                   

# Fill 2nd level geometries with the fuel pins = add_CARCEL function

# Adding cells to bottom left corner = SW = C01
CELL_MAP_SW = CELL_MAP_C01.flatten()
for cell_idx in range(len(CELL_MAP_SW)):
    if CELL_MAP_SW[cell_idx][1] == "7" or CELL_MAP_SW[cell_idx][1] == "8":
        cell_radii = Gd_radii
    else:
        cell_radii = UOX_radii
    BWR_4x4_SW.add_CARCEL(CELL_MAP_SW[cell_idx], host_region=cell_idx+1, radii = cell_radii, meshx = np.array([0.0, cell_pitch]), meshy= np.array([0.0, cell_pitch]))

# Adding cells to bottom middle = SM = C02
CELL_MAP_SM = CELL_MAP_C02.flatten()
for cell_idx in range(len(CELL_MAP_SM)):
    if CELL_MAP_SM[cell_idx][1] == "7" or CELL_MAP_SM[cell_idx][1] == "8":
        cell_radii = Gd_radii
    else:
        cell_radii = UOX_radii
    BWR_3x4_SM.add_CARCEL(CELL_MAP_SM[cell_idx], host_region=cell_idx+1, radii = cell_radii, meshx = np.array([0.0, cell_pitch]), meshy= np.array([0.0, cell_pitch]))

# Adding cells to bottom right corner = SE = C03
CELL_MAP_SE = CELL_MAP_C03.flatten()
for cell_idx in range(len(CELL_MAP_SE)):
    if CELL_MAP_SE[cell_idx][1] == "7" or CELL_MAP_SE[cell_idx][1] == "8":
        cell_radii = Gd_radii
    else:
        cell_radii = UOX_radii
    BWR_3x4_SE.add_CARCEL(CELL_MAP_SE[cell_idx], host_region=cell_idx+1, radii = cell_radii, meshx = np.array([0.0, cell_pitch]), meshy= np.array([0.0, cell_pitch]))

# Adding cells to middle right = ME = C04
CELL_MAP_ME = CELL_MAP_C04.flatten()
for cell_idx in range(len(CELL_MAP_ME)):
    if CELL_MAP_ME[cell_idx][1] == "7" or CELL_MAP_ME[cell_idx][1] == "8":
        cell_radii = Gd_radii
    else:
        cell_radii = UOX_radii
    BWR_3x3_ME.add_CARCEL(CELL_MAP_ME[cell_idx], host_region=cell_idx+1, radii = cell_radii, meshx = np.array([0.0, cell_pitch]), meshy= np.array([0.0, cell_pitch]))

# Adding cells to Top right corner = NE = C05
CELL_MAP_NE = CELL_MAP_C05.flatten()
for cell_idx in range(len(CELL_MAP_NE)):
    if CELL_MAP_NE[cell_idx][1] == "7" or CELL_MAP_NE[cell_idx][1] == "8":
        cell_radii = Gd_radii
    else:
        cell_radii = UOX_radii
    BWR_3x3_NE.add_CARCEL(CELL_MAP_NE[cell_idx], host_region=cell_idx+1, radii = cell_radii, meshx = np.array([0.0, cell_pitch]), meshy= np.array([0.0, cell_pitch]))

# Fill on the transpose :
# Adding cells to top left corner = NW, transpose of SE = C03
CELL_MAP_NW = CELL_MAP_C03.transpose().flatten()
for cell_idx in range(len(CELL_MAP_NW)):
    if CELL_MAP_NW[cell_idx][1] == "7" or CELL_MAP_NW[cell_idx][1] == "8":
        cell_radii = Gd_radii
    else:
        cell_radii = UOX_radii
    BWR_4x3_NW.add_CARCEL(CELL_MAP_NW[cell_idx], host_region=cell_idx+1, radii = cell_radii, meshx = np.array([0.0, cell_pitch]), meshy= np.array([0.0, cell_pitch]))
# Adding cells to middle left = MW, transpose of SM = C02
CELL_MAP_MW = CELL_MAP_C02.transpose().flatten()
for cell_idx in range(len(CELL_MAP_MW)):
    if CELL_MAP_MW[cell_idx][1] == "7" or CELL_MAP_MW[cell_idx][1] == "8":
        cell_radii = Gd_radii
    else:
        cell_radii = UOX_radii
    BWR_4x3_MW.add_CARCEL(CELL_MAP_MW[cell_idx], host_region=cell_idx+1, radii = cell_radii, meshx = np.array([0.0, cell_pitch]), meshy= np.array([0.0, cell_pitch]))

# Adding cells to top middle = NM, transpose of ME = C04
CELL_MAP_NM = CELL_MAP_C04.transpose().flatten()
for cell_idx in range(len(CELL_MAP_NM)):
    if CELL_MAP_NM[cell_idx][1] == "7" or CELL_MAP_NM[cell_idx][1] == "8":
        cell_radii = Gd_radii
    else:
        cell_radii = UOX_radii
    BWR_3x3_NM.add_CARCEL(CELL_MAP_NM[cell_idx], host_region=cell_idx+1, radii = cell_radii, meshx = np.array([0.0, cell_pitch]), meshy= np.array([0.0, cell_pitch]))

BWR_Lattice.add_cell(BWR_4x4_SW, host_region=1)
BWR_Lattice.add_cell(BWR_3x4_SM, host_region=2)
BWR_Lattice.add_cell(BWR_3x4_SE, host_region=3)
BWR_Lattice.add_cell(BWR_4x3_MW, host_region=4)
BWR_Lattice.add_cell(BWR_CHANNEL_BOX, host_region=5)
BWR_Lattice.add_cell(BWR_3x3_ME, host_region=6)
BWR_Lattice.add_cell(BWR_4x3_NW, host_region=7)
BWR_Lattice.add_cell(BWR_3x3_NM, host_region=8)
BWR_Lattice.add_cell(BWR_3x3_NE, host_region=9)

#BWR_Lattice.describeGeo()
#BWR_Lattice.buildConnectivityMap()
BWR_Lattice.plotter(color_map=color_map)












