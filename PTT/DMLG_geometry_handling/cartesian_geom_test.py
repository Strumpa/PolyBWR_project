# Python3 script to creates geometries and test out the SYBILT prototype
# Author: R. Guasch
# Date : 18 July 2024
# Dependencies : GEOM.py

from GEOM import GEO
from CARCEL import CARCEL
from CAR2D import CAR2D
import numpy as np

# Create a 2D cartesian cell geometry
                                                           # 0.0, Rcomb1,   Rcomb2,   Rcomb3,  Rcomb4, Rgap,   Rclad     
AT10_C1 = CARCEL("C1", level = 3, nr = 6, radii = np.array([0.0, 0.313602, 0.396678, 0.43227, 0.4435, 0.4520, 0.5140]), meshx = np.array([0.0, 1.295]), meshy = np.array([0.0, 1.295]))
AT10_C21 = CARCEL("C2", level = 3, nr = 6, radii = np.array([0.0, 0.313602, 0.396678, 0.43227, 0.4435, 0.4520, 0.5140]), meshx = np.array([0.0, 1.295]), meshy = np.array([0.0, 1.295]))
AT10_C22 = CARCEL("C2", level = 3, nr = 6, radii = np.array([0.0, 0.313602, 0.396678, 0.43227, 0.4435, 0.4520, 0.5140]), meshx = np.array([0.0, 1.295]), meshy = np.array([0.0, 1.295]))
AT10_C31 = CARCEL("C3", level = 3, nr = 6, radii = np.array([0.0, 0.313602, 0.396678, 0.43227, 0.4435, 0.4520, 0.5140]), meshx = np.array([0.0, 1.295]), meshy = np.array([0.0, 1.295]))
AT10_C32 = CARCEL("C3", level = 3, nr = 6, radii = np.array([0.0, 0.313602, 0.396678, 0.43227, 0.4435, 0.4520, 0.5140]), meshx = np.array([0.0, 1.295]), meshy = np.array([0.0, 1.295]))
AT10_C4 = CARCEL("C4", level = 3, nr = 6, radii = np.array([0.0, 0.313602, 0.396678, 0.43227, 0.4435, 0.4520, 0.5140]), meshx = np.array([0.0, 1.295]), meshy = np.array([0.0, 1.295]))
AT10_C6 = CARCEL("C6", level = 3, nr = 6, radii = np.array([0.0, 0.313602, 0.396678, 0.43227, 0.4435, 0.4520, 0.5140]), meshx = np.array([0.0, 1.295]), meshy = np.array([0.0, 1.295]))

                                                           # 0.0, Rcomb1,  Rcomb2,  Rcomb3,  Rcomb4,  Rcomb5,  Rcomb6, Rgap,   Rclad                                                 
AT10_C71 = CARCEL("C7", level = 3, nr = 8, radii = np.array([0.0, 0.19834, 0.28049, 0.34353, 0.39668, 0.43227, 0.4435, 0.4520, 0.5140]), meshx = np.array([0.0, 1.295]), meshy = np.array([0.0, 1.295]))
AT10_C72 = CARCEL("C7", level = 3, nr = 8, radii = np.array([0.0, 0.19834, 0.28049, 0.34353, 0.39668, 0.43227, 0.4435, 0.4520, 0.5140]), meshx = np.array([0.0, 1.295]), meshy = np.array([0.0, 1.295]))


# Create a 2D cartesian geometry

#AT10_3x3_cartesian = CAR2D("AT10_3x3", level = 1, lx = 3, ly = 3, lz = 1, meshx = np.array([0.0, 1.295, 2.0*1.295, 3.0*1.295]), meshy = np.array([0.0, 1.0*1.295, 2.0*1.295, 3.0*1.295]), meshz = np.array([0.0, 1.0]))
#AT10_3x3_cartesian.describeGeo()

Cartesian_3x3 = CAR2D("3x3", level=2, nx=3, ny=3, nz=1, meshx=np.array([0.0, 1.295, 2.0*1.295, 3*1.295]), meshy=np.array([0.0, 1.0*1.295, 2*1.295, 3*1.295]), meshz=np.array([0.0, 1.0]))
#Cartesian_3x3.describeGeo()

# Add cells to describe the fuel lattice

Cartesian_3x3.add_cell(AT10_C1, host_region=1)
Cartesian_3x3.add_cell(AT10_C21, host_region=2)
Cartesian_3x3.add_cell(AT10_C31, host_region=3)
Cartesian_3x3.add_cell(AT10_C22, host_region=4)
Cartesian_3x3.add_cell(AT10_C4, host_region=5)
Cartesian_3x3.add_cell(AT10_C71, host_region=6)
Cartesian_3x3.add_cell(AT10_C32, host_region=7)
Cartesian_3x3.add_cell(AT10_C72, host_region=8)
Cartesian_3x3.add_cell(AT10_C6, host_region=9)

Cartesian_3x3.describeGeo()
Cartesian_3x3.buildConnectivityMap()
Cartesian_3x3.getSurfaceGeometricalData(11)


#Cartesian_3x3.plotGeo()

#Cartesian_1x2 = CAR2D("1x2", level=1, nx=1, ny=2, nz=1, meshx=np.array([0.0, 1.295]), meshy=np.array([0.0, 1.0*1.295, 2*1.295]), meshz=np.array([0.0, 1.0]))
#Cartesian_1x2.describeGeo()

Cartesian_main_geom = CAR2D("Main", level=1, nx=1, ny=1, nz=1, meshx=np.array([0.0, 3.0*1.295]), meshy=np.array([0.0, 3.0*1.295]), meshz=np.array([0.0, 1.0]))
Cartesian_main_geom.add_cell(Cartesian_3x3, host_region=1)

Cartesian_main_geom.buildConnectivityMap()
Cartesian_main_geom.plotGeo()
