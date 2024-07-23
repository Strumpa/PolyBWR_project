# Python3 script to creates geometries and test out the SYBILT prototype
# Author: R. Guasch
# Date : 18 July 2024
# Dependencies : GEOM.py

from GEOM import *
import numpy as np

# Create a 2D cartesian cell geometry
                                                           # 0.0, Rcomb1,   Rcomb2,   Rcomb3,  Rcomb4, Rgap,   Rclad     
AT10_C11 = CARCEL("C1", level = 3, lr = 6, radii = np.array([0.0, 0.313602, 0.396678, 0.43227, 0.4435, 0.4520, 0.5140]), meshx = np.array([0.0, 1.295]), meshy = np.array([0.0, 1.295]))
AT10_C12 = CARCEL("C1", level = 3, lr = 6, radii = np.array([0.0, 0.313602, 0.396678, 0.43227, 0.4435, 0.4520, 0.5140]), meshx = np.array([0.0, 1.295]), meshy = np.array([0.0, 1.295]))
AT10_C21 = CARCEL("C2", level = 3, lr = 6, radii = np.array([0.0, 0.313602, 0.396678, 0.43227, 0.4435, 0.4520, 0.5140]), meshx = np.array([0.0, 1.295]), meshy = np.array([0.0, 1.295]))
AT10_C22 = CARCEL("C2", level = 3, lr = 6, radii = np.array([0.0, 0.313602, 0.396678, 0.43227, 0.4435, 0.4520, 0.5140]), meshx = np.array([0.0, 1.295]), meshy = np.array([0.0, 1.295]))
#AT10_C22.describeGeo()
#AT10_C22.plotGeo()


# Create a 2D cartesian geometry

#AT10_3x3_cartesian = CAR2D("AT10_3x3", level = 1, lx = 3, ly = 3, lz = 1, meshx = np.array([0.0, 1.295, 2.0*1.295, 3.0*1.295]), meshy = np.array([0.0, 1.0*1.295, 2.0*1.295, 3.0*1.295]), meshz = np.array([0.0, 1.0]))
#AT10_3x3_cartesian.describeGeo()

Cartesian_2x2 = CAR2D("2x2", level=2, lx=2, ly=2, lz=1, meshx=np.array([0.0, 1.295, 2.0*1.295]), meshy=np.array([0.0, 1.0*1.295, 2*1.295]), meshz=np.array([0.0, 1.0]))
#Cartesian_2x2.describeGeo()

# Add cells

Cartesian_2x2.add_cell(AT10_C11, host_region=1)
Cartesian_2x2.add_cell(AT10_C21, host_region=2)
Cartesian_2x2.add_cell(AT10_C22, host_region=3)
Cartesian_2x2.add_cell(AT10_C12, host_region=4)

#Cartesian_2x2.describeGeo()
Cartesian_2x2.buildConnectivityMap()

#Cartesian_1x2 = CAR2D("1x2", level=1, lx=1, ly=2, lz=1, meshx=np.array([0.0, 1.295]), meshy=np.array([0.0, 1.0*1.295, 2*1.295]), meshz=np.array([0.0, 1.0]))
#Cartesian_1x2.describeGeo()

Cartesian_main_geom = CAR2D("Main", level=1, lx=1, ly=1, lz=1, meshx=np.array([0.0, 2.0*1.295]), meshy=np.array([0.0, 2.0*1.295]), meshz=np.array([0.0, 1.0]))
Cartesian_main_geom.add_cell(Cartesian_2x2, host_region=1)

Cartesian_main_geom.buildConnectivityMap()
