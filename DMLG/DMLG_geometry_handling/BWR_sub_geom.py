# Python3 script to create subset BWR assembly geometry and test out the GEO/SYBILT prototype
# Author: R. Guasch
# Date : 04/05/2025
# Dependencies : GEOM.py

from GEOM import GEO
from CARCEL import CARCEL
from CAR2D import CAR2D
import numpy as np


### Geometric data :
UOX_radii = np.array([0.0, 0.313602, 0.396678, 0.43227, 0.4435, 0.4520, 0.5140])
Gd_radii = np.array([0.0, 0.19834, 0.28049, 0.34353, 0.39668, 0.43227, 0.4435, 0.4520, 0.5140])
cell_pitch = 1.295


### Testing 3 lvl geometry

CUOX_n3_SW = CARCEL("C_UOX_n3", level = 3, nr = 6, radii = UOX_radii, meshx = np.array([0.0, cell_pitch]), meshy = np.array([0.0, cell_pitch]))
CGd_n3_SW = CARCEL("C_Gd_n3", level = 3, nr = 8, radii = Gd_radii, meshx = np.array([0.0, cell_pitch]), meshy = np.array([0.0, cell_pitch]))

CUOX_n3_NE = CARCEL("C_UOX_n3", level = 3, nr = 6, radii = UOX_radii, meshx = np.array([0.0, cell_pitch]), meshy = np.array([0.0, cell_pitch]))
CGd_n3_NE = CARCEL("C_Gd_n3", level = 3, nr = 8, radii = Gd_radii, meshx = np.array([0.0, cell_pitch]), meshy = np.array([0.0, cell_pitch]))

BWR_4x3_n2 = CAR2D("BWR_4x3_n2", level=2, nx=4, ny=3, nz=1,
                    meshx=np.array([0.0, cell_pitch, 2*cell_pitch, 3*cell_pitch, 4*cell_pitch]), 
                    meshy=np.array([0.0, cell_pitch, 2*cell_pitch, 3*cell_pitch]), 
                    meshz=np.array([0.0, 1.0]))

BWR_4x3_n2.add_cell(CUOX_n3_SW, host_region=1)
BWR_4x3_n2.add_cell(CUOX_n3_SW, host_region=2)
BWR_4x3_n2.add_cell(CUOX_n3_SW, host_region=3)
BWR_4x3_n2.add_cell(CUOX_n3_SW, host_region=4)

BWR_4x3_n2.add_cell(CUOX_n3_SW, host_region=5)
BWR_4x3_n2.add_cell(CUOX_n3_SW, host_region=6)
BWR_4x3_n2.add_cell(CGd_n3_SW, host_region=7)
BWR_4x3_n2.add_cell(CUOX_n3_SW, host_region=8)

BWR_4x3_n2.add_cell(CUOX_n3_SW, host_region=9)
BWR_4x3_n2.add_cell(CGd_n3_SW, host_region=10)
BWR_4x3_n2.add_cell(CUOX_n3_SW, host_region=11)
BWR_4x3_n2.add_cell(CUOX_n3_SW, host_region=12)

print(f"number of sub-geometries in BWR_4x3_n2 = {len(BWR_4x3_n2.sub_geometries)}")

BWR_3x4_n2 = CAR2D("BWR_3x4_n2", level=2, nx=3, ny=4, nz=1,
                    meshx=np.array([0.0, cell_pitch, 2*cell_pitch, 3*cell_pitch]), 
                    meshy=np.array([0.0, cell_pitch, 2*cell_pitch, 3*cell_pitch, 4*cell_pitch]), 
                    meshz=np.array([0.0, 1.0]))

BWR_3x4_n2.add_cell(CUOX_n3_NE, host_region=1)
BWR_3x4_n2.add_cell(CUOX_n3_NE, host_region=2)
BWR_3x4_n2.add_cell(CUOX_n3_NE, host_region=3)
BWR_3x4_n2.add_cell(CUOX_n3_NE, host_region=4)
BWR_3x4_n2.add_cell(CGd_n3_NE, host_region=5)
BWR_3x4_n2.add_cell(CUOX_n3_NE, host_region=6)
BWR_3x4_n2.add_cell(CGd_n3_NE, host_region=7)
BWR_3x4_n2.add_cell(CUOX_n3_NE, host_region=8)
BWR_3x4_n2.add_cell(CUOX_n3_NE, host_region=9)
BWR_3x4_n2.add_cell(CUOX_n3_NE, host_region=10)
BWR_3x4_n2.add_cell(CUOX_n3_NE, host_region=11)
BWR_3x4_n2.add_cell(CUOX_n3_NE, host_region=12)

print(f"number of sub-geometries in BWR_3x4_n2 = {len(BWR_3x4_n2.sub_geometries)}")


BWR_geom = CAR2D("BWR_main_n1", level=1, nx=3, ny=3, nz=1,
                    meshx=np.array([0.0, 4*cell_pitch, 7*cell_pitch, 10*cell_pitch]), 
                    meshy=np.array([0.0, 3*cell_pitch, 6*cell_pitch, 10*cell_pitch]), 
                    meshz=np.array([0.0, 1.0]))
BWR_geom.add_cell(BWR_4x3_n2, host_region=1)
BWR_geom.add_cell(BWR_3x4_n2, host_region=9)

print(f"number of sub-geometries in BWR_geom = {len(BWR_geom.sub_geometries)}")
for i in range(len(BWR_geom.sub_geometries)):
    print(f"$$-BWR_geom : sub_geometry {i} : {BWR_geom.sub_geometries[i].name}")
    print(f"number of sub-geometries in {BWR_geom.sub_geometries[i].name} = {len(BWR_geom.sub_geometries[i].sub_geometries)}")


BWR_geom.plotter()

print(f"4* cell pitch : {4*cell_pitch}")
print(f"3* cell pitch : {3*cell_pitch}")



print(f"7* cell pitch : {7*cell_pitch}")
print(f"6* cell pitch : {6*cell_pitch}")
