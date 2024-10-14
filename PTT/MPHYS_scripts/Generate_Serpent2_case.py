# Python3 script to generate Serpent2 input case
# Purpose : generate Serpent2 3D model of pincell with variable fuel and moderator properties
# Author  : R. Guasch
# Date    : 10/10/2024

import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import DensToIsoDens 



def parse_multiPhysics_output(filename):
    """Parse the output file of the multiPhysics code"""
    data = np.loadtxt(filename)
    return data

def create_geometry(number_axial_slices, pitch, fuel_radius, gap_radius, clad_radius, height, isGd):
    """Create the geometry for the Serpent2 input file"""
    slice_thinkness = height/number_axial_slices
    x_max = pitch/2
    x_min = -pitch/2
    y_max = pitch/2
    y_min = -pitch/2
    pin_radii = computeSantamarinaradii(fuel_radius, gap_radius, clad_radius, isGd)
    print(f"% --- 1. Pin definition")
    for i in range(number_axial_slices):
        if isGd==False:
            print(f"pin fuel_slice{i+1} \n UOx_A_{i+1} {pin_radii[0]:.4f} \n UOx_B_{i+1} {pin_radii[1]:.4f} \n UOx_C_{i+1}  {pin_radii[2]:.4f} \n UOx_D_{i+1} {pin_radii[3]:.4f} \n gap  {pin_radii[4]:.4f} \n clad {pin_radii[5]:.4f} \n H2O_{i+1} \n")
        else:
            print(f"pin fuel_slice{i+1} \n Gd_A_{i+1} {pin_radii[0]:.4f} \n Gd_B_{i+1} {pin_radii[1]:.4f} \n Gd_C_{i+1}  {pin_radii[2]:.4f} \n Gd_D_{i+1} {pin_radii[3]:.4f} \n Gd_E_{i+1} {pin_radii[4]:.4f} \n Gd_F_{i+1} {pin_radii[5]:.4f} \n gap {pin_radii[6]:.4f} \n clad {pin_radii[7]:.4f} \n H2O_{i+1} \n")
    print(f"% --- 2. Surfaces definition")
    for i in range(number_axial_slices):
  
        print(f"surf {i+1} cuboid {x_min} {x_max} {y_min} {y_max} {i*slice_thinkness} {i*slice_thinkness+slice_thinkness}")
    print(f"surf {number_axial_slices+1} cuboid {x_min} {x_max} {y_min} {y_max} 0 {height}")
    print(f"% --- 3. Cell definition")
    for i in range(number_axial_slices):
        
        print(f"cell {i+1} 0 fill fuel_slice{i+1} -{i+1} -{number_axial_slices+1}")
    print(f"cell {number_axial_slices+1} 0 outside {number_axial_slices+1}")

def computeSantamarinaradii(outer_fuel_radius, gap_radius, clad_radius, isGd):
        # Helper to defined fuel region radii for fuel pins --> important for evolution calculations
        # UOX : Finer discretization close to the outer radii is important for rim effects : Pu formation
        # Gd2O3 : rim effect + Gd evolution = challenge for deterministic methods, need 6 sub regions.
        """
        r_out = float, fuel radius
        isGd = bool, adapting discretization to Gd pin --> 6 radial sub regions instead.
        A. Santamarina recommandations :
        volumes for UOX pins : 50%, 80%, 95% and 100%
        volumes for Gd2O3 pins : 20%, 40%, 60%, 80%, 95% and 100%
        """
        if isGd==False:
            pin_radii=[(0.5**0.5)*outer_fuel_radius, (0.8**0.5)*outer_fuel_radius, (0.95**0.5)*outer_fuel_radius, outer_fuel_radius, gap_radius, clad_radius]
            
        else :
            pin_radii=[(0.2**0.5)*outer_fuel_radius, (0.4**0.5)*outer_fuel_radius, (0.6**0.5)*outer_fuel_radius, 
                      (0.8**0.5)*outer_fuel_radius, (0.95**0.5)*outer_fuel_radius, outer_fuel_radius, gap_radius, clad_radius]
        return pin_radii

def create_material_volumes(Teff_fuel, TWater, iso_dens_O, iso_dens_H, isGd):
    """
    Generate Sepent2 material volumes
    Teff_fuel : list of effective temperatures in fuel
    TWater : list of temperatures in water
    iso_dens_O : list of oxygen isotopic densities (in water) obtained from DensToIsoDens(rho)
    iso_dens_H : list of hydrogen isotopic densities (in water) obtained from DensToIsoDens(rho)
    isGd : bool, True if Gd is present in the fuel pin --> then radially subdivide the fuel pin in 6 regions, if False 4 regions.
    """
    if isGd == False:
        fuel_names = ["UOx_A", "UOx_B", "UOx_C", "UOx_D"]
    else:
        fuel_names = ["Gd_A", "Gd_B", "Gd_C", "Gd_D", "Gd_E", "Gd_F"]
    print(f"% --- 4. Materials definition")
    for i in range(len(Teff_fuel)):
        for fuel_name in fuel_names:
            print(f"mat {fuel_name}_{i+1} sum burn 1 tmp {Teff_fuel[i]:.2f}")
            # For now only consider AT10_24UOX cell composition
            if 293.0<Teff_fuel[i]<550.0:
                temp_suffix = "02c"
            elif 550.0<=Teff_fuel[i]<900.0:
                temp_suffix = "05c"
            elif 900.0<=Teff_fuel[i]<1200.0:
                temp_suffix = "09c"
            elif 1200.0<=Teff_fuel[i]<1999.0:
                temp_suffix = "12c"
            elif 1999.0<=Teff_fuel[i]:
                temp_suffix = "20c"
            print(f" 92234.{temp_suffix}   5.15910E-06 \n 92235.{temp_suffix}   5.67035E-04 \n 92238.{temp_suffix}   2.27631E-02 \n 8016.{temp_suffix}   4.66705E-02")
    for i in range(len(TWater)):
        if 293.0<TWater[i]<550.0:
                temp_suffix = "02c"
        elif 550.0<=TWater[i]<900.0:
            temp_suffix = "05c"
        elif 900.0<=TWater[i]<1200.0:
            temp_suffix = "09c"
        elif 1200.0<=TWater[i]<1999.0:
            temp_suffix = "12c"
        elif 1999.0<=TWater[i]:
            temp_suffix = "20c"
        print(f"mat H2O_{i+1} sum tmp {TWater[i]:.02f} moder lwtr 1001")
        print(f" 1001.{temp_suffix}   {iso_dens_H[i]:E} \n 8016.{temp_suffix}   {iso_dens_O[i]:E}")

    for i in range(len(Teff_fuel)):
        if 293.0<Teff_fuel[i]<550.0:
                temp_suffix = "02c"
        elif 550.0<=Teff_fuel[i]<900.0:
            temp_suffix = "05c"
        elif 900.0<=Teff_fuel[i]<1200.0:
            temp_suffix = "09c"
        elif 1200.0<=Teff_fuel[i]<1999.0:
            temp_suffix = "12c"
        elif 1999.0<=Teff_fuel[i]:
            temp_suffix = "20c"
        print(f"mat U235_{i+1} 1.0 tmp {Teff_fuel[i]:.2f} 92235.{temp_suffix} 1.0")
        print(f"mat U238_{i+1} 1.0 tmp {Teff_fuel[i]:.2f} 92238.{temp_suffix} 1.0")
        print(f"mat U234_{i+1} 1.0 tmp {Teff_fuel[i]:.2f} 92234.{temp_suffix} 1.0")

def create_detectors(number_axial_slices, isGd):
    if isGd == False:
        fuel_names = ["UOx_A", "UOx_B", "UOx_C", "UOx_D"]
    else:
        fuel_names = ["Gd_A", "Gd_B", "Gd_C", "Gd_D", "Gd_E", "Gd_F"]
    print(f"% --- 5. Detectors definition")
    print(f"% --- energy grid for the detectors (1g)")
    print("ene 3 1 1.1E-11 1.9640E+1")
    print(f"% --- energy grid for the detectors (2g)")
    print("ene 2 1 1.1E-11 6.25E-7 1.9640E+1")

    

    for i in range(number_axial_slices):
        print(f"det _FLUX_2G_{i+1}") 
        print(f"de 2")
        print(f"dc {i+1} \n")

        print(f"det _FUEL_1G_{i+1}")
        print(f"de 3")
        print(f"dr 102 U235_{i+1}")
        print(f"dr 102 U238_{i+1}")
        print(f"dr 102 U234_{i+1}")
        print(f"dr -6 U235_{i+1}")
        print(f"dr -6 U238_{i+1}")
        print(f"dr -6 U234_{i+1}")
        print(f"dr 301 U235_{i+1}")
        print(f"dr 301 U238_{i+1}")
        print(f"dr 301 U234_{i+1}")
        print(f"dc {i+1} \n")
        
    

number_axial_slices = 10
pow_relax_factor = 0.9
pow_scaling_factor = 2
Teff_fuel = parse_multiPhysics_output(f"../../Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/EPRIvoidModel_Churchill_HEM1/mesh{number_axial_slices}_{pow_scaling_factor}/Data/TeffFuel_24UOX_mesh{number_axial_slices}_BiCG_EPRIvoidModel_relaxedPOW_{pow_relax_factor}_relaxedTH_0.1.txt")
#print(T_data)
TWater = parse_multiPhysics_output(f"../../Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/EPRIvoidModel_Churchill_HEM1/mesh{number_axial_slices}_{pow_scaling_factor}/Data/Twater_24UOX_mesh{number_axial_slices}_BiCG_EPRIvoidModel_relaxedPOW_{pow_relax_factor}_relaxedTH_0.1.txt")
DWater = parse_multiPhysics_output(f"../../Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/EPRIvoidModel_Churchill_HEM1/mesh{number_axial_slices}_{pow_scaling_factor}/Data/rho_24UOX_mesh{number_axial_slices}_BiCG_EPRIvoidModel_relaxedPOW_{pow_relax_factor}_relaxedTH_0.1.txt")

DWater = DWater*1e-3 # kg/m^3 to g/cm^3
iso_dens_H = []
iso_dens_O = []
for rho in DWater:
    iso_dens_H.append(DensToIsoDens.DensToIsoDens(rho)[0])
    iso_dens_O.append(DensToIsoDens.DensToIsoDens(rho)[1])


# check number of elements in each array
if len(Teff_fuel) != len(TWater) or len(TWater) != len(DWater) or len(DWater) != len(iso_dens_H) or len(iso_dens_H) != len(iso_dens_O):
    print("Error: unable to generate Serpent2 case")
    sys.exit(1)
else:
    #print("Data successfully parsed")
    #print(f"Number of axial slices to be included in the geometry: {len(Teff_fuel)}")
    #print(f"Generating geometry")
    # Create Serpent2 geometry input
    create_geometry(number_axial_slices=len(Teff_fuel), pitch=1.295, fuel_radius=0.4435, 
                    gap_radius=0.4520, clad_radius= 0.5140, height=380.0, isGd= False)
    #print(f"Number of material volumes to be generated: {len(Teff_fuel)}")
    #print(f"Generating material volumes")
    # Create Serpent2 material volumes
    create_material_volumes(Teff_fuel, TWater, iso_dens_O, iso_dens_H, isGd=False)

    #print(f"Generating detectors")
    # Create Serpent2 detectors
    create_detectors(number_axial_slices=len(Teff_fuel), isGd=False)

    






