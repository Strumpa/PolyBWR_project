### PyGan script for Computing Reactivity worth coefficients using different SelfSheilding methods
# Goal : Validation of multi-physics workflow : DRAGON5 lattice calculations vs Serpent2 reference.


# IMPORT MODULES
import re
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import lifo
import lcm
import cle2000
import serpentTools as st
import pandas as pd

# --- CELL
# GEOMETRY procedure
from cell_geometry import GEO_C_SALT
# LIBRARY creation procedure
from mix_handling_cell import LIB_24UOX
# TRACKING procedure
from tracking_operator_cell import TRK_C_SALT
# SELF-SHIELDING procedure
from self_shielding import selfShieldingUSS
# FLUX CALCULATION procedure
from flux_calculation import fluxCalculationMOC

## DRAGON5 calculation options
draglib = "J311_295"
ssh_method = "RSE"

## Generate Geometry and tracking :

GEOM, GEOM_SSH = GEO_C_SALT()
TRK, TF_EXC, TRK_SSH, TF_EXC_SSH = TRK_C_SALT(GEOM, GEOM_SSH) # CP = Collision Probability, MOC = Method of Characteristics


## IAPWS talbes for water properties
from iapws import IAPWS97
# --- DATA
TFuel_NOM = 750.0 # Fuel temperature in K
TFuel_DOP = 1800.0 # Doppler temperature in K
TCool_inlet = 543.15 # Inlet Coolant temperature in K
TCool_SAT = 561.0 # Saturation Coolant temperature in K
P_outlet = 7.2e6 # Outlet pressure in Pa
estimated_pressure_drop = 0.04e6 # Estimated pressure drop in Pa (P_inlet - P_outlet)

# Get density in inlet :
water_inlet = IAPWS97(T=TCool_inlet, P=(P_outlet+estimated_pressure_drop)/10**6) # Convert Pa to MPa
inlet_density = water_inlet.rho
print(f"inlet_density = {inlet_density}")
# Assume outlet water is at P_outlet and at saturation temperature
water_outlet = IAPWS97(x=0.75, P=P_outlet/10**6) # Convert Pa to MPa, estimate the void fraction to be 0.75 ?
outlet_density = water_outlet.rho
TCool_outlet_iapws = water_outlet.T
print(f"WARNING : Set Saturation pressure is {TCool_SAT}, Saturation pressure computed for x=0.75 and outlet pressure by IAPWS is {TCool_outlet_iapws}")
print(f"outlet_density = {outlet_density}")
### Create LIB: statepoint for inlet water density at nominal fuel temperature:
pyLIB_inlet, N_H_in, N_O_in = LIB_24UOX(draglib, ssh_method, TFuel_NOM, TCool_inlet, inlet_density)
### 
print(f"NOM LIB : TFuel = {TFuel_NOM}, TCool = {TCool_inlet}, N_H = {N_H_in}, N_O = {N_O_in}")

# perform self-shielding
pySSH_lib_inlet =  selfShieldingUSS(pyLIB_inlet, TRK_SSH, TF_EXC_SSH)

# perform flux calculation
keff_inlet, flux_inlet = fluxCalculationMOC(TRK, TF_EXC, pySSH_lib_inlet)
print(f"keff_inlet = {keff_inlet}")


## Create LIB: statepoint at estimated outlet density to estimate reactivity worth of voiding along the channel
pyLIB_void, N_H_out, N_O_out = LIB_24UOX(draglib, ssh_method, TFuel_NOM, TCool_inlet, outlet_density)
print(f"VOID LIB : TFuel = {TFuel_NOM}, TCool = {TCool_inlet}, N_H = {N_H_out}, N_O = {N_O_out}")

# perform self-shielding
pySSH_lib_void =  selfShieldingUSS(pyLIB_void, TRK_SSH, TF_EXC_SSH)

# perform flux calculation
keff_void, flux_void =  fluxCalculationMOC(TRK, TF_EXC, pySSH_lib_void)
print(f"keff_void = {keff_void}")
delta_rho_void = (1/keff_inlet - 1/keff_void)*1e5
print(f"Reactivity worth of VOID : {delta_rho_void} pcm")


## Create LIB: statepoint at inlet conditions, but TFuel - TFuel_DOP to estimate reactivity worth of Doppler effect
pyLIB_DOP, N_H_in, N_O_in = LIB_24UOX(draglib, ssh_method, TFuel_DOP, TCool_inlet, inlet_density)
print(f"DOP LIB : TFuel = {TFuel_DOP}, TCool = {TCool_inlet}, N_H = {N_H_in}, N_O = {N_O_in}")

# perform self-shielding
pySSH_lib_DOP = selfShieldingUSS(pyLIB_DOP, TRK_SSH, TF_EXC_SSH)

# perform flux calculation
keff_DOP, flux_DOP =  fluxCalculationMOC(TRK, TF_EXC, pySSH_lib_DOP)
print(f"keff_DOP = {keff_DOP}")
delta_rho_DOP = (1/keff_inlet - 1/keff_DOP)*1e5
print(f"Reactivity worth of DOP : {delta_rho_DOP} pcm")

## Create LIB: statepoint at oulet conditions, assumin uniform fuel temperature
# Assume outlet coolant temperature is saturation temperature, assuming pressure is outlet pressure and void fraction is 0.75
pyLIB_outlet, N_H_out, N_O_out = LIB_24UOX(draglib, ssh_method, TFuel_NOM, TCool_outlet_iapws, outlet_density)
print(f"OUTLET LIB : TFuel = {TFuel_NOM}, TCool = {TCool_outlet_iapws}, N_H = {N_H_out}, N_O = {N_O_out}")

# perform self-shielding
pySSH_lib_outlet =  selfShieldingUSS(pyLIB_outlet, TRK_SSH, TF_EXC_SSH)

# perform flux calculation
keff_outlet, flux_outlet =  fluxCalculationMOC(TRK, TF_EXC, pySSH_lib_outlet)
print(f"keff_outlet = {keff_outlet}")
delta_rho_outlet = (1/keff_inlet - 1/keff_outlet)*1e5
print(f"Reactivity worth of outlet vs void : {delta_rho_outlet} pcm")


### LOAD Serpent2 results

res = st.read(f"{os.environ['SERPENT_RESULTS']}/AT10_24UOX_NOM_mc_res.m")
keff_nom_S2 = res.resdata["absKeff"].T[0]

res = st.read(f"{os.environ['SERPENT_RESULTS']}/AT10_24UOX_VOID_mc_res.m")
keff_void_S2 = res.resdata["absKeff"].T[0] 

res = st.read(f"{os.environ['SERPENT_RESULTS']}/AT10_24UOX_DOP_mc_res.m")
keff_DOP_S2 = res.resdata["absKeff"].T[0]

res = st.read(f"{os.environ['SERPENT_RESULTS']}/AT10_24UOX_OUT_mc_res.m")
keff_OUT_S2 = res.resdata["absKeff"].T[0]

# Compute S2 effects on reactivity
delta_rho_void_S2 = (1/keff_nom_S2 - 1/keff_void_S2)*1e5
delta_rho_DOP_S2 = (1/keff_nom_S2 - 1/keff_DOP_S2)*1e5
delta_rho_outlet_S2 = (1/keff_nom_S2 - 1/keff_OUT_S2)*1e5

print(f"Serpent2 reactivity worth of VOID : {delta_rho_void_S2} pcm")
print(f"Serpent2 reactivity worth of DOP : {delta_rho_DOP_S2} pcm")
print(f"Serpent2 reactivity worth of OUT ie due to moderator temperature increase : {delta_rho_outlet_S2} pcm")
# in units of pcm / K or pcm / kg/m^3
print(f"Serpent2 reactivity worth of VOID : {delta_rho_void_S2 / (inlet_density-outlet_density)} pcm / (kg/m^3)")
print(f"Serpent2 reactivity worth of DOP : {delta_rho_DOP_S2 / (TFuel_DOP - TFuel_NOM)} pcm / K(fuel)")
print(f"Serpent2 reactivity worth of OUT ie due to moderator temperature increase : {delta_rho_outlet_S2 / (TCool_outlet_iapws - TCool_inlet)} pcm / K(water)")


# D5-S2 differences on keff
delta_keff_nom = (keff_inlet - keff_nom_S2) * 1e5 
delta_keff_void = (keff_void - keff_void_S2) * 1e5
delta_keff_DOP = (keff_DOP - keff_DOP_S2) * 1e5
delta_keff_OUT = (keff_outlet - keff_OUT_S2) * 1e5
print(f"Difference in keff nominal : {delta_keff_nom} pcm")
print(f"Difference in keff VOID : {delta_keff_void} pcm")
print(f"Difference in keff DOP : {delta_keff_DOP} pcm")
print(f"Difference in keff OUT : {delta_keff_OUT} pcm")

# D5-S2 differences on reactivity worth
delta_rho_void_diff = delta_rho_void - delta_rho_void_S2
delta_rho_DOP_diff = delta_rho_DOP - delta_rho_DOP_S2
delta_rho_outlet_diff = delta_rho_outlet - delta_rho_outlet_S2
print(f"Difference in reactivity worth of VOID : {delta_rho_void_diff} pcm")
print(f"Difference in reactivity worth of DOP : {delta_rho_DOP_diff} pcm")
print(f"Difference in reactivity worth of OUT-IN : {delta_rho_outlet_diff} pcm")

print(f"Difference in reactivity worth of VOID : {delta_rho_void_diff / (inlet_density-outlet_density)} pcm / (kg/m^3)")
print(f"Difference in reactivity worth of DOP : {delta_rho_DOP_diff / (TFuel_DOP - TFuel_NOM)} pcm / K")

## Consider RSE values too 

# store results in a df 
results = {
    "Parameter": ["Keff Inlet", "Keff Void", "Keff DOP", "Keff Outlet",
                  "Reactivity Void", "Reactivity DOP", "Reactivity Outlet",
                  "Delta Keff Nominal", "Delta Keff Void", "Delta Keff DOP", "Delta Keff Outlet",
                  "Delta Reactivity Void", "Delta Reactivity DOP", "Delta Reactivity Outlet"],
    "D5": [keff_inlet, keff_void, keff_DOP, keff_outlet,
           delta_rho_void, delta_rho_DOP, delta_rho_outlet,
           delta_keff_nom, delta_keff_void, delta_keff_DOP, delta_keff_OUT,
           delta_rho_void_diff, delta_rho_DOP_diff, delta_rho_outlet_diff],
    "S2": [keff_nom_S2, keff_void_S2, keff_DOP_S2, keff_OUT_S2,
           delta_rho_void_S2, delta_rho_DOP_S2, delta_rho_outlet_S2,
           None, None, None, None,
           None, None, None]
}

df_results = pd.DataFrame(results)
# Save results to a CSV file
output_file = f"df_AT10_24UOX_DOP_VOID_{draglib}_{ssh_method}_results.csv"
df_results.to_csv(output_file, index=False)

save_dir = "DRAGON5_validation_for_mphys"
if not os.path.exists(save_dir):
    os.makedirs(save_dir)


