### PyGan script for validation study of AT10_45Gd single pincell benchamrk
## Author : R. Guasch
# Date : 2025/02/25
# Purpose : test and compare performances of Jeff3.1.1 vs ENDF/B-VIII.1 libraries on the SHEM295 energy mesh
# Second purpose : validate self-consistency of ENDF/B-VIII.1 ace file generation through : 
#                  - Njoy2012_up_140
#                  - PyNjoy2016 (IRSN fork) 
#                  - NJOY2016 (official release)

# Importing libraries
import re
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import lifo
import lcm
import cle2000

# Importing PyGan procedures
from assertS import *

# --- CELL
# GEOMETRY class
from GEO_C_SALT import *
from GEO_C_NXT import *
# MIXTURES class
from MIX_C import *
from MIX_NG0 import *
# TRACKING class
from TRK_C_SALT import *
from TRK_C_SYBNXT import *
# BU DEPLETION class
from BU_C import *
from BU_C_SYBNXT import *

# --- OTHERS
from getLists import *


# --- DRAGON5-SERPENT2 comparison : post treatment
from postproc_cst_pow_evol import DRAGON_case as D5_case


# --- Selecting calculation options

# 1) Selecting the geometry tracking options
# Tracking : SALT or SYBNXT
# 
tracking_option = "SALT" #, "SYBNXT"

# 2) Selecting the evaluation used for the simulation
# ENDFb8r1_295 or J311_295
#
evaluation = "endfb8r1" # Jeff3.1.1 or ENDF/B-VIII.1
draglib_name = "endfb8r1_295"

# 3) Selecting the self-shielding method
# RSE, PT or AUTO
#
ssh_option = "RSE" # "RSE" , "PT", "SUBG", "AUTO"

# 4) Selecting the burnup calculation options
# burnup_steps = "UOx", "UOx_autop5", "UOx2_autop5", "UOx4_autop5", "UOx6_autop5" etc
#
burnup_points = "Gd_autop3" # "Gd_autop3", "Gd_autop4", "Gd"

# 5) Selecting the burnup calculation options
# Solver : "RUNG" or "KAPS"
# Saturation : "NODI" or "DIRA" : use eq 3.32 or eq 3.33 (with dirac delta) to treat saturation
# val_exp : list of values to impose saturation for isotopes with lambda * (xtf - xti) >= val_exp
#
solver_option = "RUNG"
saturation_option = "NODI" #, "DIRA"
rates_extr = "NOEX" # "EXTR", "NOEX"

# 6) Selecting the energy deposition options
# Global energy deposition : "NOGL"=only energy release in fuel is used for normalization or "GLOB" = global energy release model, 

glob_opt = "GLOB"

# 7) Select which D5 case to run
exec_D5_no_modif = False # True : run DRAGON5 calculations, False : skip DRAGON5 calculations
exec_D5_no_NG0 = True # True : run DRAGON5 calculations with NG0 depletion chain, False : skip DRAGON5 calculations

######## Options for DRAGON5-SERPENT2 comparison ########
# Select origin of SERPENT2 data + evaluation
#
S2_evaluations = "endfb8r1"#,"Jef3.1.1"]
Njoy_versions = "pynjoy2012" #"PyNjoy2016"

######## 
# Result handling and creation of the results directory
#
tracked_nuclides = ["U235","U238","Pu239","Pu240","Pu241","Pu242","Am241","Xe135","Sm149","Gd155","Gd157"]

#
# Create the results directory
path=os.getcwd()
save_dir_D5 = f"{path}/AT10_45Gd_Cst_pow_evol_results/D5_{evaluation}"
save_dir_S2 = f"{path}/AT10_45Gd_Cst_pow_evol_results/S2_{evaluation}"
save_dir_comparison = f"{path}/AT10_45Gd_Cst_pow_evol_results/Comparison_{evaluation}"
if not os.path.exists(save_dir_D5):
    os.makedirs(save_dir_D5)
if not os.path.exists(save_dir_S2):
    os.makedirs(save_dir_S2)
if not os.path.exists(save_dir_comparison):
    os.makedirs(save_dir_comparison)


print(f"Burnup points : {burnup_points}")
# Recovering ListBU ListAUTOP ListCOMPO
[ListeBU,ListeAUTOP,ListeCOMPO]=getLists(burnup_points)
BU_lists = {"BU": ListeBU, "AUTOP": ListeAUTOP, "COMPO": ListeCOMPO}
# Create Steplist for BU - SELFSHIELDING - COMPO save 
StepList = lcm.new('LCM','burnup_steps')
StepList['ListBU']    = np.array(ListeBU, dtype='f')
StepList['ListAutop'] = np.array(ListeAUTOP, dtype='f')
StepList['ListCompo'] = np.array(ListeCOMPO, dtype='f')
StepList.close() # close without erasing


### BEGIN DRAGON5 calculations with default DEPL structure ###

if tracking_option == "SALT":
    # geometry definition
    GEOM, GEOM_SS = GEO_C_SALT()
    # tracking of geometries
    TRK, TF_EXC, TRK_SS, TF_EXC_SS = TRK_C_SALT(GEOM,GEOM_SS) # CP = Collision Probability, MOC = Method of Characteristics
elif tracking_option == "SYBNXT":
    GEOM, GEOM_SS = GEO_C_NXT()
    TRK, TF_EXC, TRK_SS = TRK_C_SYBNXT(GEOM,GEOM_SS)
else:
    print("Tracking option not recognized")
    sys.exit(1)

if exec_D5_no_modif:
    print(f"Self-shielding option : {ssh_option}")
    LIB = MIX_C(draglib_name,ssh_option)

    print(f"Saturation option : {saturation_option}")
    print(f"Global energy deposition option : {glob_opt}") 
    name_compo = f"_COMPO_AT10_45Gd_{draglib_name}_{ssh_option}_{tracking_option}_{solver_option}_{saturation_option}_{glob_opt}_{burnup_points}"
    if tracking_option == "SALT":
        CPO = BU_C("COMPO", LIB, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, name_compo, ssh_option, solver_option, glob_opt, saturation_option, rates_extr, val_exp = 80.0)
    elif tracking_option == "SYBNXT": 
        CPO = BU_C_SYBNXT("COMPO", LIB, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, name_compo, ssh_option, solver_option, glob_opt, saturation_option, val_exp = 80.0)
    print(f"creating D5 case for {name_compo}")
    D5case = D5_case(CPO, draglib_name, burnup_points, ssh_option, "CORR", saturation_option, solver_option, tracked_nuclides, BU_lists, save_dir_D5)
    D5case.plot_keffs()
    for iso in tracked_nuclides:
        D5case.plot_Ni(iso)


if exec_D5_no_NG0:
    ### Begin calculations with modified DEPL structure ###
        
    # --- Call to DRAGON5 CLE-2000 procedures :
    # --- DRAGON5 microlib generation
    pyLIB_NG0 = MIX_NG0(draglib_name) # Creation of the microlib, default D5 energy deposition mode
    ssh_option = "RSE"
    #
    # names for exportation
    print(f"State of the calculation NG0 : {draglib_name} {ssh_option} {saturation_option} {solver_option}")
    compo_name = f"_CPO_{draglib_name}_NG0_{ssh_option}_{tracking_option}_{burnup_points}_{solver_option}_{saturation_option}_{rates_extr}_{glob_opt}"
    # run DRAGON5 calculation with BU evolution
    if tracking_option == "SALT":
        CPO_NG0 = BU_C("COMPO", pyLIB_NG0, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, compo_name, ssh_option, solver_option, glob_opt, rates_extr, saturation_option)
    elif tracking_option == "SYBNXT": 
        CPO_NG0 = BU_C_SYBNXT("COMPO", pyLIB_NG0, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, compo_name, ssh_option, solver_option, glob_opt, saturation_option)
    print(f"creating D5 case for {compo_name}")
    D5case_NG0 = D5_case(CPO_NG0, draglib_name, burnup_points, ssh_option, "CORR", saturation_option, solver_option, tracked_nuclides, BU_lists, save_dir_D5)
    D5case_NG0.plot_keffs()
    for iso in tracked_nuclides:
        D5case_NG0.plot_Ni(iso)
    # --- Post-processing of DRAGON5 results
        
    #### END SCRIPT #### 
