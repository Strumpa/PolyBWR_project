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
from BU_PCC import *

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
draglib_name = "endfb8r1_295" # "endfb8r1_295", "endfb81295K" : with total KERMA (MT301), "endfb81295K2" : with modified KERMA (MT301-MT318)+MT458 data.

# 3) Selecting the self-shielding method
# RSE, PT or AUTO
#
ssh_option = "PT" # "RSE" , "PT", "SUBG", "AUTO"
corr = "NOCORR" # "CORR" or "NOCORR"

# 4) Selecting the burnup calculation options
# burnup_steps = "UOx", "UOx_autop5", "UOx2_autop5", "UOx4_autop5", "UOx6_autop5" etc
#
burnup_points = "Gdf_BOC" # "Gd_autop3", "Gd_autop4", "Gd", "Gd2_autop6"
# New tests are "Gd_BOC_finest", "Gd_BOC_fine1", "Gd_BOC_t1", "Gd_BOC_test2", "Gd_BOC_test3", "Gd_BOC_test4", "Gd_BOC_test5"

# 5) Selecting the burnup calculation options
# Solver : "RUNG" or "KAPS"
# Saturation : "NODI" or "DIRA" : use eq 3.32 or eq 3.33 (with dirac delta) to treat saturation
# val_exp : list of values to impose saturation for isotopes with lambda * (xtf - xti) >= val_exp
#
solver_option = "KAPS"
saturation_option = "NODI" #, "DIRA"
rates_extr = "EXTR" # "EXTR", "NOEX"
# Options for PCC scheme :
bos_rates_opt = "EXTR" # "EXTR", "NOEX"
eos_rates_opt = "NOEX" # "EXTR", "NOEX"

# 6) Selecting the energy deposition options
# Global energy deposition : "NOGL"=only energy release in fuel is used for normalization or "GLOB" = global energy release model, 

glob_opt = "GLOB" # "GLOB", "NOGL"

# 7) Select which D5 case to run
exec_D5_no_modif = False # True : run DRAGON5 calculations, False : skip DRAGON5 calculations
exec_D5_no_NG0 = True # True : run DRAGON5 calculations with NG0 depletion chain, False : skip DRAGON5 calculations
exec_PCC_case = True # True : run PCC calculations, False : skip PCC calculations


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
    name_compo = f"_CPO_{draglib_name}_{ssh_option}_{tracking_option}_{solver_option}_{saturation_option}_{glob_opt}_{rates_extr}_{burnup_points}"
    if tracking_option == "SALT":
        if exec_PCC_case:
            name_compo = f"_CPO_{draglib_name}_{ssh_option}_{tracking_option}_{solver_option}_{glob_opt}_pcc{bos_rates_opt}_{eos_rates_opt}_{burnup_points}"
            CPO = BU_PCC("COMPO", LIB, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, name_compo, solver_option, glob_opt, bos_rates_opt, eos_rates_opt)
        else:
            name_compo = f"_CPO_{draglib_name}_{ssh_option}_{corr}_{tracking_option}_{solver_option}_{saturation_option}_{glob_opt}_{rates_extr}_{burnup_points}"
            CPO = BU_C("COMPO", LIB, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, name_compo, ssh_option, solver_option, glob_opt, saturation_option, rates_extr)
    elif tracking_option == "SYBNXT": 
        name_compo = f"_CPO_{draglib_name}_{ssh_option}_{corr}_{tracking_option}_{solver_option}_{saturation_option}_{glob_opt}_{rates_extr}_{burnup_points}"
        CPO = BU_C_SYBNXT("COMPO", LIB, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, name_compo, ssh_option, solver_option, glob_opt, saturation_option)
    print(f"creating D5 case for {name_compo}")
    D5case = D5_case(CPO, draglib_name, burnup_points, ssh_option, "CORR", saturation_option, solver_option, tracked_nuclides, BU_lists, save_dir_D5)
    D5case.plot_keffs()
    for iso in tracked_nuclides:
        D5case.plot_Ni(iso)


if exec_D5_no_NG0:
    ### Begin calculations with modified DEPL structure ###
        
    # --- Call to DRAGON5 CLE-2000 procedures :
    # --- DRAGON5 microlib generation
    pyLIB_NG0 = MIX_NG0(ssh_option,corr) # Creation of the microlib, default D5 energy deposition mode
    #
    # names for exportation
    if corr == "CORR":
        corr_name = "C"
    elif corr == "NOCORR":
        corr_name = "N"
    print(f"State of the calculation NG0 : {draglib_name} {ssh_option} {saturation_option} {solver_option}")
    compo_name = f"_CPO_{draglib_name}_NG0_{ssh_option}_{tracking_option}_{burnup_points}_{solver_option}_{saturation_option}_{rates_extr}_{glob_opt}"
    # run DRAGON5 calculation with BU evolution
    if tracking_option == "SALT":
        if exec_PCC_case:
            compo_name = f"_CPO_{draglib_name}_NG0_{ssh_option}_{corr_name}_{tracking_option}_{burnup_points}_{solver_option}_pcc{bos_rates_opt}_{eos_rates_opt}_{glob_opt}"
            CPO_NG0 = BU_PCC("COMPO", pyLIB_NG0, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, compo_name, solver_option, glob_opt, bos_rates_opt, eos_rates_opt)
        else:
            compo_name = f"_CPO_{draglib_name}_NG0_{ssh_option}_{corr_name}_{tracking_option}_{burnup_points}_{solver_option}_{saturation_option}_{rates_extr}_{glob_opt}"
            CPO_NG0 = BU_C("COMPO", pyLIB_NG0, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, compo_name, ssh_option, solver_option, glob_opt, rates_extr, saturation_option)
    elif tracking_option == "SYBNXT": 
        compo_name = f"_CPO_{draglib_name}_NG0_{ssh_option}_{corr_name}_{tracking_option}_{burnup_points}_{solver_option}_{saturation_option}_{rates_extr}_{glob_opt}"
        CPO_NG0 = BU_C_SYBNXT("COMPO", pyLIB_NG0, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, compo_name, ssh_option, solver_option, glob_opt, saturation_option)
    print(f"creating D5 case for {compo_name}")
    D5case_NG0 = D5_case(CPO_NG0, draglib_name, burnup_points, ssh_option, "CORR", saturation_option, solver_option, tracked_nuclides, BU_lists, save_dir_D5)
    D5case_NG0.plot_keffs()
    for iso in tracked_nuclides:
        D5case_NG0.plot_Ni(iso)
    # --- Post-processing of DRAGON5 results
        
    #### END SCRIPT #### 
