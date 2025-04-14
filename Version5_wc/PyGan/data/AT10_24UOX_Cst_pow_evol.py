### PyGan script for validation study of AT10_24UOX single pincell benchamrk
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
from MIX_E0 import *
# TRACKING class
from TRK_C_SALT import *
from TRK_C_SYBNXT import *
# BU DEPLETION class
from BU_C import *
from BU_C_SYBNXT import *

# --- OTHERS
from getLists import *

# --- DRAGOM5-SERPENT2 comparison : post treatment
from postproc_cst_pow_evol import DRAGON_case as D5_case
from postproc_cst_pow_evol import Serpent2_case as S2_case
from postproc_cst_pow_evol import multiD5S2_comparisons as multiD5S2
from postproc_cst_pow_evol import D5multiS2_comparisons as D5multiS2 


# --- Selecting calculation options
# 0) Select if executing stupid tests with ENDFb8r1 DEPL structure manually modified --> not compatible with Jeff311 lib 
run_NG0 = False
run_E0 = False 

# 1) Selecting the geometry tracking options
# Tracking : SALT or SYBNXT
# 
tracking_options = ["SALT"] #, "SYBNXT"]

# 2) Selecting the evaluation used for the simulation
# ENDFb8r1_295 or J311_295
#
evaluation = "endfb8r1" # "endfb8r1" # "Jeff311"
draglibs = ["endfb8r1_295", "endfb81295K", "endfb81295K2"]  # ["endfb8r1_295", "endfb81295K", "endfb81295K2"] # ["J311_295", "J311_295K"] #, "ENDFb8r1_295"]

# 3) Selecting the self-shielding method
# RSE, PT or AUTO
#
ssh_options = ["PT"]

# 4) Selecting the burnup calculation options
# burnup_steps = "UOx", "UOx_autop5", "UOx2_autop5", "UOx4_autop5", "UOx6_autop5" etc
#
burnup_steps_to_test = ["UOx2_autop5"]

# ["UOx_autop1", "UOx_autop2", "UOx_autop4", "UOx_autop5"]
# ["UOx2_autop1", "UOx2_autop3", "UOx2_autop4", "UOx2_autop5", "UOx2_autop6"]
# ["UOx4_autop1", "UOx4_autop3", "UOx4_autop4", "UOx4_autop5"]
# ["UOx6_autop1", "UOx6_autop3", "UOx6_autop4", "UOx6_autop6", "UOx6_autop7", "UOx6_autop8"]
# ["UOx8_autop8", "UOx8_autop9"]

# 5) Selecting the burnup calculation options
# Solver : "RUNG" or "KAPS"
# Saturation : "NODI" or "DIRA" : use eq 3.32 or eq 3.33 (with dirac delta) to treat saturation
# val_exp : list of values to impose saturation for isotopes with lambda * (xtf - xti) >= val_exp
#
solver_options = ["RUNG"]
saturation_options = ["NODI"]

# 6) Selecting the energy deposition options
# Global energy deposition : "NOGL"=only energy release in fuel is used for normalization or "GLOB" = global energy release model, 

glob_opt = "GLOB" # ["GLOB"] #, "NOGL"]

######## Options for DRAGON5-SERPENT2 comparison ########
# Select origin of SERPENT2 data + evaluation
#
S2_evaluations = ["endfb8r1","J311"]
Njoy_versions = ["pynjoy2012","NJOY2016", "PyNjoy2016"]

######## 
# Result handling and creation of the results directory
#
tracked_nuclides = ["U235","U238","Pu239","Pu240","Pu241","Pu242","Am241","Xe135","Sm149","Gd155","Gd157"]

#
# Create the results directory
path=os.getcwd()
save_dir_D5 = f"{path}/AT10_24UOX_Cst_pow_evol_results/D5_{evaluation}"
save_dir_S2 = f"{path}/AT10_24UOX_Cst_pow_evol_results/S2_{evaluation}"
save_dir_comparison = f"{path}/AT10_24UOX_Cst_pow_evol_results/Comparison_{evaluation}"
if not os.path.exists(save_dir_D5):
    os.makedirs(save_dir_D5)
if not os.path.exists(save_dir_S2):
    os.makedirs(save_dir_S2)
if not os.path.exists(save_dir_comparison):
    os.makedirs(save_dir_comparison)




### BEGIN DRAGON5 calculations ###
D5_cases = []

for trk_opt in tracking_options:
    if trk_opt == "SALT":
        # geometry definition
        GEOM, GEOM_SS = GEO_C_SALT()
        # tracking of geometries
        TRK, TF_EXC, TRK_SS, TF_EXC_SS = TRK_C_SALT(GEOM,GEOM_SS) # CP = Collision Probability, MOC = Method of Characteristics
    elif trk_opt == "SYBNXT":
        GEOM, GEOM_SS = GEO_C_NXT()
        TRK, TF_EXC, TRK_SS = TRK_C_SYBNXT(GEOM,GEOM_SS)
    else:
        print("Tracking option not recognized")
        sys.exit(1)
    for dlib_name in draglibs:
        for ssh_opt in ssh_options:
            LIB = MIX_C(dlib_name,ssh_opt)
            for burnup_points in burnup_steps_to_test:
                # Recovering ListBU ListAUTOP ListCOMPO
                [ListeBU,ListeAUTOP,ListeCOMPO]=getLists(burnup_points)
                BU_lists = {"BU": ListeBU, "AUTOP": ListeAUTOP, "COMPO": ListeCOMPO}
                # Create Steplist for BU - SELFSHIELDING - COMPO save 
                StepList = lcm.new('LCM','burnup_steps')
                StepList['ListBU']    = np.array(ListeBU, dtype='f')
                StepList['ListAutop'] = np.array(ListeAUTOP, dtype='f')
                StepList['ListCompo'] = np.array(ListeCOMPO, dtype='f')
                StepList.close() # close without erasing

                for solver_opt in solver_options:
                    for sat_opt in saturation_options: 
                        name_compo = f"_CPO_AT10_24UOX_{dlib_name}_{ssh_opt}_{trk_opt}_{solver_opt}_{sat_opt}_{glob_opt}_{burnup_points}"
                        if trk_opt == "SALT":
                            CPO = BU_C("COMPO", LIB, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, name_compo, ssh_opt, solver_opt, glob_opt, sat_opt)
                        elif trk_opt == "SYBNXT": 
                            CPO = BU_C_SYBNXT("COMPO", LIB, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, name_compo, ssh_opt, solver_opt, glob_opt, sat_opt)
                        print(f"DRAGON5 calculation for {name_compo}")
                        D5case = D5_case(CPO, dlib_name, burnup_points, ssh_opt, "noCORR", sat_opt, solver_opt, tracked_nuclides, BU_lists, save_dir_D5)
                        D5case.plot_keffs()
                        for iso in tracked_nuclides:
                            D5case.plot_Ni(iso)
                        print(f"appendind D5case for dl {dlib_name}, burnup {burnup_points}, ssh {ssh_opt}, trk {trk_opt}, solver {solver_opt}, sat {sat_opt}")
                        D5_cases.append(D5case)


if run_NG0:
    D5_NG0_cases = []

    ### BEGIN DRAGON5 calculations for NG Q-values = 0 ###
    D5_NG0_cases = []
    for dlib_name in draglibs:
        for ssh_opt in ssh_options:
            LIB_NG0 = MIX_NG0(dlib_name) ### Only modified for ENDFb8r1 chain!!!!! ---> fix this or Im gonna jump
            for burnup_points in burnup_steps_to_test:
                # Recovering ListBU ListAUTOP ListCOMPO
                [ListeBU,ListeAUTOP,ListeCOMPO]=getLists(burnup_points)
                BU_lists = {"BU": ListeBU, "AUTOP": ListeAUTOP, "COMPO": ListeCOMPO}
                # Create Steplist for BU - SELFSHIELDING - COMPO save 
                StepList = lcm.new('LCM','burnup_steps')
                StepList['ListBU']    = np.array(ListeBU, dtype='f')
                StepList['ListAutop'] = np.array(ListeAUTOP, dtype='f')
                StepList['ListCompo'] = np.array(ListeCOMPO, dtype='f')
                StepList.close() # close without erasing

                for solver_opt in solver_options:
                    for sat_opt in saturation_options:
                        name_compo_NG0 = f"_CPO_AT10_24UOX_{dlib_name}NG0_{ssh_opt}_{trk_opt}_{solver_opt}_{sat_opt}_{glob_opt}_{burnup_points}"
                        if trk_opt == "SALT":
                            CPO_NG0 = BU_C("COMPO_NG0", LIB_NG0, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, name_compo_NG0, ssh_opt, solver_opt, glob_opt, sat_opt)
                        elif trk_opt == "SYBNXT": 
                            CPO_NG0 = BU_C_SYBNXT("COMPO_NG0", LIB_NG0, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, name_compo_NG0, ssh_opt, solver_opt, glob_opt, sat_opt)
                        D5NG0case = D5_case(CPO_NG0, f"{dlib_name}NG0", burnup_points, ssh_opt, "noCORR", sat_opt, solver_opt, tracked_nuclides, BU_lists, save_dir_D5)
                        D5NG0case.plot_keffs()
                        for iso in tracked_nuclides:
                            D5NG0case.plot_Ni(iso)
                        D5_NG0_cases.append(D5NG0case)


# --- END OF SCRIPT --- #

