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

# 1) Selecting the geometry tracking options
# Tracking : SALT or SYBNXT
# 
tracking_options = ["SALT"] #, "SYBNXT"]

# 2) Selecting the evaluation used for the simulation
# ENDFb8r1_295 or J311_295
#
draglibs = ["ENDFb8r1_295"]

# 3) Selecting the self-shielding method
# RSE, PT or AUTO
#
ssh_options = ["RSE"]

# 4) Selecting the burnup calculation options
# burnup_steps = "UOx", "UOx_autop5", "UOx2_autop5", "UOx4_autop5", "UOx6_autop5" etc
#
burnup_steps_to_test = ["UOx2_autop1", "UOx2_autop3", "UOx2_autop4", "UOx2_autop5", "UOx2_autop6"]

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

ene_deposition_options = ["GLOB"]

######## Options for DRAGON5-SERPENT2 comparison ########
# Select origin of SERPENT2 data + evaluation
#
S2_evaluations = ["endfb8r1"]#,"Jef3.1.1"]
Njoy_versions = ["pynjoy2012","NJOY2016"], #"PyNjoy2016"

######## 
# Result handling and creation of the results directory
#
tracked_nuclides = ["U235","U238","Pu239","Pu240","Pu241","Pu242","Am241","Xe135","Sm149","Gd155","Gd157"]

#
# Create the results directory
path=os.getcwd()
save_dir_D5 = f"{path}/AT10_24UOX_Cst_pow_evol_results/D5_PASS4_GRMIN35"
save_dir_S2 = f"{path}/AT10_24UOX_Cst_pow_evol_results/S2"
save_dir_comparison = f"{path}/AT10_24UOX_Cst_pow_evol_results/Comparison_USS4PASS_GRMIN35"
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
                        for glob_opt in ene_deposition_options: 
                            name_compo = f"COMPO_AT10_24UOX_{dlib_name}_{ssh_opt}_{trk_opt}_{solver_opt}_{sat_opt}_{glob_opt}_{burnup_points}"
                            if trk_opt == "SALT":
                                CPO = BU_C("COMPO", LIB, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, name_compo, ssh_opt, solver_opt, glob_opt, sat_opt, val_exp = 80.0)
                            elif trk_opt == "SYBNXT": 
                                CPO = BU_C_SYBNXT("COMPO", LIB, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, name_compo, ssh_opt, solver_opt, glob_opt, sat_opt, val_exp = 80.0)
                            D5case = D5_case(CPO, dlib_name, burnup_points, ssh_opt, "noCORR", sat_opt, solver_opt, tracked_nuclides, BU_lists, save_dir_D5)
                            D5case.plot_keffs()
                            for iso in tracked_nuclides:
                                D5case.plot_Ni(iso)
                            D5_cases.append(D5case)

### BEGIN SERPENT2 post treatment
                            


S2_endfb8r1_edep0_pcc1 = S2_case("AT10_24UOX", "endfb8r1_pynjoy2012", 0, False, 1, 38.6, tracked_nuclides, save_dir_S2)
S2_endfb8r1_edep0_pcc1.plot_keff()
for iso in tracked_nuclides:
    S2_endfb8r1_edep0_pcc1.plot_concentrations([iso])

S2_endfb8r1_edep0_pcc2 = S2_case("AT10_24UOX", "endfb8r1_pynjoy2012", 0, False, 2, 38.6, tracked_nuclides, save_dir_S2)
S2_endfb8r1_edep0_pcc2.plot_keff()
for iso in tracked_nuclides:
    S2_endfb8r1_edep0_pcc2.plot_concentrations([iso])
"""
S2_endfb8r1_edep1 = S2_case("AT10_24UOX", "endfb8r1_pynjoy2012", 1, False, 1, 38.6, tracked_nuclides, save_dir_S2)
S2_endfb8r1_edep0_pcc1.plot_keff()
for iso in tracked_nuclides:
    S2_endfb8r1_edep0_pcc1.plot_concentrations([iso])
"""

### Compare DRAGON5 and SERPENT2 results

D5_cases_to_endfb8r1_edep0_pcc1 = multiD5S2(f"24UOX evolution {burnup_steps_to_test[0].split('_')[0]} {ssh_options[0]} vs Serpent2 edep0 pcc1", D5_cases, S2_endfb8r1_edep0_pcc1, tracked_nuclides, save_dir_comparison)
D5_cases_to_endfb8r1_edep0_pcc2 = multiD5S2(f"24UOX evolution {burnup_steps_to_test[0].split('_')[0]} {ssh_options[0]} vs Serpent2 edep0 pcc2", D5_cases, S2_endfb8r1_edep0_pcc2, tracked_nuclides, save_dir_comparison)

D5_cases_to_endfb8r1_edep0_pcc1.compare_keffs()
D5_cases_to_endfb8r1_edep0_pcc2.compare_keffs()
D5_cases_to_endfb8r1_edep0_pcc1.compare_Ni()
D5_cases_to_endfb8r1_edep0_pcc2.compare_Ni()
D5_cases_to_endfb8r1_edep0_pcc1.plot_delta_Keff()
D5_cases_to_endfb8r1_edep0_pcc2.plot_delta_Keff()
D5_cases_to_endfb8r1_edep0_pcc1.plot_delta_Ni()
D5_cases_to_endfb8r1_edep0_pcc2.plot_delta_Ni()


