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
draglibs_to_test = ["endfb8r1_295", "endfb81295K", "endfb81295K2"] # "endfb8r1_295", "endfb81295K" : with total KERMA (MT301), "endfb81295K2" : with modified KERMA (MT301-MT318)+MT458 data.

# 3) Selecting the self-shielding method
# RSE, PT or AUTO
#
ssh_options = ["PT"] # "RSE" , "PT", "SUBG", "AUTO"
corr_options = ["NOCORR"] # "CORR" or "NOCORR"

# 4) Selecting the burnup calculation options
# burnup_steps = "UOx", "UOx_autop5", "UOx2_autop5", "UOx4_autop5", "UOx6_autop5" etc
#
burnup_points_to_test = ["Gd_autop3", "Gd2_autop6"] #"Gd_autop3" # "Gd_autop4", "Gd_autop5", "Gd", "Gd2_autop6"
# New tests are "Gd_BOC_finest", "Gd_BOC_fine1", "Gd_BOC_t1", "Gd_BOC_test2", "Gd_BOC_test3", "Gd_BOC_test4", "Gd_BOC_test5"

# 5) Selecting the burnup calculation options
# Solver : "RUNG" or "KAPS"
# Saturation : "NODI" or "DIRA" : use eq 3.32 or eq 3.33 (with dirac delta) to treat saturation
# val_exp : list of values to impose saturation for isotopes with lambda * (xtf - xti) >= val_exp
#
solver_option = "KAPS"
saturation_option = "NODI" #, "DIRA"

# 6) Selecting the energy deposition options
# Global energy deposition : "NOGL"=only energy release in fuel is used for normalization or "GLOB" = global energy release model, 

glob_opt = "NOGL" # "GLOB", "NOGL"

# 7) Select which D5 case to run
exec_D5_no_modif = False # True : run DRAGON5 calculations, False : skip DRAGON5 calculations
exec_D5_no_NG0 = True # True : run DRAGON5 calculations with NG0 depletion chain, False : skip DRAGON5 calculations

#time_integrator = "CECM" # "CECM", "CECE", "EXTR", "NOEX", "EXTR2"
time_integrator_list = ["EXTR"]
#time_integrator_list = ["CECE", "EXTR", "NOEX", "EXTR2", "CECM"]

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


### BEGIN DRAGON5 calculations with default DEPL structure ###

if tracking_option == "SALT":
    # geometry definition
    GEOM, GEOM_SS = GEO_C_SALT()
    # tracking of geometries
    TRK, TF_EXC, TRK_SS, TF_EXC_SS = TRK_C_SALT(GEOM,GEOM_SS) 
elif tracking_option == "SYBNXT":
    GEOM, GEOM_SS = GEO_C_NXT()
    TRK, TF_EXC, TRK_SS = TRK_C_SYBNXT(GEOM,GEOM_SS)
else:
    print("Tracking option not recognized")
    sys.exit(1)

if exec_D5_no_modif:
    for draglib_name in draglibs_to_test:
        for ssh_option in ssh_options:
            for corr in corr_options:
                if corr == "CORR":
                    corr_name = "C"
                elif corr == "NOCORR":
                    corr_name = "N"
                else:
                    print("Correlation option not recognized")
                    sys.exit(1)
                # Materials library
                LIB = MIX_C(draglib_name,ssh_option,corr)
                for burnup_points in burnup_points_to_test:
                    for time_integrator in time_integrator_list:
                        [BUList,SSHList,COMPOList]=getLists(burnup_points)
                        # Create Steplist for BU - SELFSHIELDING - COMPO save 
                        StepList = lcm.new('LCM','burnup_steps')
                        StepList['ListBU']    = np.array(BUList, dtype='f')
                        StepList['ListAutop'] = np.array(SSHList, dtype='f')
                        StepList['ListCompo'] = np.array(COMPOList, dtype='f')
                        StepList.close() # close without erasing
                        # names for exportation

                        name_compo = f"_CPO_{draglib_name}_{ssh_option}_{corr_name}_{tracking_option}_{burnup_points}_{solver_option}_{saturation_option}_{time_integrator}_{glob_opt}"
                        if time_integrator == "CECM":
                            CPO = CECM("COMPO", LIB, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, name_compo, solver_option, glob_opt)
                        elif time_integrator == "CECE":
                            CPO = CECE("COMPO", LIB, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, name_compo, solver_option, glob_opt)
                        elif time_integrator == "EXTR":
                            CPO = BU_C("COMPO", LIB, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, name_compo, ssh_option, solver_option, glob_opt, saturation_option, rates_extr="EXTR")
                        elif time_integrator == "NOEX":
                            CPO = BU_C("COMPO", LIB, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, name_compo, ssh_option, solver_option, glob_opt, saturation_option, rates_extr="NOEX")
                        elif time_integrator == "EXTR2":
                            CPO = BU_EXTR2("COMPO", LIB, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, name_compo, ssh_option, solver_option, glob_opt, saturation_option)

if exec_D5_no_NG0:
    ### Begin calculations with modified DEPL structure ###
    for draglib_name in draglibs_to_test:
        for ssh_option in ssh_options:
            for corr in corr_options:
                if corr == "CORR":
                    corr_name = "C"
                elif corr == "NOCORR":
                    corr_name = "N"
                else:
                    print("Correlation option not recognized")
                    sys.exit(1)
                # Materials library
                pyLIB_NG0 = MIX_NG0(ssh_option,corr)
                for burnup_points in burnup_points_to_test:
                    for time_integrator in time_integrator_list:
                        #
                        print(f"State of the calculation NG0 : {draglib_name} {ssh_option} {saturation_option} {solver_option}")
            
                        [BUList,SSHList,COMPOList]=getLists(burnup_points)
                        # Create Steplist for BU - SELFSHIELDING - COMPO save 
                        StepList = lcm.new('LCM','burnup_steps')
                        print(f"StepList : {BUList}")
                        StepList['ListBU']    = np.array(BUList, dtype='f')
                        StepList['ListAutop'] = np.array(SSHList, dtype='f')
                        StepList['ListCompo'] = np.array(COMPOList, dtype='f')
                        StepList.close() # close without erasing
                        compo_name = f"_CPO_{draglib_name}_NG0_{ssh_option}_{corr_name}_{tracking_option}_{burnup_points}_{solver_option}_{saturation_option}_{time_integrator}_{glob_opt}"
                        # run DRAGON5 calculation with BU evolution
                        if time_integrator == "CECM":
                            CPO_NG0 = CECM("COMPO", pyLIB_NG0, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, compo_name, solver_option, glob_opt)
                        elif time_integrator == "CECE":
                            CPO_NG0 = CECE("COMPO", pyLIB_NG0, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, compo_name, solver_option, glob_opt)
                        elif time_integrator == "EXTR":
                            rates_extr = "EXTR"
                            CPO_NG0 = BU_C("COMPO", pyLIB_NG0, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, compo_name, ssh_option, solver_option, glob_opt, rates_extr, saturation_option)
                        elif time_integrator == "NOEX":
                            rates_extr = "NOEX"
                            CPO_NG0 = BU_C("COMPO", pyLIB_NG0, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, compo_name, ssh_option, solver_option, glob_opt, rates_extr, saturation_option)
                        elif time_integrator == "EXTR2":
                            CPO_NG0 = BU_EXTR2("COMPO", pyLIB_NG0, TRK, TF_EXC, TRK_SS, TF_EXC_SS, StepList, compo_name, ssh_option, solver_option, glob_opt, saturation_option)    
#### END SCRIPT #### 
