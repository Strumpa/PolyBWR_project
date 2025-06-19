# PyGan script for burnup evolution of homogeneous fuel cell with evolving Gd157
# Cell : AT10 45Gd homogenized
# Constant FLUX normalization of rates over time stages in EVO: -- > keyword FLUX 1.0E12 
# Author: R. Guasch
# Date: 2025-01-29
# Hypothesis to test : Difference (D5-S2) in N_Gd157 is due to differente energy deposition models between DRAGON5 and SERPENT2
# 

import sys
import numpy as np
import os, shutil
import lifo
import lcm
import cle2000

# POST-PROCESSING
from postproc_cst_flx_evol import DRAGON_case as D5_case
from postproc_cst_flx_evol import Serpent2_case as S2_case
from postproc_cst_flx_evol import D5S2_comparisons as D5S2_comparisons

# --- StepLists
from getLists import *

# Cle-2000 PROCEDURES
from HOM_Gd157_cflx import *

### Execution parameters
draglib = "J311_295" # or "J311_295"
ssh_module = "USS"
ssh_method = "RSE"
correlation = "CORR" # "CORR", "noCORR"
iso_saturation_options = ["SAT", "NSAT"]
solution_algorithms = ["RUNG", "KAPS"]

burnup_points = 'BOC_fine_autop5'

tracked_nuclides = ["U235","U238","Pu239","Pu240","Pu241","Pu242", 
                    "Gd154","Gd155","Gd156","Gd157","Gd158","Gd160",
                    "Xe135","Sm149","Sm154","Eu155","Eu156","Eu157"]

# Recovering ListBU ListAUTOP ListCOMPO
[ListeBU,ListeAUTOP,ListeCOMPO]=getLists(burnup_points)
BU_lists = {"BU": ListeBU, "AUTOP": ListeAUTOP, "COMPO": ListeCOMPO}
#
# Create Steplist for BU - SELFSHIELDING - COMPO save 
StepList = lcm.new('LCM','burnup_steps')
StepList['ListBU']    = np.array(ListeBU, dtype='f')
StepList['ListAutop'] = np.array(ListeAUTOP, dtype='f')
StepList['ListCompo'] = np.array(ListeCOMPO, dtype='f')
StepList.close() # close without erasing
pyCOMPOS_dict = {}

# Creation of results directory
path=os.getcwd()

a=os.path.exists(f"HOM_Gd157_Cst_flx_evol")
if a==False:
	os.mkdir("HOM_Gd157_Cst_flx_evol")

save_dir_D5 = f"{path}/HOM_Gd157_Cst_flx_evol/D5"
save_dir_S2 = f"{path}/HOM_Gd157_Cst_flx_evol/S2"
save_dir_comparison = f"{path}/HOM_Gd157_Cst_flx_evol/Comparison"
if not os.path.exists(save_dir_D5):
    os.makedirs(save_dir_D5)
if not os.path.exists(save_dir_S2):
    os.makedirs(save_dir_S2)
if not os.path.exists(save_dir_comparison):
    os.makedirs(save_dir_comparison)

list_of_D5_cases = []
ipLifo2 = lifo.new()
UOX_Gd157_Cst_flx_evol = cle2000.new('UOX_Gd157',ipLifo2,1)
# Execution of procedures
for sat in iso_saturation_options:
    for depl_sol in solution_algorithms:
        
        suffixe = f"cst_flx_{sat}_{depl_sol}"
        COMPO = HOM_Gd157_cflx(f"COMPO_Gd157_{suffixe}",StepList,f"./_COMPO_HOM_UOX_Gd157_{suffixe}_{ssh_module}_{ssh_method}",ssh_module,ssh_method,sat,depl_sol,draglib)
        # Post process constant flux evolution cases
        D5_test_case = D5_case(COMPO, "J311_295", burnup_points, ssh_module, ssh_method, correlation, sat, depl_sol, tracked_nuclides, BU_lists, save_dir_D5)
        D5_test_case.plot_keffs()
        for iso in tracked_nuclides:
            D5_test_case.plot_Ni(iso)
        
        list_of_D5_cases.append(D5_test_case)


        
S2_ref_case = S2_case("HOM_UOX_Gd157", "PyNjoy2016", 1, tracked_nuclides, save_dir_S2)
S2_ref_case.plot_keff()
S2_ref_case.plot_concentrations(["Gd157"])

# Compare D5 cases for different depletion and sturation options with S2 reference case
Compare = D5S2_comparisons(f"HOM_UOX_Gd157_{draglib}_cst_flx", list_of_D5_cases, S2_ref_case, tracked_nuclides, save_dir_comparison)
Compare.plot_delta_Keff()
Compare.plot_delta_Ni()