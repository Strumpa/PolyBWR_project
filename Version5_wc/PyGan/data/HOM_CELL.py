
##################################################################################
#                                                                                #
# Description : PyGan script for BWR simulation with DRAGON5                     #
# Author      : R. Guasch                                                        #
# Date        : 2024                                                             #
# Purpose     : Study isotopic evolution in the presence of UOX+Gd               #
# Case : Homogenized AT10_45Gd : C7 cell                                         #
#                                                                                #
##################################################################################
#
########################################################################################################################################################################################
#                                #
#             IMPORTS            #
#                                #
##################################
import sys
import numpy as np
import os, shutil
import lifo
import lcm
import cle2000
import matplotlib

# --- DRAGON5 procedures
from assertS import *
# POST-PROCESSING class
from POSTPROC_hom import *
from MULTI_SERP_POSTPROC import *

# --- HOMOGENEOUS CELL PROCEDURES

from HOM_U5 import *
from HOM_U5_U8 import *
from HOM_UOX import *
from UOXCladnoZr import *
from UOX_Clad import *
from UOX_Gd155 import *
from UOX_Gd157 import *
from UOX_no155157 import *
from C7_hom import *
from CHAIN_Gd157 import *

# --- OTHERS
from getLists import *



########################################################################################################################################################################################
#                                                        #
#   PARAMETERS OF FLUX CALCULATION : TO BE MODIFIED      #
#                                                        #
##########################################################
#
# Library = name of the Library used to build LIB: object
Library = 'J311_295'


name_geom = 'HOM_CELL'
cases_to_run = ["HOM_UOX_Gd157"] # "HOM_U5", "HOM_U5_U8", "HOM_UOX", "HOM_UOX_clad_noZr", "HOM_UOX_clad", "HOM_UOX_Gd155", "HOM_UOX_Gd157", "HOM_UOX_no155157", "HOM_UOXGd"


iso_chain_tests = [] #["NO_NG_toGd158","NO_NG_NO_ENE_toGd158"] # "NO_NP_toGd158", "NO_NG_toGd158"
iso_chain = "Nominal" # "Nominal", "NO_NP_toGd158", "NO_NG_toGd158", "NO_NG_NO_ENE_toGd158"

S2_libs = ["PyNjoy2016"]#,"oldlib"]
# Serpent2 energy deposition mode : set edep 0 = "Constant energy deposition per fission" --> "PyNjoy2016" and "oldlib" / default setting
# 											 1 = "Local energy deposition based on ENDF MT458 data"  --> "PyNjoy2016" and "oldlib" 
# 											 2 = "Local photon energy deposition" --> "PyNjoy2016" only
set_edep_mode_to_treat = {"PyNjoy2016":[0,1,2],"oldlib":[0]}

ssh_module = "AUTO" #"USS", "AUTO"
ssh_method = "SUBG" #"PT", "RSE", "SUBG" all supported for USS: but AUTO: only takes SUBG
# for RSE method : test with eps_RSE = 1.0E-1, 5.0E-2, 1.0E-2, 5.0E-3, 1.0E-3, 1.0E-4
correlation = "noCORR" # "CORR", "noCORR"

sat = "" # "NSAT", "SAT", "" for short-lived nuclide saturation or not
depl_sol = "KAPS" # "RUNG" : 5th order Runge-Kutta (Cash-Karp) or "KAPS": 4th order Kaps-Rentrop

#
# burnup_points = string to use predefine burnup points matching with SERPENT2 results

burnup_points = 'BOC_fine_autop5' #'Gd_VBOC_fine2_autop9' #'Gd_BOC_fine', 'Gd_autop4', 'Gd_autop3', 'Gd_VBOC_fine', 'Gd_BOC_fine'
# suffixe = suffixe added to name_geom for creation of figures, MULTICOMPO and BU vector
suffixe = burnup_points+"_"+iso_chain
#
Nmin = 0
#
# Visualisation options for post-processing functions
# 	= 1 if you want to plot and save the results
# 	= 0 either
# NB : if visu_XXX_Keff==1 or visu_XXX_BU==1 or visu_XXX_ISOTOPOESDENS==1 then visu_XXX must be equal to 1
visu_DRAGON=1
visu_SERPENT=1
visu_COMP=0
visu_DELTA=1
#
########################################################################################################################################################################################
#                                                 #
#         LOCAL PARAMETERS : DO NOT MODIFY        #
#                                                 #
###################################################

case =  'HOM_CELL' 

#
# names for exportation
if sat == "SAT":
	SAT = "_SAT"
elif sat == "SATOFF":
	SAT = "_SATOFF"
elif sat == "NSAT":
	SAT = "_NSAT"
else:
	SAT = ""
if correlation == "CORR":
	correlation_name = "_CORR"
elif correlation == "noCORR":
	correlation_name = "_noCORR"
else:
	correlation_name = ""

name_BUvector = f"{name_geom}_{suffixe}"

#
# Visualisation parameters stored in VISU_param
VISU_param=[visu_DRAGON,visu_SERPENT,visu_COMP,visu_DELTA]
#
# Creation of results directory
path=os.getcwd()

a=os.path.exists(f"BWRresults_PyGan_{name_geom}")
if a==False:
	os.mkdir(f"BWRresults_PyGan_{name_geom}")

SAVE_DIR=f'BWRresults_PyGan_{name_geom}/{suffixe}_postprocess/{ssh_module}_{ssh_method}/'
a=os.path.exists(SAVE_DIR)
if a==False:
	os.makedirs(SAVE_DIR)

#
# Recovering ListBU ListAUTOP ListCOMPO
[ListeBU,ListeAUTOP,ListeCOMPO]=getLists(burnup_points)
#
# Create Steplist for BU - SELFSHIELDING - COMPO save 
StepList = lcm.new('LCM','burnup_steps')
StepList['ListBU']    = np.array(ListeBU, dtype='f')
StepList['ListAutop'] = np.array(ListeAUTOP, dtype='f')
StepList['ListCompo'] = np.array(ListeCOMPO, dtype='f')
StepList.close() # close without erasing
#
# Save BU vector
os.chdir(path+'/'+'BWRresults_PyGan_'+name_geom)
np.savetxt(name_BUvector+'_BUvector.txt',ListeCOMPO)
os.chdir(path)

########################################################################################################################################################################################
#                                                 #
#             EXECUTION OF PROCEDURES             #
#                                                 #
###################################################

#
# --------------------------------------
#             HOMOGENEOUS CELL                 
# --------------------------------------

if ssh_module == "AUTO":
	ssh_method = "SUBG"
	correlation = "noCORR"
pyCOMPOs = {}
print(cases_to_run)
if "HOM_U5" in cases_to_run:
	# -- Test Case 1 :
	pyCOMPO_HOM_U5 = HOM_U5("COMPO_U5",StepList,f"./_COMPO_HOM_U5_{suffixe}_{depl_sol}{SAT}_{ssh_module}_{ssh_method}",ssh_module,ssh_method,sat,depl_sol)
	pyCOMPOs["HOM_U5"] = pyCOMPO_HOM_U5

if "HOM_U5_U8" in cases_to_run:
	# -- Test Case 2 :
	pyCOMPO_HOM_U5_U8 = HOM_U5_U8("COMPO_U5U8",StepList,f"./_COMPO_HOM_U5_U8_{suffixe}_{depl_sol}{SAT}_{ssh_module}_{ssh_method}{correlation_name}",ssh_module,ssh_method,correlation,sat,depl_sol)
	pyCOMPOs["HOM_U5_U8"] = pyCOMPO_HOM_U5_U8

if "HOM_UOX" in cases_to_run:
	# -- Test Case 3 :
	pyCOMPO_HOM_UOX = HOM_UOX("COMPO_UOX",StepList,f"./_COMPO_HOM_UOX_{suffixe}_{depl_sol}{SAT}_{ssh_module}_{ssh_method}{correlation_name}",ssh_module,ssh_method,correlation,sat,depl_sol)
	pyCOMPOs["HOM_UOX"] = pyCOMPO_HOM_UOX

if "HOM_UOX_clad_noZr" in cases_to_run:
	# -- Test Case 4 :
	pyCOMPO_HOM_UOX_clad_noZr = UOXCladnoZr("COMPO_CnZ",StepList,f"./_COMPO_HOM_UOX_clad_noZr_{suffixe}_{depl_sol}{SAT}_{ssh_module}_{ssh_method}",ssh_module,ssh_method,sat,depl_sol)
	pyCOMPOs["HOM_UOX_clad_noZr"] = pyCOMPO_HOM_UOX_clad_noZr

if "HOM_UOX_clad" in cases_to_run:
	# -- Test Case 5 :
	pyCOMPO_HOM_UOX_clad = UOX_Clad("COMPO_C",StepList,f"./_COMPO_HOM_UOX_clad_{suffixe}_{depl_sol}{SAT}_{ssh_module}_{ssh_method}",ssh_module,ssh_method,sat,depl_sol)
	pyCOMPOs["HOM_UOX_clad"] = pyCOMPO_HOM_UOX_clad

if "HOM_UOX_Gd155" in cases_to_run:
	# -- Test Case 6 :
	pyCOMPO_HOM_UOX_Gd155 = UOX_Gd155("COMPO_Gd155",StepList,f"./_COMPO_HOM_UOX_Gd155_{suffixe}_{depl_sol}{SAT}_{ssh_module}_{ssh_method}",ssh_module,ssh_method,sat,depl_sol)
	pyCOMPOs["HOM_UOX_Gd155"] = pyCOMPO_HOM_UOX_Gd155

if "HOM_UOX_Gd157" in cases_to_run:
	# -- Test Case 7 :
	print("In test case 7")
	if "NO_NP_toGd158" in iso_chain_tests:
		pyCOMPO_HOM_UOX_Gd157_test_noNP = CHAIN_Gd157("COMPO_Gd157_test_noNP",StepList,f"./_COMPO_HOM_UOX_Gd157_NONP_{suffixe}_{depl_sol}{SAT}",ssh_module,sat,depl_sol,chain_modif="NO_NP_toGd158")
		pyCOMPOs["HOM_UOX_Gd157_NO_NP_toGd158"] = pyCOMPO_HOM_UOX_Gd157_test_noNP
	if "NO_NG_toGd158" in iso_chain_tests:
		pyCOMPO_HOM_UOX_Gd157_test_noNG = CHAIN_Gd157("COMPO_Gd157_test_noNG",StepList,f"./_COMPO_HOM_UOX_Gd157_NONG_{suffixe}_{depl_sol}{SAT}",ssh_module,sat,depl_sol,chain_modif="NO_NG_toGd158")
		pyCOMPOs["HOM_UOX_Gd157_NO_NG_toGd158"] = pyCOMPO_HOM_UOX_Gd157_test_noNG
	if "NO_NG_NO_ENE_toGd158" in iso_chain_tests:
		pyCOMPO_HOM_UOX_Gd157_test_noNG_noENE = CHAIN_Gd157("COMPO_Gd157_test_noNG_noENE",StepList,f"./_COMPO_HOM_UOX_Gd157_NONGNOENE_{suffixe}_{depl_sol}{SAT}",ssh_module,sat,depl_sol,chain_modif="NO_NG_NO_ENE_toGd158")
		pyCOMPOs["HOM_UOX_Gd157_NO_NG_NO_ENE_toGd158"] = pyCOMPO_HOM_UOX_Gd157_test_noNG_noENE
	
	pyCOMPO_HOM_UOX_Gd157 = UOX_Gd157("COMPO_Gd157",StepList,f"./_COMPO_HOM_UOX_Gd157_{suffixe}_{depl_sol}{SAT}_{ssh_module}_{ssh_method}",ssh_module,ssh_method,sat,depl_sol)
	pyCOMPOs["HOM_UOX_Gd157"] = pyCOMPO_HOM_UOX_Gd157

if "HOM_UOX_no155157" in cases_to_run:
	# -- Test Case 8 :
	pyCOMPO_HOM_UOX_no155157 = UOX_no155157("COMPO_no5557",StepList,f"./_COMPO_HOM_UOX_no155157_{suffixe}_{depl_sol}{SAT}_{ssh_module}_{ssh_method}",ssh_module,ssh_method,sat,depl_sol)
	pyCOMPOs["HOM_UOXGd_no155157"] = pyCOMPO_HOM_UOX_no155157

if "HOM_UOXGd" in cases_to_run:
	# -- Test Case 9 :
	pyCOMPO_HOM_UOXGd = C7_hom("COMPO_UOXGd", StepList, f"./_COMPO_HOM_UOXGd_{suffixe}_{depl_sol}{SAT}_{ssh_module}_{ssh_method}", ssh_module, ssh_method, "Nominal", sat, depl_sol)
	pyCOMPOs["HOM_UOXGd"] = pyCOMPO_HOM_UOXGd
# --------------------------------------
# 		   POST-PROCESSING

print(pyCOMPOs)
#POSTPROC_hom(pyCOMPO, ListeCOMPO, ListeAUTOP, name_geom, name_mix, suffixe, VISU_param, Nmin, GdCompo, S2_libs, ssh_module, ssh_method, correlation, sat, depl_sol)
MULTI_SERP_POSTPROC(pyCOMPOs, ListeCOMPO, ListeAUTOP, name_geom, suffixe, VISU_param, Nmin, S2_libs, ssh_module, ssh_method, correlation, depl_sol, sat, set_edep_mode_to_treat)