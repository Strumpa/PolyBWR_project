##################################################################################
#                                                                                #
# Description : PyGan script for BWR simulation with DRAGON5                     #
# Author      : R. Guasch                                                        #
# Date        : 2024                                                             #
# Purpose     : Study isotopic evolution in the presence of UOX+Gd               #
#                and compare eulerian vs predictor-corrector methods             #
# Case : Homogenized AT10 Cell with Gd157                                        #
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


# POST-PROCESSING 
from POSTPROC_Gd157 import *

# --- StepLists 
from getLists import *

# Cle-2000 PROCEDURES
from MixHGd157 import * # Builds the LIBRARY object
from PCC0_Gd157 import * # Performs eulerian step predictor BU method with Constant Extrapolation only
from PCC1_Gd157 import * # Performs predicor-corrector BU method with Constant Extrapolation and linear interpolation (?)
from PCC2_Gd157 import * # Performs predicor step with Linear Extrapolation
from PCC3_Gd157 import * # Performs predicor-corrector Linear Extrapolation and corrective step
from PCC3b_Gd157 import * # Same as PCC3 but with a second EVO: run to ensure convergence


########################################################################################################################################################################################
#                                                   			     #
#   PARAMETERS OF FLUX/EVOLUTION CALCULATION : TO BE MODIFIED        #
#                                                        			 #
######################################################################
#
# Library = name of the Library used to build LIB: object
Library = 'J311_315' # J311_281, J311_295, J311_315, J311_172
# If Library = J311_172, then a PyNjoy2016 generated library is used
name_study = 'HOM_UOX_Gd157_PCC' # name of the study, for exportation purposes only

S2_libs = ["PyNjoy2016","oldlib"]
S2_PCC = ["PCC0","PCC1","PCC2","PCC3","PCC4","PCC6"] # all ran for PyNjoy2016, waiting on 3, 4 and 6 for oldlib

ssh_module = "AUTO" #"USS", "AUTO", "SHI"
ssh_method = "SUBG" #"PT", "RSE", "SUBG" all supported for USS: but AUTO: only takes SUBG
if ssh_module == "AUTO":
	ssh_method = "SUBG"
# for RSE method : test with eps_RSE = 1.0E-1, 5.0E-2, 1.0E-2, 5.0E-3, 1.0E-3, 1.0E-4
correlation = "noCORR" # "CORR", "noCORR"

# burnup_points = string to use predefine burnup points matching with SERPENT2 results

burnup_points = 'BOC_fine_autop5' #'Gd_VBOC_fine2_autop9' #'Gd_BOC_fine', 'Gd_autop4', 'Gd_autop3', 'Gd_VBOC_fine', 'Gd_BOC_fine'
# suffixe = suffixe added to name_geom for creation of figures, MULTICOMPO and BU vector
suffixe = burnup_points
#
Nmin = 0
#
# Visualisation options for post-processing functions
# 	= 1 if you want to plot and save the results
# 	= 0 either
# NB : if visu_XXX_Keff==1 or visu_XXX_BU==1 or visu_XXX_ISOTOPOESDENS==1 then visu_XXX must be equal to 1
visu_DRAGON=0
visu_SERPENT=0
visu_COMP=1
visu_DELTA=1

run_DRAGON = True
post_process_Gd157 = False
post_treat_S2 = False

########################################################################################################################################################################################
#                                                 #
#         LOCAL PARAMETERS : DO NOT MODIFY        #
#                                                 #
###################################################

if run_DRAGON:
	if correlation == "CORR":
		correlation_name = "_CORR"
	elif correlation == "noCORR":
		correlation_name = "_noCORR"
	else:
		correlation_name = ""

	name_BUvector = f"{name_study}_{suffixe}"

	#
	# Visualisation parameters stored in VISU_param
	VISU_param=[visu_DRAGON,visu_SERPENT,visu_COMP,visu_DELTA]
	#
	# Creation of results directory
	path=os.getcwd()

	a=os.path.exists(f"PyGan_results_{name_study}")
	print(f"PyGan_results_{name_study} exists: {a}")
	if a==False:
		os.mkdir(f"PyGan_results_{name_study}")

	SAVE_DIR=f'PyGan_results_{name_study}/{suffixe}_postprocess/{ssh_module}_{ssh_method}/'
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
	os.chdir(path+'/'+'PyGan_results_'+name_study)
	np.savetxt(name_BUvector+'_BUvector.txt',ListeCOMPO)
	os.chdir(path)


	########################################################################################################################################################################################
	#                                                 #
	#      EXECUTION OF CLE-2000 PROCEDURES           #         
	#                                                 #
	###################################################

	pyCOMPOs = {}
	pyMIX = MixHGd157("LIBRARY", Library, ssh_method, correlation)
	pyCOMPO_HOM_Gd157_PCC0 = PCC0_Gd157("COMPO_Gd157_PCC0", pyMIX, StepList, f"./_COMPO_HOM_UOX_Gd157_PCC0_{Library}_{suffixe}_{ssh_module}_{ssh_method}_{correlation}", ssh_module)
	pyCOMPOs["HOM_Gd157_PCC0"] = pyCOMPO_HOM_Gd157_PCC0

	pyCOMPO_HOM_Gd157_PCC1 = PCC1_Gd157("COMPO_Gd157_PCC1", pyMIX, StepList, f"./_COMPO_HOM_UOX_Gd157_PCC1_{Library}_{suffixe}_{ssh_module}_{ssh_method}_{correlation}", ssh_module)
	pyCOMPOs["HOM_Gd157_PCC1"] = pyCOMPO_HOM_Gd157_PCC1

	pyCOMPO_HOM_Gd157_PCC2 = PCC2_Gd157("COMPO_Gd157_PCC2", pyMIX, StepList, f"./_COMPO_HOM_UOX_Gd157_PCC2_{Library}_{suffixe}_{ssh_module}_{ssh_method}_{correlation}", ssh_module)
	pyCOMPOs["HOM_Gd157_PCC2"] = pyCOMPO_HOM_Gd157_PCC2

	pyCOMPO_HOM_Gd157_PCC3 = PCC3_Gd157("COMPO_Gd157_PCC3", pyMIX, StepList, f"./_COMPO_HOM_UOX_Gd157_PCC3_{Library}_{suffixe}_{ssh_module}_{ssh_method}_{correlation}", ssh_module)
	pyCOMPOs["HOM_Gd157_PCC3"] = pyCOMPO_HOM_Gd157_PCC3

	pyCOMPO_HOM_Gd157_PCC3b = PCC3b_Gd157("COMPO_Gd157_PCC3b", pyMIX, StepList, f"./_COMPO_HOM_UOX_Gd157_PCC3b_{Library}_{suffixe}_{ssh_module}_{ssh_method}_{correlation}", ssh_module)
	pyCOMPOs["HOM_Gd157_PCC3b"] = pyCOMPO_HOM_Gd157_PCC3b


	########################################################################################################################################################################################
	#                                                 #
	#             POST-PROCESSING                     #
	#                                                 #
	###################################################
if post_process_Gd157:
	print("Starting post-processing")
	print(f"pyCompos.keys() = {pyCOMPOs.keys()}")
	POSTPROC_Gd157(pyCOMPOs,ListeCOMPO,ListeAUTOP,name_study,suffixe,VISU_param,Nmin,S2_libs,S2_PCC,ssh_module,ssh_method,correlation)


	########################################################################################################################################################################################



		