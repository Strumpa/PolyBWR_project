
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
from C7_hom import *

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
name_mix = "HOM_U5" # "AT10_C7_hom", "AT10_C7_noGd_hom", "AT10_C7_Gd154_hom", "AT10_C7_Gd155_hom", "AT10_C7_Gd156_hom", "AT10_C7_Gd157_hom", "AT10_C7_Gd158_hom", "AT10_C7_Gd160_hom",
# HOM_U5, HOM_U5_U8, HOM_UOX, HOM_UOX_Clad_noZr, HOM_UOX_Clad 

S2_libs = ["oldlib", "PyNjoy2016"]

ssh_module = "USS" #"USS", "AUTO"
ssh_method = "RSE" #"PT", "RSE", "SUBG" all supported for USS: but AUTO: only takes SUBG
# for RSE method : test with eps_RSE = 1.0E-1, 5.0E-2, 1.0E-2, 5.0E-3, 1.0E-3, 1.0E-4
correlation = "" # "CORR", "noCORR"

sat = "SATOFF" # "SAT", "" for short-lived nuclide saturation or not
depl_sol = "KAPS" # "RUNG" : 5th order Runge-Kutta (Cash-Karp) or "KAPS": 4th order Kaps-Rentrop

#
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
visu_SERPENT=1
visu_COMP=0
visu_DELTA=1
#
########################################################################################################################################################################################
#                                                 #
#         LOCAL PARAMETERS : DO NOT MODIFY        #
#                                                 #
###################################################
#
# case = 'CELL' / 'ASSEMBLY' - string used to determine wich python class must be called
case =  'HOM_CELL' #'CELL'

#
# names for exportation
if sat == "SAT":
	SAT = "_SAT"
elif sat == "SATOFF":
	SAT = "_SATOFF"
else:
	SAT = ""
if correlation == "CORR":
	correlation_name = "_CORR"
elif correlation == "noCORR":
	correlation_name = "_noCORR"
else:
	correlation_name = ""
name_compo = f"./_COMPO_{name_mix}_{suffixe}_{depl_sol}{SAT}_{ssh_module}_{ssh_method}{correlation_name}"
name_BUvector = f"{name_geom}_{name_mix}_{suffixe}"
name_fig = f"{name_geom}_{name_mix}{suffixe}_2L.ps"


#
# Visualisation parameters stored in VISU_param
VISU_param=[visu_DRAGON,visu_SERPENT,visu_COMP,visu_DELTA]
#
# Creation of results directory
path=os.getcwd()

a=os.path.exists(f"BWRresults_PyGan_{name_geom}")
if a==False:
	os.mkdir(f"BWRresults_PyGan_{name_geom}")

SAVE_DIR=f'BWRresults_PyGan_{name_geom}/{name_mix}/{suffixe}_postprocess/{ssh_module}_{ssh_method}/'
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
np.savetxt(name_BUvector+'_'+name_mix+'_BUvector.txt',ListeCOMPO)
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

if name_mix == "AT10_C7_hom":
	pyCOMPO = C7_hom("COMPO",StepList,name_compo,ssh_module,ssh_method,sat,depl_sol)
elif name_mix == "HOM_U5":
	pyCOMPO = HOM_U5("COMPO",StepList,name_compo,ssh_module,ssh_method,sat,depl_sol)
elif name_mix == "HOM_U5_U8":
	pyCOMPO = HOM_U5_U8("COMPO",StepList,name_compo,ssh_module,ssh_method,correlation,sat,depl_sol)

#POSTPROC_hom(pyCOMPO, ListeCOMPO, ListeAUTOP, name_geom, name_mix, suffixe, VISU_param, Nmin, GdCompo, S2_libs, ssh_module, ssh_method, correlation, sat, depl_sol)
MULTI_SERP_POSTPROC(pyCOMPO, ListeCOMPO, ListeAUTOP, name_geom, name_mix, suffixe, VISU_param, Nmin, S2_libs, ssh_module, ssh_method, correlation, depl_sol,  sat)