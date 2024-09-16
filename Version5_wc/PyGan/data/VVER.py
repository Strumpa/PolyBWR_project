##################################################################################
#                                                                                #
# Description : PyGan scritp for VVER simulation with DRAGON5                    #
# Author      : L. Fede                                                          #
# Date        : 2023                                                             #
# Purpose     : main                                                             #
#                                                                                #
##################################################################################
#
########################################################################################################################################################################################
#                                #
#             IMPORTS            #
#                                #
##################################
import os, shutil, sys
import lifo
import lcm
import cle2000
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import matplotlib.pyplot as plt
from assertS import *

# POST-PROCESSING class
from POSTPROC import *

# --- CELL
# GEOMETRY class
from GEO_C import *
# MIXTURES class
from MIX_C import *
# TRACKING class
from TRK_C import *
# BU DEPLETION class
from BU_C import *

# --- HYBRIDE
# GEOMETRY class
from GEO_A_H import *
# MIXTURES class
from MIX_A_H import *
# TRACKING class
from TRK_A_H import *
# BU DEPLETION class
from BU_A_H import *

# --- ALAMOS
# MIXTURES 
from MIX_A_A import *
# TRACKING 
from TRK_A_A import *
# BU DEPLETION class
from BU_A_A import *

# --- NATIVE
# GEOMETRY class
from GEO_A_N import *
# MIXTURES 
from MIX_A_N import *
# TRACKING 
from TRK_A_N import *
# BU DEPLETION class
from BU_A_N import *

# --- OTHERS
from getLists import *
#
########################################################################################################################################################################################
#                                                        #
#         PARAMETERS OF CALCULUS : TO BE MODIFIED        #
#                                                        #
##########################################################
#
# Library = name of the used Library
Library = 'DLIB_J3p1p1'
#        - DLIB_J3p1p1
#        - DLIBENDFB7p1
#
# TYPE = type of calculus
#           - 'ALAMOS'
#           - 'HYBRIDE'
#           - 'NATIVE'
TYPE = 'HYBRIDE'
#
# name_geom = name_geom of the simulated case 
# - cells :
#           - 'KML2_22AU_22UOx'
#           - 'KML2_40AU_40UOx'
#           - 'KML2_398GO_44UOx'
#           - 'KML2_398GO_36Gd'
#           - 'KML2_39AWU_40UOx'
#           - 'KML2_39AWU_36UOx'
#           - 'KML2_39AWU_33Gd'
# - assemblies : 
#           - 'KML2_13AU_ASS'
#           - 'KML2_22AU_ASS'
#           - 'KML2_40AU_ASS'
#           - 'KML2_398GO_ASS'
#           - 'KML2_39AWU_ASS'
name_geom = 'KML2_22AU_22UOx'
#
# Multicompo = 1 if you want to generate a MULTICOMPO object 
#            = 0 either
Multicompo = 1
#
# suffixe = suffixe added to name_geom for creation of figures, MULTICOMPO and BU vector
suffixe = '_LIVRABLE_TEST'
#
# burnup_points = string to use predefine burnup points matching with SERPENT2 results
#        - 'test' : short depletion to test the code
#        - 'UOx'  : used for UOx fuel without Gd poison
#        - 'Gd'   : used for UOx fuel with Gd poison
#        - 'free' : modify the burnup points as you wish
burnup_points = 'UOx2_autop5'
#
# AUTOP: self-shielding parameter
#        - 'ALL'
#        - 'REGI'
#        - 'NONE'
AUTOP = 'REGI' ;
#
# SALT: tracking parameters for REL2005 calculus
an2d = 18 ;
#
# form = format of the figures generated with vverCOMP.py
#        - 'eps' (recommended)
#        - 'png'
form = 'eps'
#
# Nmin = minimum BU step for DELTA calculation and plot
Nmin = 0
#
# Visualisation options for post-processing functions
# 	= 1 if you want to plot and save the results
# 	= 0 either
# NB : if visu_XXX_Keff==1 or visu_XXX_BU==1 or visu_XXX_ISOTOPOESDENS==1 then visu_XXX must be equal to 1
visu_DRAGON=0
visu_SERPENT=0
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
case = 'CELL'
if name_geom == 'KML2_22AU_ASS' or name_geom == 'KML2_40AU_ASS' or name_geom == 'KML2_398GO_ASS' or name_geom == 'KML2_39AWU_ASS' or name_geom == '13UAX' or name_geom == 'V1_S1_ASS' or name_geom == 'V2_S1_ASS' or name_geom == 'KML2_13AU_ASS' :
	case = 'ASSEMBLY'
#
# name_alamos : name used to recover the ALAMOS geometry in PyGan/data folder
if case == 'ASSEMBLY' :
	if name_geom == 'KML2_22AU_ASS' or name_geom == 'KML2_40AU_ASS' or name_geom == 'KML2_13AU_ASS' :
		name_alamos_SSH_ap3 = './KML2_ARO_22AU_SSH.ap3'
		name_alamos_MOC_ap3 = './KML2_ARO_22AUX_MOC.ap3'
		name_alamos_SSH_za = './KML2_ARO_22AU_SSH.za'
		name_alamos_MOC_za = './KML2_ARO_22AUX_MOC.za'
	elif name_geom == 'KML2_398GO_ASS' :
		name_alamos_SSH_ap3 = './KML2_ARO_398GO_SSH.ap3'
		name_alamos_MOC_ap3 = './KML2_ARO_398GOX_MOC.ap3'
		name_alamos_SSH_za = './KML2_ARO_398GO_SSH.za'
		name_alamos_MOC_za = './KML2_ARO_398GOX_MOC.za'
	elif name_geom == 'KML2_39AWU_ASS' :
		name_alamos_SSH_ap3 = './KML2_ARO_39AWU_SSH.ap3'
		name_alamos_MOC_ap3 = './KML2_ARO_39AWUX_MOC.ap3'
		name_alamos_SSH_za = './KML2_ARO_39AWU_SSH.za'
		name_alamos_MOC_za = './KML2_ARO_39AWUX_MOC.za'
	elif name_geom == '13UAX' :
		name_alamos_MOC_ap3 = './13UAX_ARO_MOC.ap3'
		name_alamos_MOC_za = './13UAX_ARO_MOC.za'
		name_alamos_SSH_ap3 = './13UAX_ARO_MOC.ap3'
		name_alamos_SSH_za = './13UAX_ARO_MOC.za'
#
# names for exportation
name_compo = './_COMPO_'+name_geom+suffixe
name_BUvector = name_geom+suffixe
name_fig = name_geom+suffixe+'_2L.ps'
name_fig_SS = name_geom+suffixe+'_SS.ps'
name_fig_1L = name_geom+suffixe+'_1L.ps'
namFIG_SS_NXT=name_geom+suffixe+'_SS_NXT.ps'
namFIG_1L_NXT=name_geom+suffixe+'_1L_NXT.ps'
#
# Visualisation parameters stored in VISU_param
VISU_param=[visu_DRAGON,visu_SERPENT,visu_COMP,visu_DELTA]
#
# Creation of results directory
path=os.getcwd()

a=os.path.exists('VVERresults_PyGan_'+name_geom)
if a==False:
	os.mkdir('VVERresults_PyGan_'+name_geom)

SAVE_DIR='VVERresults_PyGan_'+name_geom+'/'+suffixe+'_postprocess'
a=os.path.exists(SAVE_DIR)
if a==False:
	os.mkdir(SAVE_DIR)

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
# Save BU evctor
os.chdir(path+'/'+'VVERresults_PyGan_'+name_geom)
np.savetxt(name_BUvector+'_BUvector.txt',ListeCOMPO)
os.chdir(path)
#
#
#
########################################################################################################################################################################################
#                                                 #
#             EXECUTION OF PROCEDURES             #
#                                                 #
###################################################
#
# --------------------------------------
#                 CELL                 
# --------------------------------------
if case == 'CELL' :

	# ---------- GEOMETRY
	[pyGEOM,pyGEOM_SS] = GEO_C("GEOM","GEOM_SS",name_geom)
	
	# ---------- MIXTURES
	pyMIX = MIX_C("LIBRARY",Library,name_geom)
	
	# ---------- TRACKING
	[pyUOX_TBH,pyFIG,pyTRACK,pyTF_EXC,pyUOX_TBH_SS,pyFIG_SS,pyTRACK_SS,pyTF_EXC_SS] = TRK_C("UOX_TBH",name_fig,"TRACK","TF_EXC","UOX_TBH_SS",name_fig_SS,"TRACK_SS","TF_EXC_SS",pyGEOM,pyGEOM_SS)
	
	# ---------- BU DEPLETION
	pyCOMPO = BU_C("COMPO",pyMIX,pyTRACK,pyTF_EXC,pyTRACK_SS,pyTF_EXC_SS,StepList,Multicompo,name_geom,name_compo)

# --------------------------------------
#               ASSEMBLIES                
# --------------------------------------
elif case == 'ASSEMBLY':

# --------------------------
#                    REL2005
# --------------------------
	if TYPE =='HYBRIDE' :

		# ---------- GEOMETRY 
		[pyGEOM_1L,pyGEOM_SS,pyGEOM_1L_NXT,pyGEOM_SS_NXT] = GEO_A_H("GEOM_1L","GEOM_SS","GEOM_1L_NXT","GEOM_SS_NXT",name_geom)
		
		# ---------- MIXTURES
		pyMIX = MIX_A_H("LIBRARY",Library,name_geom)
		
		# ---------- TRACKING
		[pyTRACK_2L,pyTF_EXC,pyTRACK_SS,pyFIG_2L,pyFIG_SS,pyTRACK_1L,pyFIG_1L,pyFIG_1L_NXT,pyFIG_SS_NXT] = TRK_A_H("TRACK_2L","TF_EXC","TRACK_SS",name_fig,name_fig_SS,"TRACK_1L",name_fig_1L,pyGEOM_SS,pyGEOM_1L,pyGEOM_SS_NXT,pyGEOM_1L_NXT,namFIG_SS_NXT,namFIG_1L_NXT,name_alamos_SSH_ap3,name_alamos_MOC_ap3,name_alamos_SSH_za,name_alamos_MOC_za,an2d)
		
		# ---------- BU DEPLETION
		pyCOMPO = BU_A_H("COMPO",pyMIX,pyTRACK_2L,pyTF_EXC,pyTRACK_1L,pyTRACK_SS,StepList,Multicompo,name_geom,name_compo,AUTOP)


# --------------------------
#                    ALAMOS
# --------------------------
	elif TYPE =='ALAMOS' :

		# ---------- MIXTURES
		pyMIX = MIX_A_A("LIBRARY",Library,name_geom)
		
		# ---------- TRACKING
		[pyTRACK,pyTF_EXC,pyTRACK_SS,pyTF_EXC_SS,pyFIG,pyFIG_SS] = TRK_A_A("TRACK","TF_EXC","TRACK_SS","TF_EXC_SS",name_fig,name_fig_SS,name_alamos_MOC_ap3,name_alamos_MOC_za,name_alamos_SSH_ap3,name_alamos_SSH_za)
		
		# ---------- BU DEPLETION
		pyCOMPO = BU_A_A("COMPO",pyMIX,pyTRACK,pyTF_EXC,pyTRACK_SS,pyTF_EXC_SS,StepList,Multicompo,name_geom,name_compo,AUTOP)

# --------------------------
#                    NATIVE
# --------------------------
	elif TYPE =='NATIVE' :

		# ---------- GEOMETRY 
		[pyGEOM_1L,pyGEOM_SS] = GEO_A_N("GEOM_1L","GEOM_SS",name_geom)
		
		# ---------- MIXTURES
		pyMIX = MIX_A_N("LIBRARY",Library,name_geom)
		
		# ---------- TRACKING
		[pyTRACK_1L,pyTRACK_SS,pyFIG_1L,pyFIG_SS] = TRK_A_N("TRACK_1L","TRACK_SS",name_fig,name_fig_SS,pyGEOM_1L,pyGEOM_SS,an2d)
		
		# ---------- BU DEPLETION
		pyCOMPO = BU_A_N("COMPO",pyMIX,pyTRACK_1L,pyTRACK_SS,StepList,Multicompo,name_geom,name_compo)

# -----------------------------------
#    POST-PROCESSING OF GLOBAL VALUES
# -----------------------------------

POSTPROC(pyCOMPO,ListeCOMPO,name_geom,suffixe,VISU_param,form,Nmin)


_
