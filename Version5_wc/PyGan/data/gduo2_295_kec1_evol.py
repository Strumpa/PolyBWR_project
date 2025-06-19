# PyGan script for gduo2_295_kec1 evolution study
# Authors : R.Guasch, A. Hébert
# Date : 2025-01-09
# Hypothesis : RSE method with CORR should elimiate mutual self shielding effect between U8 and Gd157 in group 100.

# Importing libraries
import sys
import numpy as np
import os, shutil
import lifo
import lcm
import cle2000
import matplotlib
import matplotlib.pyplot as plt

# Importing Procedures
# --- Geometry
from geom_gduo2 import *
# --- Materials
from lib_gduo2 import *
# --- Traking
from trk_gduo2 import *
# --- Burnup
from bu_gduo2 import *

# --- OTHERS
from getLists import *


# Input data
burnup_points = ["Gd_autop3", "Gd_autop4", "Gd_autop5"] # 'Gd_autop4' or 'Gd_autop3', refine if necessary

# Evaluation name
dlib_name = 'endfb8r1_295' # 'J311_295' or 'endfb8r1_295'
dlib_name_kerma = 'endfb81295K' # 'endfb81295K' or 'endfb81295K2'  or 'J311_295K
dlib_name_kerma2 = 'endfb81295K2'
self_shielding_methods = ['PT', 'RSE'] # 'PT' or 'RSE'
correlations = ['NOCORR', 'CORR'] # 'CORR' or 'NOCORR'

tracked_nuclides = ["U235", "U234", "U238", "Pu239", "Pu240", "Pu241", "Pu242", "Gd155", "Gd157", "Xe135", "Sm149"]

# Time integrator
time_integrators = ["EXTR"] # ["NOEX", "EXTR2", "CECE", "EXTR"] # "EXTR" + "CECM" --> to debug

glob_opt = "GLOB" # "GLOB" or "NOGL"

# Creation of results directory
path=os.getcwd()

a=os.path.exists(f'{path}/gduo2_295_kec1_evol_results/')
if a==False:
	os.mkdir(f'{path}/gduo2_295_kec1_evol_results/')

# Begin calculation procedure for gduo2_kec1 evolution study

# Geometry
pyGEOM = geom_gduo2()
# Tracking
pyUOX_TBH, pyTRACK, pyTF_EXC, pyTRACK_SS = trk_gduo2(pyGEOM)

for ssh_opt in self_shielding_methods:
	for correlation in correlations:
		if correlation == 'CORR':
			corr_name = 'C'
		elif correlation == 'NOCORR':
			corr_name = 'N'

		# Materials library
		
		pyLIB = lib_gduo2(dlib_name, ssh_opt, correlation, NG0=False)
		pyLIB_NG0 = lib_gduo2(dlib_name, ssh_opt, correlation, NG0=True)
		pyLIB_KERMA = lib_gduo2(dlib_name_kerma, ssh_opt, correlation, NG0=False)
		pyLIB_KERMA2 = lib_gduo2(dlib_name_kerma2, ssh_opt, correlation,  NG0=False)
		for time_integrator in time_integrators:
			for bu_pts in burnup_points:
				[BUList,SSHList,COMPOList] = getLists(bu_pts)
				# Burnup evolution
				name_compo_default = f"_CPO_{dlib_name}_{ssh_opt}_{corr_name}_SALT_{bu_pts}_KAPS_NODI_{time_integrator}_{glob_opt}"
				# StepList
				StepList = lcm.new('LCM','burnup_steps')
				StepList['BUlist']    = np.array(BUList, dtype='f')
				StepList['SSHlist'] = np.array(SSHList, dtype='f')
				StepList['COMPOlist'] = np.array(COMPOList, dtype='f')
				StepList.close() # close without erasing
				print(f"name_compo_default = {name_compo_default}")
				pyCOMPO = bu_gduo2("COMPO", pyLIB, pyTRACK, pyTF_EXC, pyTRACK_SS, StepList, name_compo_default, time_integrator, glob_opt)

				name_compo_NG0 = f"_CPO_{dlib_name}_NG0_{ssh_opt}_{corr_name}_SALT_{bu_pts}_KAPS_NODI_{time_integrator}_{glob_opt}"
				# StepList 2
				StepList2 = lcm.new('LCM','burnup_steps')
				StepList2['BUlist']    = np.array(BUList, dtype='f')
				StepList2['SSHlist'] = np.array(SSHList, dtype='f')
				StepList2['COMPOlist'] = np.array(COMPOList, dtype='f')
				StepList2.close() # close without erasing
				pyCOMPO_NG0 = bu_gduo2("COMPO", pyLIB_NG0, pyTRACK, pyTF_EXC, pyTRACK_SS, StepList2, name_compo_NG0, time_integrator, glob_opt)

				name_compo_KERMA = f"_CPO_{dlib_name_kerma}_{ssh_opt}_{corr_name}_SALT_{bu_pts}_KAPS_NODI_{time_integrator}_{glob_opt}"
				# StepList 3
				StepList3 = lcm.new('LCM','burnup_steps')
				StepList3['BUlist']    = np.array(BUList, dtype='f')
				StepList3['SSHlist'] = np.array(SSHList, dtype='f')
				StepList3['COMPOlist'] = np.array(COMPOList, dtype='f')
				StepList3.close() # close without erasing
				pyCOMPO_KERMA = bu_gduo2("COMPO", pyLIB_KERMA, pyTRACK, pyTF_EXC, pyTRACK_SS, StepList3, name_compo_KERMA,  time_integrator, glob_opt)
				
				name_compo_KERMA2 = f"_CPO_{dlib_name_kerma2}_{ssh_opt}_{corr_name}_SALT_{bu_pts}_KAPS_NODI_{time_integrator}_{glob_opt}"
				# StepList 4
				StepList4 = lcm.new('LCM','burnup_steps')
				StepList4['BUlist']    = np.array(BUList, dtype='f')
				StepList4['SSHlist'] = np.array(SSHList, dtype='f')
				StepList4['COMPOlist'] = np.array(COMPOList, dtype='f')
				StepList4.close() # close without erasing

				pyCOMPO_KERMA2 = bu_gduo2("COMPO", pyLIB_KERMA, pyTRACK, pyTF_EXC, pyTRACK_SS, StepList4, name_compo_KERMA2, time_integrator, glob_opt)


# End of calculation procedure for gduo2_kec1 evolution study