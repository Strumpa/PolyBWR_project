# PyGan procedure to import compos from the BUScheme_Gd157 case and post-process / compare them to Serpent2 BU evolution files

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import lifo
import lcm
import cle2000

from Serpent2_edep_pcc import Serpent2_case 
from Serpent2_edep_pcc import SERPENT2_comparisons
from POSTPROC_BUScheme import DRAGON_case

post_treat_S2 = True

if post_treat_S2:
	CASES = []
	PyNjoy_cases = []
	oldlib_cases = []
	lib_name = "PyNjoy2016"
	S2_libs = ["PyNjoy2016","oldlib"]
	for lib_name in S2_libs:
		if lib_name == "PyNjoy2016":
			pcc_ids = [1,2]
			edep_ids = [0,1,2]
		elif lib_name == "oldlib":
			pcc_ids = [1]
			edep_ids = [0,1]
		for pcc in pcc_ids:
			for edep in edep_ids:
				CASES.append(Serpent2_case("HOM_UOX_Gd157", edep, pcc, lib_name, ["U235", "U238", "Gd157", "Gd158"]))
				if lib_name == "PyNjoy2016":
					PyNjoy_cases.append(CASES[-1])
				elif lib_name == "oldlib":
					oldlib_cases.append(CASES[-1])


	for S2_case in CASES:
		S2_case.plot_keff()
		S2_case.plot_concentrations(["Gd157", "Gd158"])
		#S2_case.rescale_BU("MWd/tU")	

	S2_comp_PyNjoy = SERPENT2_comparisons("Gd157 evolution schemes Pynjoy", PyNjoy_cases)
	S2_comp_PyNjoy.set_accepted_params_for_lib("PyNjoy2016", [0,1,2], [1,2,3,4,6])
	S2_comp_PyNjoy.set_ref_case(ref_lib="PyNjoy2016", ref_edep=0, ref_pcc=1)
	S2_comp_PyNjoy.plot_keffs()
	S2_comp_PyNjoy.compute_delta_keffs()
	S2_comp_PyNjoy.compute_delta_Niso(["U235", "U238", "Gd157", "Gd158"])
	S2_comp_PyNjoy.plot_delta_keffs()
	S2_comp_PyNjoy.plot_delta_Niso(["U235", "U238", "Gd157", "Gd158"])

	S2_comp_oldlib = SERPENT2_comparisons("Gd157 evolution schemes", oldlib_cases)
	S2_comp_oldlib.set_ref_case(ref_lib="oldlib", ref_edep=0, ref_pcc=1)
	S2_comp_oldlib.plot_keffs()
	S2_comp_oldlib.compute_delta_keffs()
	S2_comp_oldlib.compute_delta_Niso(["U235", "U238", "Gd157", "Gd158"])
	S2_comp_oldlib.plot_delta_keffs()
	S2_comp_oldlib.plot_delta_Niso(["U235", "U238", "Gd157", "Gd158"])


	print("Serpent2 comparisons done")



