### post treat HOM Gd157 Very Beginning of Cycle results
# open compos from Dragon/Linux_*/
# open Serpent2 results

import serpentTools as st
import lcm
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

from postproc_cst_pow_evol import Serpent2_case as S2_case
from postproc_cst_pow_evol import DRAGON_case as D5_case
from postproc_cst_pow_evol import OpenMC_case as OMC_case


#### Running the post-treatment : options

eval_name = "endfb8r1" # or "J311" "endfb8r1"
origin_of_S2_data = "pynjoy2012_kerma" #"pynjoy2012_kerma" # or "pynjoy2016" or "NJOY2016"
origin_of_OMC_data = "NJOY2016" # no other option for now, check consistency with Serpent2 data and DRAGON data

# Serpent2 results
post_process_S2 = True # True or False
post_process_OMC = True # True or False
post_process_D5 = True # True or False


# Creation of results directory
path=os.getcwd()

a=os.path.exists(f'{path}/HOM_Gd157_VBOC_results_{eval_name}')
if a==False:
	os.mkdir(f'{path}/HOM_Gd157_VBOC_results_{eval_name}')

save_dir = f'{path}/HOM_Gd157_VBOC_results_{eval_name}'

save_dir_Serpent2 = f'{save_dir}/Serpent2'
a=os.path.exists(save_dir_Serpent2)
if a==False:
    os.makedirs(save_dir_Serpent2)

save_dir_DRAGON = f'{save_dir}/DRAGON'
a=os.path.exists(save_dir_DRAGON)
if a==False:
    os.makedirs(save_dir_DRAGON)

save_dir_OpenMC = f'{save_dir}/OpenMC'
a=os.path.exists(save_dir_OpenMC)
if a==False:
    os.makedirs(save_dir_OpenMC)

## Discriminate between burnup points used in DRAGON5 : VBOC_OMC and VBOC
# VBOC_OMC : 14 burnup points for comparison with OpenMC
# VBOC : 40 burnup points for comparison with Serpent2
    
save_dir_comparisons_VBOC_OMC = f'{save_dir}/Comparisons_BU_VBOC_OMC'
a=os.path.exists(save_dir_comparisons_VBOC_OMC)
if a==False:
    os.makedirs(save_dir_comparisons_VBOC_OMC)

save_dir_comparisons_VBOC = f'{save_dir}/Comparisons_BU_VBOC'
a=os.path.exists(save_dir_comparisons_VBOC)
if a==False:
    os.makedirs(save_dir_comparisons_VBOC)


tracked_nuclides = ["Gd157", "Gd158", "U235", "U238", "Pu239", "Xe135", "Sm149"]# ["U235", "U234", "U238", "Pu239", "Pu240", "Pu241", "Pu242", "Gd157", "Gd158", "Xe135", "Sm149"]

if post_process_S2:
    # Serpent2 results
    # ----------------
    Serpent2_case_edep0 = S2_case(case_name = "HOM_Gd157_VBOC_OMC", lib_name=f"{eval_name}_{origin_of_S2_data}", 
                            edep_id=0, areQfissSet=False, isEcaptSet=False, 
                            pcc_id=0, 
                            specific_power=38.6, 
                            tracked_nuclides=tracked_nuclides, save_dir=save_dir_Serpent2)	


    Serpent2_case_edep0_setQfiss = S2_case(case_name = "HOM_Gd157_VBOC_OMC", lib_name=f"{eval_name}_{origin_of_S2_data}",
                            edep_id=0, areQfissSet=True, isEcaptSet=False, 
                            pcc_id=0, 
                            specific_power=38.6, 
                            tracked_nuclides=tracked_nuclides, save_dir=save_dir_Serpent2)

if post_process_OMC:
    ## Create OpenMC cases for fissq normalization, with and without setting the fission q-values
    OMC_case_edep0 = OMC_case(case_name = "HOM_Gd157_VBOC", lib_name=f"{eval_name}_{origin_of_OMC_data}",
                            edep_id="fissq", areQfissSet=False, integrator = "Predictor", specific_power=38.6,
                            tracked_nuclides=tracked_nuclides, save_dir=save_dir_OpenMC)
    #OMC_case_fissq.plot_keff()
    OMC_case_edep0_setQfiss = OMC_case(case_name = "HOM_Gd157_VBOC", lib_name=f"{eval_name}_{origin_of_OMC_data}",
                            edep_id="fissq", areQfissSet=True, integrator = "Predictor", specific_power=38.6,
                            tracked_nuclides=tracked_nuclides, save_dir=save_dir_OpenMC)
    #OMC_case_fissq_set_values.plot_keff()

if post_process_D5:
    # DRAGON results
    # --------------
        
    ## Compo on 14 VBOC BU points for OpenMC comparison ("VBOC_OMC"),
    # PCC0 = no linear extrapolation of rates over a time stage, if not specified, using EXTR (default, linear extrapolation)
    name_compo_OMC = f"COMPO_HOM_Gd157_VBOC_OMC_{eval_name}_295"
    name_compo_OMC_PCC0 = f"COMPO_HOM_Gd157_VBOC_OMC_{eval_name}_295_PCC0"
    ## Compos obtained with modified DEPL structure
    name_compo_OMC_NG0 = f"COMPO_HOM_Gd157_VBOC_OMC_{eval_name}_295_NG0"
    name_compo_OMC_NG0_PCC0 = f"COMPO_HOM_Gd157_VBOC_OMC_{eval_name}_295_NG0_PCC0"

    ## COMPOs on 40 VBOC BU points for finer D5 calculation : to compare with S2/OMC results for 14 VBOC BU points
    compo_VBOC = f"COMPO_HOM_Gd157_VBOC_{eval_name}_295"
    compo_VBOC_PCC0 = f"COMPO_HOM_Gd157_VBOC_{eval_name}_295_PCC0"
    compo_VBOC_NG0 = f"COMPO_HOM_Gd157_VBOC_{eval_name}_295_NG0"
    compo_VBOC_NG0_PCC0 = f"COMPO_HOM_Gd157_VBOC_{eval_name}_295_NG0_PCC0"

    
    ## Load the COMPOs
    os.chdir("D5_COMPOs_path/")
    #print(os.listdir())
    #print(name_compo)
    pyCOMPO_VBOC_OMC = lcm.new('LCM_INP', name_compo_OMC, impx=0)
    pyCOMPO_VBOC_OMC_PCC0 = lcm.new('LCM_INP', name_compo_OMC_PCC0, impx=0)
    pyCOMPO_VBOC_OMC_NG0 = lcm.new('LCM_INP', name_compo_OMC_NG0, impx=0)
    pyCOMPO_VBOC_OMC_NG0_PCC0 = lcm.new('LCM_INP', name_compo_OMC_NG0_PCC0, impx=0)

    pyCOMPO_VBOC = lcm.new('LCM_INP', compo_VBOC, impx=0)
    pyCOMPO_VBOC_PCC0 = lcm.new('LCM_INP', compo_VBOC_PCC0, impx=0)
    pyCOMPO_VBOC_NG0 = lcm.new('LCM_INP', compo_VBOC_NG0, impx=0)
    pyCOMPO_VBOC_NG0_PCC0 = lcm.new('LCM_INP', compo_VBOC_NG0_PCC0, impx=0)

    
    print(f"pyCOMPO['EDIBU']['GLOBAL']['pval000000001'] are = {pyCOMPO_VBOC_OMC['EDIBU']['GLOBAL']['pval00000001']}")
    os.chdir(path)
    VBOC = {"BU":[0.0, 25.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0,
                    750.0, 800.0, 850.0, 900.0, 950.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0,    
                    5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0 ],
                "AUTOP":[25.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0,
                    750.0, 800.0, 850.0, 900.0, 950.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0,    
                    5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0 ],
                "COMPO":[0.0, 25.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0,
                    750.0, 800.0, 850.0, 900.0, 950.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0,    
                    5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0 ]}
    VBOC_OMC = {"BU" : [0.0, 25.0, 50.0, 100.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 2000.0, 4000.0, 6000.0, 8000.0, 10000.0],
                    "AUTOP" : [25.0, 50.0, 100.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 2000.0, 4000.0, 6000.0, 8000.0, 10000.0],
                    "COMPO" : [0.0, 25.0, 50.0, 100.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 2000.0, 4000.0, 6000.0, 8000.0, 10000.0]}
    
    # Create D5 cases for the VBOC_OMC burnup points cases
    DRAGON_case_OMC = D5_case(pyCOMPO=pyCOMPO_VBOC_OMC, dlib_name=f"{eval_name}_295", bu_points="VBOC_OMC",
                            ssh_opt="PT", correlation="CORR", sat="def", depl_sol="def",
                            tracked_nuclides=tracked_nuclides, BU_lists=VBOC_OMC, save_dir=save_dir_DRAGON)
    
    DRAGON_case_OMC_PCC0 = D5_case(pyCOMPO=pyCOMPO_VBOC_OMC_PCC0, dlib_name=f"{eval_name}_295", bu_points="VBOC_OMC", 
                        ssh_opt="PT", correlation="CORR", sat="def", depl_sol="def", 
                        tracked_nuclides=tracked_nuclides, BU_lists=VBOC_OMC, save_dir=save_dir_DRAGON)

    DRAGON_case_OMC_NG0 = D5_case(pyCOMPO=pyCOMPO_VBOC_OMC_NG0, dlib_name=f"{eval_name}_295_NG0", bu_points="VBOC_OMC",
                            ssh_opt="PT", correlation="CORR", sat="def", depl_sol="def",
                            tracked_nuclides=tracked_nuclides, BU_lists=VBOC_OMC, save_dir=save_dir_DRAGON)
    
    DRAGON_case_OMC_NG0_PCC0 = D5_case(pyCOMPO=pyCOMPO_VBOC_OMC_NG0_PCC0, dlib_name=f"{eval_name}_295_NG0", bu_points="VBOC_OMC",
                            ssh_opt="PT", correlation="CORR", sat="def", depl_sol="def",
                            tracked_nuclides=tracked_nuclides, BU_lists=VBOC_OMC, save_dir=save_dir_DRAGON)
    
    # Create D5 cases for the VBOC burnup points cases
    DRAGON_case_VBOC = D5_case(pyCOMPO=pyCOMPO_VBOC, dlib_name=f"{eval_name}_295", bu_points="VBOC",
                            ssh_opt="PT", correlation="CORR", sat="def", depl_sol="def",
                            tracked_nuclides=tracked_nuclides, BU_lists=VBOC, save_dir=save_dir_DRAGON)
    
    DRAGON_case_VBOC_PCC0 = D5_case(pyCOMPO=pyCOMPO_VBOC_PCC0, dlib_name=f"{eval_name}_295", bu_points="VBOC",
                            ssh_opt="PT", correlation="CORR", sat="def", depl_sol="def",
                            tracked_nuclides=tracked_nuclides, BU_lists=VBOC, save_dir=save_dir_DRAGON)
    
    DRAGON_case_VBOC_NG0 = D5_case(pyCOMPO=pyCOMPO_VBOC_NG0, dlib_name=f"{eval_name}_295_NG0", bu_points="VBOC",
                            ssh_opt="PT", correlation="CORR", sat="def", depl_sol="def",
                            tracked_nuclides=tracked_nuclides, BU_lists=VBOC, save_dir=save_dir_DRAGON)
    
    DRAGON_case_VBOC_NG0_PCC0 = D5_case(pyCOMPO=pyCOMPO_VBOC_NG0_PCC0, dlib_name=f"{eval_name}_295_NG0", bu_points="VBOC",
                            ssh_opt="PT", correlation="CORR", sat="def", depl_sol="def",
                            tracked_nuclides=tracked_nuclides, BU_lists=VBOC, save_dir=save_dir_DRAGON)
    

### Analysis performed on the VBOC_OMC burnup points for ALL three codes
# TO DO 1 : compare DRAGON5 results with and without (n,gamma) Q-values, with and without PCC0     

# Comparison between DRAGON5 results
# ---------------------------------
# Compare keffs
delta_keff_NG0_withPCC0 = (DRAGON_case_OMC_NG0_PCC0.DRAGON_Keff - DRAGON_case_OMC_PCC0.DRAGON_Keff)*1e5 # pcm
delta_keff_NG0 = (DRAGON_case_OMC_NG0.DRAGON_Keff - DRAGON_case_OMC.DRAGON_Keff)*1e5 # pcm
# plot delta keffs

plt.figure()
plt.plot(DRAGON_case_OMC.DRAGON_BU, delta_keff_NG0, label = f"D5 : Q(n,\\gamma) set to 0 - default DEPL, EXTR", marker = "x", linestyle = "--")
plt.plot(DRAGON_case_OMC.DRAGON_BU, delta_keff_NG0_withPCC0, label = f"D5 : Q(n,\\gamma) set to 0 - default DEPL, NOEX", marker = "x", linestyle = "--")
plt.title(f"$\\Delta$ Keff between DRAGON5 cases with and without (n,\\gamma) Q-values, {eval_name} \n on VBOC_OMC burnup points")
plt.xlabel("Burnup [MWd/kgU]")
plt.ylabel("Keff difference [pcm]")
plt.legend()
plt.grid()
plt.savefig(f"{save_dir_DRAGON}/DRAGON_Keff_diff_{DRAGON_case_OMC.draglib_name}_{DRAGON_case_OMC_NG0.draglib_name}.png")
plt.close()

# compare isotopes densities
for iso in tracked_nuclides:
    delta_Niso_NG0 = [(DRAGON_case_OMC_NG0.DRAGON_ISOTOPESDENS[iso][idx] - DRAGON_case_OMC.DRAGON_ISOTOPESDENS[iso][idx]) * 100 / DRAGON_case_OMC.DRAGON_ISOTOPESDENS[iso][idx]
                    if DRAGON_case_OMC.DRAGON_ISOTOPESDENS[iso][idx] != 0 else 0
                    for idx in range(len(DRAGON_case_OMC.DRAGON_ISOTOPESDENS[iso]))]
    delta_Niso_NG0_withPCC0 = [(DRAGON_case_OMC_NG0_PCC0.DRAGON_ISOTOPESDENS[iso][idx] - DRAGON_case_OMC_PCC0.DRAGON_ISOTOPESDENS[iso][idx]) * 100 / DRAGON_case_OMC_PCC0.DRAGON_ISOTOPESDENS[iso][idx]
                    if DRAGON_case_OMC_PCC0.DRAGON_ISOTOPESDENS[iso][idx] != 0 else 0
                    for idx in range(len(DRAGON_case_OMC_PCC0.DRAGON_ISOTOPESDENS[iso]))]
    #print(f"For iso = {iso}, deltaNiso = {delta_Niso_NG0} %")
    fig, ax1 = plt.subplots()
    ax1.plot(DRAGON_case_OMC.DRAGON_BU, DRAGON_case_OMC.DRAGON_ISOTOPESDENS[iso], label = f"D5 : Q(n,\\gamma) set to 0 - default DEPL, EXTR", marker = "x", linestyle = "--")
    ax1.plot(DRAGON_case_OMC.DRAGON_BU, DRAGON_case_OMC_NG0.DRAGON_ISOTOPESDENS[iso], label = f"D5 : Q(n,\\gamma) set to 0 - default DEPL, NOEX", marker = "x", linestyle = "--")
    ax1.set_xlabel("Burnup [MWd/kgU]")
    ax1.set_ylabel(f"{iso} density [iso/cm*b]", color='tab:blue')

    ax2 = ax1.twinx()
    ax2.plot(DRAGON_case_OMC.DRAGON_BU, delta_Niso_NG0, label = f"D5 : Q(n,\\gamma) set to 0 - default DEPL, EXTR", marker = "x", linestyle = "--", color='tab:red')
    ax2.plot(DRAGON_case_OMC.DRAGON_BU, delta_Niso_NG0_withPCC0, label = f"D5 : Q(n,\\gamma) set to 0 - default DEPL, NOEX", marker = "x", linestyle = "--", color='tab:green')
    ax2.set_ylabel(f"{iso} density difference [%]", color='tab:red')
    plt.title(f"$\\Delta$ N{iso} between DRAGON5 cases with and without (n,\\gamma) Q-values \n {eval_name}, on VBOC_OMC burnup points")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_DRAGON}/DRAGON_{iso}_diff_{DRAGON_case_OMC.draglib_name}_{DRAGON_case_OMC_NG0.draglib_name}.png")
    plt.close()

        
# TO DO 2 : compare Seprent2 edep0 with set Qfiss + OpenMC with set Qfiss with DRAGON5 with and without (n,gamma) Q-values, with and without PCC0
    
# Comparison between Serpent2 and DRAGON5 results
# ----------------------------------------------
# S2 with PCC0 : compare impact of setting Qfiss and impact of EXTR or NOEX on DRAGON5 results
delta_keff_D5_S2 = (DRAGON_case_OMC.DRAGON_Keff - Serpent2_case_edep0.keffs)*1e5 # pcm
delta_keff_D5_NG0_S2 = (DRAGON_case_OMC_NG0.DRAGON_Keff - Serpent2_case_edep0.keffs)*1e5 # pcm

delta_keff_D5_S2_setQfiss = (DRAGON_case_OMC.DRAGON_Keff - Serpent2_case_edep0_setQfiss.keffs)*1e5 # pcm
delta_keff_D5_NG0_S2_setQfiss = (DRAGON_case_OMC_NG0.DRAGON_Keff - Serpent2_case_edep0_setQfiss.keffs)*1e5 # pcm

delta_keff_D5_S2_PCC0 = (DRAGON_case_OMC_PCC0.DRAGON_Keff - Serpent2_case_edep0.keffs)*1e5 # pcm
delta_keff_D5_NG0_S2_PCC0 = (DRAGON_case_OMC_NG0_PCC0.DRAGON_Keff - Serpent2_case_edep0.keffs)*1e5 # pcm

delta_keff_D5_S2_PCC0_setQfiss = (DRAGON_case_OMC_PCC0.DRAGON_Keff - Serpent2_case_edep0_setQfiss.keffs)*1e5 # pcm
delta_keff_D5_NG0_S2_PCC0_setQfiss = (DRAGON_case_OMC_NG0_PCC0.DRAGON_Keff - Serpent2_case_edep0_setQfiss.keffs)*1e5 # pcm

# plot delta keffs
plt.figure(figsize=(10,6))
plt.plot(Serpent2_case_edep0.BU, delta_keff_D5_S2, label = f"D5 - S2", marker = "x", linestyle = "--") # normal DEPL, EXTR - default Qfiss, PCC0
plt.plot(Serpent2_case_edep0.BU, delta_keff_D5_NG0_S2, label = f"D5 (NG0) - S2", marker = "D", linestyle = "--") # NG0 DEPL, EXTR - default Qfiss, PCC0

plt.plot(Serpent2_case_edep0.BU, delta_keff_D5_S2_setQfiss, label = f"D5 - S2 setQfiss", marker = "x", linestyle = "--") # normal DEPL, EXTR - setQfiss, PCC0
plt.plot(Serpent2_case_edep0.BU, delta_keff_D5_NG0_S2_setQfiss, label = f"D5 (NG0) - S2 setQfiss", marker = "D", linestyle = "--") # NG0 DEPL, EXTR - setQfiss, PCC0

plt.plot(Serpent2_case_edep0.BU, delta_keff_D5_S2_PCC0, label = f"D5 PCC0 - S2", marker = "x", linestyle = "--") # normal DEPL, NOEX - default Qfiss, PCC0
plt.plot(Serpent2_case_edep0.BU, delta_keff_D5_NG0_S2_PCC0, label = f"D5 (NG0) PCC0 - S2", marker = "D", linestyle = "--") # NG0 DEPL, NOEX - default Qfiss, PCC0

plt.plot(Serpent2_case_edep0.BU, delta_keff_D5_S2_PCC0_setQfiss, label = f"D5 PCC0 - S2 setQfiss", marker = "x", linestyle = "--") # normal DEPL, NOEX - setQfiss, PCC0
plt.plot(Serpent2_case_edep0.BU, delta_keff_D5_NG0_S2_PCC0_setQfiss, label = f"D5 (NG0) PCC0- S2 setQfiss", marker = "D", linestyle = "--") # NG0 DEPL, NOEX - setQfiss, PCC0, this should be the most promising candidate

plt.plot(Serpent2_case_edep0.BU, 300*np.ones_like(Serpent2_case_edep0.BU),  color = "red", linestyle = "-.")
plt.plot(Serpent2_case_edep0.BU, -300*np.ones_like(Serpent2_case_edep0.BU), color = "red", linestyle = "-.")
plt.title(f"$\\Delta$ Keff between DRAGON5 and Serpent2  \n with and without (n,$\\gamma$) Q-values, {eval_name} \n on VBOC_OMC burnup points")
plt.xlabel("Burnup [MWd/kgU]")
plt.ylabel("Keff difference [pcm]")
plt.legend()
plt.grid()
plt.savefig(f"{save_dir_comparisons_VBOC_OMC}/D5_S2_Keff_diff_{Serpent2_case_edep0.case_name}_{DRAGON_case_OMC.draglib_name}.png")
plt.close()

# Comparison between OpenMC and DRAGON5 results
# ----------------------------------------------
# OpenMC results with Predictor = PCC0 : compare impact of setting Qfiss and impact of EXTR or NOEX on DRAGON5 results
delta_keff_D5_OMC = (DRAGON_case_OMC.DRAGON_Keff - OMC_case_edep0.keffs)*1e5 # pcm
delta_keff_D5_NG0_OMC = (DRAGON_case_OMC_NG0.DRAGON_Keff - OMC_case_edep0.keffs)*1e5 # pcm

delta_keff_D5_OMC_setQfiss = (DRAGON_case_OMC.DRAGON_Keff - OMC_case_edep0_setQfiss.keffs)*1e5 # pcm
delta_keff_D5_NG0_OMC_setQfiss = (DRAGON_case_OMC_NG0.DRAGON_Keff - OMC_case_edep0_setQfiss.keffs)*1e5 # pcm

delta_keff_D5_OMC_PCC0 = (DRAGON_case_OMC_PCC0.DRAGON_Keff - OMC_case_edep0.keffs)*1e5 # pcm
delta_keff_D5_NG0_OMC_PCC0 = (DRAGON_case_OMC_NG0_PCC0.DRAGON_Keff - OMC_case_edep0.keffs)*1e5 # pcm

delta_keff_D5_OMC_PCC0_setQfiss = (DRAGON_case_OMC_PCC0.DRAGON_Keff - OMC_case_edep0_setQfiss.keffs)*1e5 # pcm
delta_keff_D5_NG0_OMC_PCC0_setQfiss = (DRAGON_case_OMC_NG0_PCC0.DRAGON_Keff - OMC_case_edep0_setQfiss.keffs)*1e5 # pcm

plt.figure(figsize=(10,6))
plt.plot(OMC_case_edep0.BU, delta_keff_D5_OMC, label = f"D5 - OMC", marker = "x", linestyle = "--")
plt.plot(OMC_case_edep0.BU, delta_keff_D5_NG0_OMC, label = f"D5 (NG0) - OMC", marker = "D", linestyle = "--")

plt.plot(OMC_case_edep0.BU, delta_keff_D5_OMC_setQfiss, label = f"D5 - OMC setQfiss", marker = "x", linestyle = "--")
plt.plot(OMC_case_edep0.BU, delta_keff_D5_NG0_OMC_setQfiss, label = f"D5 (NG0) - OMC setQfiss", marker = "D", linestyle = "--")

plt.plot(OMC_case_edep0.BU, delta_keff_D5_OMC_PCC0, label = f"D5 PCC0 - OMC", marker = "x", linestyle = "--")
plt.plot(OMC_case_edep0.BU, delta_keff_D5_NG0_OMC_PCC0, label = f"D5 (NG0) PCC0 - OMC", marker = "D", linestyle = "--")

plt.plot(OMC_case_edep0.BU, delta_keff_D5_OMC_PCC0_setQfiss, label = f"D5 PCC0 - OMC setQfiss", marker = "x", linestyle = "--")
plt.plot(OMC_case_edep0.BU, delta_keff_D5_NG0_OMC_PCC0_setQfiss, label = f"D5 (NG0) PCC0 - OMC setQfiss", marker = "D", linestyle = "--")

plt.plot(OMC_case_edep0.BU, 300*np.ones_like(OMC_case_edep0.BU),  color = "red", linestyle = "-.")
plt.plot(OMC_case_edep0.BU, -300*np.ones_like(OMC_case_edep0.BU), color = "red", linestyle = "-.")
plt.title(f"$\\Delta$ Keff between DRAGON5 and OpenMC \n with and without (n,$\\gamma$) Q-values, {eval_name} \n on VBOC_OMC burnup points")
plt.xlabel("Burnup [MWd/kgU]")
plt.ylabel("Keff difference [pcm]")
plt.legend()
plt.grid()
plt.savefig(f"{save_dir_comparisons_VBOC_OMC}/D5_OMC_Keff_diff_{OMC_case_edep0.case_name}_{DRAGON_case_OMC.draglib_name}.png")
plt.close()


## Select the best candidate for the comparison between Serpent2 and DRAGON5 results, plot the results
plt.figure(figsize=(10,6))
# plot delta keffs for D5 EXTR NG0 vs S2 setQfiss
#plt.plot(Serpent2_case_edep0.BU, delta_keff_D5_NG0_S2_setQfiss, label = f"D5 (NG0) - S2 setQfiss", marker = "D", linestyle = "--") # NG0 DEPL, EXTR - setQfiss, pcc0
# plot delta keffs for D5 NOEX vs S2 setQfiss
plt.plot(Serpent2_case_edep0.BU, delta_keff_D5_NG0_S2_PCC0_setQfiss, label = f"D5 (NG0) NOEX - S2 setQfiss", marker = "D", linestyle = "--") # NG0 DEPL, NOEX - setQfiss, pcc0
## Select the best candidate for the comparison between OpenMC and DRAGON5 results, plot the results
# plot delta keffs for D5 EXTR NG0 vs OMC setQfiss
plt.plot(OMC_case_edep0.BU, delta_keff_D5_NG0_OMC_PCC0, label = f"D5 (NG0) NOEX - OMC", marker = "D", linestyle = "--") # NG0 DEPL, NOEX - Predictor
# plot delta keffs for D5 NOEX vs OMC setQfiss
plt.plot(OMC_case_edep0.BU, delta_keff_D5_NG0_OMC_PCC0_setQfiss, label = f"D5 (NG0) NOEX - OMC setQfiss", marker = "D", linestyle = "--") # NG0 DEPL, NOEX - setQfiss, Predictor

plt.plot(Serpent2_case_edep0.BU, 300*np.ones_like(Serpent2_case_edep0.BU),  color = "red", linestyle = "-.")
plt.plot(Serpent2_case_edep0.BU, -300*np.ones_like(Serpent2_case_edep0.BU), color = "red", linestyle = "-.")

plt.title(f"$\\Delta$ Keff between D5 and S2/OMC cases \n modified DEPL structure, {eval_name} \n on VBOC_OMC burnup points")
plt.xlabel(f"Burnup [{Serpent2_case_edep0.unitsBU}]")
plt.ylabel("Keff difference [pcm]")
plt.legend()
plt.grid()
plt.savefig(f"{save_dir_comparisons_VBOC_OMC}/BestCandidates_Keff_diff_{Serpent2_case_edep0.case_name}_{DRAGON_case_OMC_NG0.draglib_name}.png")


# For the best candidate, plot the isotopes densities

for iso in tracked_nuclides:

    # Compare iso number densities
    plt.figure()
    plt.plot(DRAGON_case_OMC_NG0_PCC0.DRAGON_BU, DRAGON_case_OMC_NG0_PCC0.DRAGON_ISOTOPESDENS[iso], label = f"D5 (NG0) PCC0", marker = "D", linestyle = "--")
    plt.plot(OMC_case_edep0_setQfiss.BU, OMC_case_edep0_setQfiss.Ni[iso], label = f"OMC setQfiss", marker = "o", linestyle = "--")
    plt.plot(Serpent2_case_edep0_setQfiss.BU, Serpent2_case_edep0_setQfiss.Ni[iso], label = f"S2 setQfiss", marker = "v", linestyle = "--")
    plt.title(f"{iso} number density comparison between D5 and S2/OMC cases, {eval_name}")
    plt.xlabel(f"Burnup [{Serpent2_case_edep0.unitsBU}]")
    plt.ylabel("Number density [atom/b-cm]")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_comparisons_VBOC_OMC}/COMP_{iso}_Niso_{Serpent2_case_edep0.case_name}_{DRAGON_case_OMC_NG0.draglib_name}.png")
    plt.close()

    # Compute relative and absolute differences
    rel_delta_Niso_OMC = [(DRAGON_case_OMC_NG0_PCC0.DRAGON_ISOTOPESDENS[iso][idx] - OMC_case_edep0_setQfiss.Ni[iso][idx]) * 100 / OMC_case_edep0_setQfiss.Ni[iso][idx]
                        if OMC_case_edep0_setQfiss.Ni[iso][idx] != 0 else 0
                        for idx in range(len(OMC_case_edep0_setQfiss.Ni[iso]))]
    rel_delta_Niso_S2 = [(DRAGON_case_OMC_NG0_PCC0.DRAGON_ISOTOPESDENS[iso][idx] - Serpent2_case_edep0_setQfiss.Ni[iso][idx]) * 100 / Serpent2_case_edep0_setQfiss.Ni[iso][idx]
                        if Serpent2_case_edep0_setQfiss.Ni[iso][idx] != 0 else 0
                        for idx in range(len(Serpent2_case_edep0_setQfiss.Ni[iso]))]

    abs_delta_Niso_OMC = DRAGON_case_OMC_NG0_PCC0.DRAGON_ISOTOPESDENS[iso] - OMC_case_edep0_setQfiss.Ni[iso]

    abs_delta_Niso_S2 = DRAGON_case_OMC_NG0_PCC0.DRAGON_ISOTOPESDENS[iso] - Serpent2_case_edep0_setQfiss.Ni[iso]

    ## Plot relative differences
    plt.figure()
    plt.plot(OMC_case_edep0_setQfiss.BU, rel_delta_Niso_OMC, label = f"OMC setQfiss", marker = "o", linestyle = "--")
    plt.plot(Serpent2_case_edep0_setQfiss.BU, rel_delta_Niso_S2, label = f"S2 setQfiss", marker = "v", linestyle = "--")
    plt.title(f"Relative differences in {iso} number density \n between D5 and S2/OMC cases, {eval_name}")
    plt.plot(DRAGON_case_OMC_NG0_PCC0.DRAGON_BU, 2.0*np.ones_like(DRAGON_case_OMC_NG0_PCC0.DRAGON_BU), color = "red", linestyle = "-.")
    plt.plot(DRAGON_case_OMC_NG0_PCC0.DRAGON_BU, -2.0*np.ones_like(DRAGON_case_OMC_NG0_PCC0.DRAGON_BU), color = "red", linestyle = "-.")
    plt.xlabel(f"Burnup [{Serpent2_case_edep0.unitsBU}]")
    plt.ylabel("Relative difference [%]")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_comparisons_VBOC_OMC}/COMP_{iso}_Niso_rel_{Serpent2_case_edep0.case_name}_{DRAGON_case_OMC_NG0.draglib_name}.png")
    plt.close()

    ## Plot absolute differences
    plt.figure()
    plt.plot(OMC_case_edep0_setQfiss.BU, abs_delta_Niso_OMC, label = f"OMC setQfiss", marker = "o", linestyle = "--")
    plt.plot(Serpent2_case_edep0_setQfiss.BU, abs_delta_Niso_S2, label = f"S2 setQfiss", marker = "v", linestyle = "--")
    #plt.plot(DRAGON_case_OMC_NG0_PCC0.DRAGON_BU, 2.0*np.ones_like(DRAGON_case_OMC_NG0_PCC0.DRAGON_BU), color = "red", linestyle = "-.")
    #plt.plot(DRAGON_case_OMC_NG0_PCC0.DRAGON_BU, -2.0*np.ones_like(DRAGON_case_OMC_NG0_PCC0.DRAGON_BU), color = "red", linestyle = "-.")
    plt.title(f"Absolute differences in {iso} number density \n between D5 and S2/OMC cases, {eval_name}")
    plt.xlabel(f"Burnup [{Serpent2_case_edep0.unitsBU}]")
    plt.ylabel("Absolute difference [atom/b-cm]")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_comparisons_VBOC_OMC}/COMP_{iso}_abs_{Serpent2_case_edep0.case_name}_{DRAGON_case_OMC_NG0.draglib_name}.png")
    plt.close()

# Now compare when using 40 VBOC BU points for DRAGON5 and 14 VBOC BU points for S2/OMC
    
# Prioritize "best" candidate for the comparison with S2/OpenMC
# ------------------------------------------------------------
# Interpolate the DRAGON5 results on the S2/OMC burnup points
keffs_D5_NG0_PCC0_interp = np.interp(VBOC_OMC["BU"], DRAGON_case_VBOC_NG0_PCC0.DRAGON_BU, DRAGON_case_VBOC_NG0_PCC0.DRAGON_Keff)
keffs_D5_NG0_interp = np.interp(VBOC_OMC["BU"], DRAGON_case_VBOC_NG0.DRAGON_BU, DRAGON_case_VBOC_NG0.DRAGON_Keff)

# Compare keffs

delta_keff_D5_NG0_S2_setQfiss = (keffs_D5_NG0_interp - Serpent2_case_edep0_setQfiss.keffs)*1e5 # pcm
delta_keff_D5_NG0_OMC_setQfiss = (keffs_D5_NG0_interp - OMC_case_edep0_setQfiss.keffs)*1e5 # pcm

delta_keff_D5_NG0_PCC0_S2_setQfiss = (keffs_D5_NG0_PCC0_interp - Serpent2_case_edep0_setQfiss.keffs)*1e5 # pcm
delta_keff_D5_NG0_PCC0_OMC_setQfiss = (keffs_D5_NG0_PCC0_interp - OMC_case_edep0_setQfiss.keffs)*1e5 # pcm

## Interpolate NGd157 number densities

NGd157_D5_NG0_interp = np.interp(VBOC_OMC["BU"], DRAGON_case_VBOC_NG0_PCC0.DRAGON_BU, DRAGON_case_VBOC_NG0_PCC0.DRAGON_ISOTOPESDENS["Gd157"])
NGd157_D5_NG0_PCC0_interp = np.interp(VBOC_OMC["BU"], DRAGON_case_VBOC_NG0_PCC0.DRAGON_BU, DRAGON_case_VBOC_NG0_PCC0.DRAGON_ISOTOPESDENS["Gd157"])

# Compare NGd157

delta_NGd157_D5_NG0_S2_setQfiss = (NGd157_D5_NG0_interp - Serpent2_case_edep0_setQfiss.Ni["Gd157"]) * 100 / Serpent2_case_edep0_setQfiss.Ni["Gd157"]
delta_NGd157_D5_NG0_PCC0_S2_setQfiss = (NGd157_D5_NG0_PCC0_interp - Serpent2_case_edep0_setQfiss.Ni["Gd157"]) * 100 / Serpent2_case_edep0_setQfiss.Ni["Gd157"]

delta_NGd157_D5_NG0_OMC_setQfiss = (NGd157_D5_NG0_interp - OMC_case_edep0_setQfiss.Ni["Gd157"]) * 100 / OMC_case_edep0_setQfiss.Ni["Gd157"]
delta_NGd157_D5_NG0_PCC0_OMC_setQfiss = (NGd157_D5_NG0_PCC0_interp - OMC_case_edep0_setQfiss.Ni["Gd157"]) * 100 / OMC_case_edep0_setQfiss.Ni["Gd157"]

# plot delta keffs
plt.figure()
plt.plot(Serpent2_case_edep0_setQfiss.BU, delta_keff_D5_NG0_S2_setQfiss, label = f"D5 (NG0) - S2 setQfiss", marker = "D", linestyle = "--") # NG0 DEPL, EXTR - setQfiss, PCC0
plt.plot(OMC_case_edep0_setQfiss.BU, delta_keff_D5_NG0_OMC_setQfiss, label = f"D5 (NG0) - OMC setQfiss", marker = "D", linestyle = "--") # NG0 DEPL, EXTR - setQfiss, PCC0
plt.title("Keff differences between D5 and S2/OMC cases \n with NG0 on VBOC burnup points")
plt.xlabel("Burnup [MWd/tU]")
plt.ylabel("Keff difference [pcm]")
plt.legend()
plt.grid()
plt.savefig(f"{save_dir_comparisons_VBOC}/BestCandidates_Keff_diff_HOM_Gd157_VBOC_vs_S2_OMC_VBOC_OMC.png")
plt.close()

# plot delta NGd157
plt.figure()
plt.plot(Serpent2_case_edep0_setQfiss.BU, delta_NGd157_D5_NG0_S2_setQfiss, label = f"D5 (NG0) - S2 setQfiss", marker = "D", linestyle = "--") # NG0 DEPL, EXTR - setQfiss, PCC0
plt.plot(OMC_case_edep0_setQfiss.BU, delta_NGd157_D5_NG0_OMC_setQfiss, label = f"D5 (NG0) - OMC setQfiss", marker = "D", linestyle = "--") # NG0 DEPL, EXTR - setQfiss, PCC0
plt.title("Relative differences in Gd157 number density between D5 and S2/OMC cases \n with NG0 on VBOC burnup points")
plt.xlabel("Burnup [MWd/tU]")
plt.ylabel("Relative difference [%]")
plt.legend()
plt.grid()
plt.savefig(f"{save_dir_comparisons_VBOC}/BestCandidates_Gd157_Niso_rel_diff_HOM_Gd157_VBOC_vs_S2_OMC_VBOC_OMC.png")
plt.close()








