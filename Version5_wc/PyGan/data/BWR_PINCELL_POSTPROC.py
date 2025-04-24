### POST processing script for BWR pin cell

# Author: R. Guasch
# Date: 2025/04/14


# Importing libraries
import re
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import lcm

from postproc_cst_pow_evol import DRAGON_case as D5_case
from postproc_cst_pow_evol import Serpent2_case as S2_case
from getLists import getLists


# --- Selecting what case to post-treat
case_name = "AT10_24UOX_Cst_pow_evol"
cell_name = "AT10_24UOX"

ssh_method = "PT" # "RSE", "PT", "SUBG", "AUTO"
tracking_opt = "SALT" # "SALT", "SYBNXT"
# --- Selecting the evaluation used for the simulation
evaluation = "endfb8r1" #
# selecting solver / paramters used in EVO:
sol_depl = "RUNG"
dirac_opt = "NODI" # "DIRA", "NODI"
glob_opt = "GLOB"
# --- Selecting the burnup mesh
burnup_points = "UOx2_autop5" # "UOx", "UOx_autop5", "UOx2_autop5", "UOx4_autop5", "UOx6_autop5"

draglibs = ["endfb8r1_295", "endfb81295K", "endfb81295K2"]

#
path = os.getcwd()
save_dir_24UOX = f"{path}/BWR_PINCELL_POSTPROC/AT10_24UOX_evol/"
# create the directories if they do not exist
if not os.path.exists(save_dir_24UOX):
    os.makedirs(save_dir_24UOX)
path_to_PyGan_results = f"PYGAN_COMPOS_path/{case_name}_results/"
os.chdir(path_to_PyGan_results)
name_CPO_1 = f"CPO_{cell_name}_endfb8r1_295_{ssh_method}_{tracking_opt}_{sol_depl}_{dirac_opt}_{glob_opt}_{burnup_points}"
name_CPO_2 = f"CPO_{cell_name}_endfb81295K_{ssh_method}_{tracking_opt}_{sol_depl}_{dirac_opt}_{glob_opt}_{burnup_points}"
name_CPO_3 = f"CPO_{cell_name}_endfb81295K2_{ssh_method}_{tracking_opt}_{sol_depl}_{dirac_opt}_{glob_opt}_{burnup_points}"
name_CPO_4 = f"CPO_{cell_name}_EB81K2_fix_{ssh_method}_{tracking_opt}_{sol_depl}_{dirac_opt}_{glob_opt}_{burnup_points}"

# NOGL :
glob_opt = "NOGL"
name_CPO_1_NOGL = f"CPO_{cell_name}_endfb8r1_295_{ssh_method}_{tracking_opt}_{sol_depl}_{dirac_opt}_{glob_opt}_{burnup_points}"
name_CPO_2_NOGL = f"CPO_{cell_name}_endfb81295K_{ssh_method}_{tracking_opt}_{sol_depl}_{dirac_opt}_{glob_opt}_{burnup_points}"
#name_CPO_3_NOGL = f"CPO_{cell_name}_endfb81295K2_{ssh_method}_{tracking_opt}_{sol_depl}_{dirac_opt}_{glob_opt}_{burnup_points}"
name_CPO_4_NOGL = f"CPO_{cell_name}_EB81K2_fix_{ssh_method}_{tracking_opt}_{sol_depl}_{dirac_opt}_{glob_opt}_{burnup_points}"


# create LCM object
CPO_24UOX_NOKERMA = lcm.new('LCM_INP', name_CPO_1, impx=0)
CPO_24UOX_KERMAtot = lcm.new('LCM_INP', name_CPO_2, impx=0)
CPO_24UOX_KERMA2 = lcm.new('LCM_INP', name_CPO_3, impx=0)
CPO_24UOX_KERMA2_fix = lcm.new('LCM_INP', name_CPO_4, impx=0)

CPO_24UOX_NOKERMA_NOGL = lcm.new('LCM_INP', name_CPO_1_NOGL, impx=0)
CPO_24UOX_KERMAtot_NOGL = lcm.new('LCM_INP', name_CPO_2_NOGL, impx=0)
#CPO_24UOX_KERMA2_NOGL = lcm.new('LCM_INP', name_CPO_3_NOGL, impx=0)
CPO_24UOX_KERMA2_fix_NOGL = lcm.new('LCM_INP', name_CPO_4_NOGL, impx=0)
os.chdir(path)

# --- DRAGON5-SERPENT2 comparison : post treatment
# DRAGON case
# isotopes to be tracked 
tracked_nuclides = ["U235","U238","Pu239","Pu240","Pu241","Pu242","Am241","Xe135","Sm149"]
D5_AT10_24UOX_noKERMA = D5_case(pyCOMPO = CPO_24UOX_NOKERMA,
                            dlib_name = "endfb8r1_295",
                            bu_points = burnup_points,
                            ssh_opt = ssh_method,
                            correlation = "NOCORR",
                            sat = "",
                            depl_sol = "RUNG",
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists(burnup_points),
                            save_dir = save_dir_24UOX)
D5_AT10_24UOX_KERMAtot = D5_case(pyCOMPO = CPO_24UOX_KERMAtot,
                            dlib_name = "endfb8r1_295K",
                            bu_points = burnup_points,
                            ssh_opt = ssh_method,
                            correlation = "NOCORR",
                            sat = "",
                            depl_sol = "RUNG",
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists(burnup_points),
                            save_dir = save_dir_24UOX)
D5_AT10_24UOX_KERMA2 = D5_case(pyCOMPO = CPO_24UOX_KERMA2,
                            dlib_name = "endfb8r1_295K2",
                            bu_points = burnup_points,
                            ssh_opt = ssh_method,
                            correlation = "NOCORR",
                            sat = "",
                            depl_sol = "RUNG",
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists(burnup_points),
                            save_dir = save_dir_24UOX)
D5_AT10_24UOX_KERMA2_fix = D5_case(pyCOMPO = CPO_24UOX_KERMA2_fix,
                            dlib_name = "endfb8r1_295K2_fix",
                            bu_points = burnup_points,
                            ssh_opt = ssh_method,
                            correlation = "NOCORR",
                            sat = "",
                            depl_sol = "RUNG",
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists(burnup_points),
                            save_dir = save_dir_24UOX)

D5_AT10_24UOX_noKERMA_NOGL = D5_case(pyCOMPO = CPO_24UOX_NOKERMA_NOGL,
                            dlib_name = "endfb8r1_295",
                            bu_points = burnup_points,
                            ssh_opt = ssh_method,
                            correlation = "NOCORR",
                            sat = "",
                            depl_sol = "RUNG",
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists(burnup_points),
                            save_dir = save_dir_24UOX)

D5_AT10_24UOX_KERMAtot_NOGL = D5_case(pyCOMPO = CPO_24UOX_KERMAtot_NOGL,
                            dlib_name = "endfb8r1_295K",
                            bu_points = burnup_points,
                            ssh_opt = ssh_method,
                            correlation = "NOCORR",
                            sat = "",
                            depl_sol = "RUNG",
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists(burnup_points),
                            save_dir = save_dir_24UOX)
D5_AT10_24UOX_KERMA2_fix_NOGL = D5_case(pyCOMPO = CPO_24UOX_KERMA2_fix_NOGL,
                            dlib_name = "endfb8r1_295K2_fix",
                            bu_points = burnup_points,
                            ssh_opt = ssh_method,
                            correlation = "NOCORR",
                            sat = "",
                            depl_sol = "RUNG",
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists(burnup_points),
                            save_dir = save_dir_24UOX)





### Compare DRAGON5 Cases for 24UOX
# take no kerma as reference and compare with the other two
## Compute difference on keff [pcm]

delta_keff_24UOX_KERMAtot = (D5_AT10_24UOX_KERMAtot.keff - D5_AT10_24UOX_noKERMA.keff) * 1e5
delta_keff_24UOX_KERMA2 = (D5_AT10_24UOX_KERMA2.keff - D5_AT10_24UOX_noKERMA.keff) * 1e5
delta_keff_24UOX_KERMA2_fix = (D5_AT10_24UOX_KERMA2_fix.keff - D5_AT10_24UOX_noKERMA.keff) * 1e5

# plot delta_keffs

plt.figure(figsize=(10, 6))
plt.plot(D5_AT10_24UOX_noKERMA.BU, delta_keff_24UOX_KERMAtot, label='KERMAtot', marker='x')
plt.plot(D5_AT10_24UOX_noKERMA.BU, delta_keff_24UOX_KERMA2_fix, label='KERMA, MT301-MT318+MT458', marker='x')
plt.xlabel('Burnup (MWd/tU)')
plt.ylabel('Delta keff (pcm)')
plt.title('Delta keff for 24UOX cases')
plt.legend()
plt.grid()
plt.savefig(f'{save_dir_24UOX}/delta_keff_24UOX_KERMA_tot.png')
plt.close()


plt.figure(figsize=(10, 6))
#plt.plot(D5_AT10_24UOX_noKERMA.BU, delta_keff_24UOX_KERMAtot, label='KERMAtot', marker='x')
plt.plot(D5_AT10_24UOX_noKERMA.BU, delta_keff_24UOX_KERMA2, label='KERMA, MT301-MT318+MT458', marker='x')
plt.xlabel('Burnup (MWd/tU)')
plt.ylabel('Delta keff (pcm)')
plt.title('Delta keff for 24UOX cases')
plt.legend()
plt.grid()
plt.savefig(f'{save_dir_24UOX}/delta_keff_24UOX_KERMA2.png')
plt.close()


plt.figure(figsize=(10, 6))
plt.plot(D5_AT10_24UOX_noKERMA.BU, D5_AT10_24UOX_noKERMA.keff, label='No KERMA', marker='x')
plt.plot(D5_AT10_24UOX_KERMAtot.BU, D5_AT10_24UOX_KERMAtot.keff, label='KERMAtot', marker='x')
plt.plot(D5_AT10_24UOX_KERMA2.BU, D5_AT10_24UOX_KERMA2.keff, label='KERMA, MT301-MT318+MT458', marker='x')
plt.plot(D5_AT10_24UOX_KERMA2_fix.BU, D5_AT10_24UOX_KERMA2_fix.keff, label='KERMA2 fix', marker='x')
plt.xlabel('Burnup (MWd/tU)')
plt.ylabel('keff')
plt.title('keff for 24UOX cases')
plt.legend()
plt.grid()
plt.savefig(f'{save_dir_24UOX}/keff_24UOX_KERMA_test.png')
plt.close()

print(f"Initial keff for 24UOX no KERMA: {D5_AT10_24UOX_noKERMA.keff[0]}")
print(f"Initial keff for 24UOX KERMAtot: {D5_AT10_24UOX_KERMAtot.keff[0]}")
print(f"Initial keff for 24UOX KERMA2: {D5_AT10_24UOX_KERMA2.keff[0]}")

# --- SERPENT2 case

# case_name, lib_name, edep_id, areQfissSet, isEcaptSet, pcc_id, specific_power, tracked_nuclides, save_dir

S2_AT10_24UOX_edep2 = S2_case(case_name = "AT10_24UOX",
                            lib_name = "endfb8r1_pynjoy2012_kerma",
                            edep_id = 2,
                            areQfissSet = False,
                            isEcaptSet = False,
                            pcc_id = 1,
                            specific_power = 0.0,
                            tracked_nuclides = tracked_nuclides,
                            save_dir = save_dir_24UOX)

delta_keff_noKERMA = (D5_AT10_24UOX_noKERMA.keff - S2_AT10_24UOX_edep2.keff) * 1e5
delta_keff_KERMAtot = (D5_AT10_24UOX_KERMAtot.keff - S2_AT10_24UOX_edep2.keff) * 1e5
delta_keff_KERMA2 = (D5_AT10_24UOX_KERMA2_fix.keff - S2_AT10_24UOX_edep2.keff) * 1e5

delta_keff_noKERMA_NOGL = (D5_AT10_24UOX_noKERMA_NOGL.keff - S2_AT10_24UOX_edep2.keff) * 1e5
delta_keff_KERMAtot_NOGL = (D5_AT10_24UOX_KERMAtot_NOGL.keff - S2_AT10_24UOX_edep2.keff) * 1e5
delta_keff_KERMA2_NOGL = (D5_AT10_24UOX_KERMA2_fix_NOGL.keff - S2_AT10_24UOX_edep2.keff) * 1e5


# plot delta_keffs
plt.figure(figsize=(10, 6))
plt.plot(D5_AT10_24UOX_noKERMA.BU, delta_keff_noKERMA, label='No KERMA, GLOB', marker='x')
plt.plot(D5_AT10_24UOX_noKERMA.BU, delta_keff_KERMAtot, label='total KERMA, GLOB', marker='x')
plt.plot(D5_AT10_24UOX_noKERMA.BU, delta_keff_KERMA2, label='KERMA, MT301-MT318+MT458, GLOB', marker='x')
plt.xlabel('Burnup (MWd/tU)')
plt.ylabel('Delta keff (pcm)')
plt.title('Delta keff for 24UOX cases')
plt.legend()
plt.grid()
plt.savefig(f'{save_dir_24UOX}/delta_keff_24UOX_KERMA_S2.png')
plt.close()

plt.figure(figsize=(10, 6))
plt.plot(D5_AT10_24UOX_noKERMA.BU, delta_keff_noKERMA_NOGL, label='No KERMA, NOGL', marker='x')
plt.plot(D5_AT10_24UOX_noKERMA.BU, delta_keff_KERMAtot_NOGL, label='total KERMA, NOGL', marker='x')
plt.plot(D5_AT10_24UOX_noKERMA.BU, delta_keff_KERMA2_NOGL, label='KERMA, MT301-MT318+MT458, NOGL', marker='x')
plt.xlabel('Burnup (MWd/tU)')
plt.ylabel('Delta keff (pcm)')
plt.title('Delta keff for 24UOX cases, NOGL')
plt.legend()
plt.grid()
plt.savefig(f'{save_dir_24UOX}/delta_keff_24UOX_KERMA_S2_NOGL.png')
plt.close()


