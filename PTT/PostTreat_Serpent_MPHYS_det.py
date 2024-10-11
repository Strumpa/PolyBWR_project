# Python3 script to post treat serpent2 rates output
# Author : R. Guasch
# Date : 11/10/2024
# Purpose : Post treat serpent2 rates output to validate power distribution in Multi-physics simulations

import numpy as np
import matplotlib.pyplot as plt
import serpentTools as st
import os
import sys


case_name = "AT10_24UOX_3D_MPHYS_mc" # name of the serpent2 case
path_to_S2results = "/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/" # path to the serpent2 results
# Read the serpent output file
det = st.read(f'{path_to_S2results}/{case_name}_det0.m')
number_axial_slices = 20
n_mat_per_slice = 4
print(det.detectors)
#detectors = []
fiss_rates = {"U235":[], "U238":[], "U234":[]}
n_gamma_rates = {"U235":[], "U238":[], "U234":[]}
heat_production = {"U235":[], "U238":[], "U234":[]}
for i in range(number_axial_slices):
    #detectors.append(det.detectors[f"_FUEL_1G_{i+1}"])
    for mat_number in range(n_mat_per_slice):
        for iso_number in range(len(fiss_rates)):
            n_gamma_rates[list(fiss_rates.keys())[iso_number]].append(det.detectors[f"_FUEL_1G_{i+1}"].tallies[mat_number][iso_number])
            fiss_rates[list(fiss_rates.keys())[iso_number]].append(det.detectors[f"_FUEL_1G_{i+1}"].tallies[mat_number][iso_number+len(fiss_rates)])
        #print(det.detectors[f"_FUEL_1G_{i+1}"].tallies[0])
sum_fiss_rates = np.zeros(int(len(fiss_rates["U235"])/4)) 
sum_n_gamma_rates = np.zeros(int(len(n_gamma_rates["U235"])/4))           
for i in range(len(sum_fiss_rates)):
    # sum every 4 consecutive values to get the total fiss rate in each slice
    sum_fiss_rates[i]=sum(fiss_rates["U235"][i*4:i*4+4])+sum(fiss_rates["U238"][i*4:i*4+4])+sum(fiss_rates["U234"][i*4:i*4+4])   
    sum_n_gamma_rates[i]=sum(n_gamma_rates["U235"][i*4:i*4+4])+sum(n_gamma_rates["U238"][i*4:i*4+4])+sum(n_gamma_rates["U234"][i*4:i*4+4])


# normalize the fission rates
sum_fiss_rates = sum_fiss_rates/sum(sum_fiss_rates)
# Plot the results
fig, ax = plt.subplots()
ax.plot(range(number_axial_slices), sum_fiss_rates, label="Fission rates")
fig.savefig(f"{case_name}_fiss_rates.png")
# retrieve DONJON5 power distribution
donjon5_power = np.loadtxt(f"/home/p117902/working_dir/PolyBWR_project/Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/EPRIvoidModel_Churchill_HEM1/mesh20/Data/Power_Distrib_24UOX_mesh20_BiCG_EPRIvoidModel_relaxedPOW_0.1_relaxedTH_0.1.txt")

donjon5_power2 = np.loadtxt(f"/home/p117902/working_dir/PolyBWR_project/Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/EPRIvoidModel_Churchill_HEM1/mesh20/Data/Power_Distrib_24UOX_mesh20_BiCG_EPRIvoidModel_relaxedPOW_0.5_relaxedTH_0.1.txt")

donjon_power3 = np.loadtxt(f"/home/p117902/working_dir/PolyBWR_project/Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/EPRIvoidModel_Churchill_HEM1/mesh20/Data/Power_Distrib_24UOX_mesh20_BiCG_EPRIvoidModel_relaxedPOW_0.5_relaxedTH_0.1.txt")
# normalize the donjon5 power distribution
donjon5_power = donjon5_power/sum(donjon5_power)
donjon5_power2 = donjon5_power2/sum(donjon5_power2)
fig2, ax2 = plt.subplots()
ax2.plot(range(number_axial_slices), donjon5_power, label="Donjon5 power")
fig2.savefig(f"{case_name}_donjon5_power.png")

# compare the results
fig3, ax3 = plt.subplots()
ax3.plot(range(number_axial_slices), sum_fiss_rates, label="Fission rates")
ax3.plot(range(number_axial_slices), donjon5_power, label="Donjon5 power")
fig3.savefig(f"{case_name}_comparison.png")

# calculate the relative error
relative_error = np.abs(donjon5_power-sum_fiss_rates)*100/sum_fiss_rates
relative_error2 = np.abs(donjon5_power2-sum_fiss_rates)*100/sum_fiss_rates
fig4, ax4 = plt.subplots()
ax4.plot(range(number_axial_slices), relative_error, label="Relative error 0.1 relax")
ax4.plot(range(number_axial_slices), relative_error2, label="Relative error 0.5 relax")
ax4.set_ylabel("Relative error (%)")
ax4.set_xlabel("Axial slice")
fig4.legend()
fig4.savefig(f"{case_name}_relative_error.png")


