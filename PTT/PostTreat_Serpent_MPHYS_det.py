# Python3 script to post treat serpent2 rates output
# Author : R. Guasch
# Date : 11/10/2024
# Purpose : Post treat serpent2 rates output to validate power distribution in Multi-physics simulations

import numpy as np
import matplotlib.pyplot as plt
import serpentTools as st
import os
import sys


def parse_detector_rates(number_axial_slices, det):
    #print(det.detectors["_FUEL_1G_1"].tallies.shape)
    fiss_rates = {"U235":[], "U238":[], "U234":[]}
    n_gamma_rates = {"U235":[], "U238":[], "U234":[]}
    heat_production = {"U235":[], "U238":[], "U234":[]}
    for i in range(number_axial_slices):
        #detectors.append(det.detectors[f"_FUEL_1G_{i+1}"])
        n_gamma_rates["U235"].append(det.detectors[f"_FUEL_1G_{i+1}"].tallies[0])
        n_gamma_rates["U238"].append(det.detectors[f"_FUEL_1G_{i+1}"].tallies[1])
        n_gamma_rates["U234"].append(det.detectors[f"_FUEL_1G_{i+1}"].tallies[2])
        fiss_rates["U235"].append(det.detectors[f"_FUEL_1G_{i+1}"].tallies[3])
        fiss_rates["U238"].append(det.detectors[f"_FUEL_1G_{i+1}"].tallies[4])
        fiss_rates["U234"].append(det.detectors[f"_FUEL_1G_{i+1}"].tallies[5])
        heat_production["U235"].append(det.detectors[f"_FUEL_1G_{i+1}"].tallies[6])
        heat_production["U238"].append(det.detectors[f"_FUEL_1G_{i+1}"].tallies[7])
        heat_production["U234"].append(det.detectors[f"_FUEL_1G_{i+1}"].tallies[8])
    print(fiss_rates["U235"])           

    sum_fiss_rates = np.zeros(number_axial_slices)
    sum_n_gamma_rates = np.zeros(number_axial_slices)
    sum_heat_production = np.zeros(number_axial_slices)
    for i in range(number_axial_slices):
        for iso in fiss_rates.keys():
            sum_fiss_rates[i] += fiss_rates[iso][i]
            sum_n_gamma_rates[i] += n_gamma_rates[iso][i]
            sum_heat_production[i] += heat_production[iso][i]
    # normalize the fission rates
    fiss_rates = sum_fiss_rates/sum(sum_fiss_rates)
    n_gamma_rates = sum_n_gamma_rates/sum(sum_n_gamma_rates)
    heat_production = sum_heat_production/sum(sum_heat_production)

    return fiss_rates, n_gamma_rates, heat_production

def parse_detector_fluxes(number_axial_slices, det):

    flux_G1 = []
    flux_G2 = []

    for i in range(number_axial_slices):
        flux_G1.append(det.detectors[f"_FLUX_2G_{i+1}"].tallies[0])
        flux_G2.append(det.detectors[f"_FLUX_2G_{i+1}"].tallies[1])

    return flux_G1, flux_G2

def parse_DONJON_power():

    # retrieve DONJON5 power distribution
    donjon5_power = np.loadtxt(f"/home/p117902/working_dir/PolyBWR_project/Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/EPRIvoidModel_Churchill_HEM1/mesh20/Data/Power_Distrib_24UOX_mesh20_BiCG_EPRIvoidModel_relaxedPOW_0.1_relaxedTH_0.1.txt")

    donjon5_power2 = np.loadtxt(f"/home/p117902/working_dir/PolyBWR_project/Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/EPRIvoidModel_Churchill_HEM1/mesh20/Data/Power_Distrib_24UOX_mesh20_BiCG_EPRIvoidModel_relaxedPOW_0.5_relaxedTH_0.1.txt")

    donjon_power3 = np.loadtxt(f"/home/p117902/working_dir/PolyBWR_project/Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/EPRIvoidModel_Churchill_HEM1/mesh20/Data/Power_Distrib_24UOX_mesh20_BiCG_EPRIvoidModel_relaxedPOW_0.9_relaxedTH_0.1.txt")
    # normalize the donjon5 power distribution
    donjon5_power = donjon5_power/sum(donjon5_power)
    donjon5_power2 = donjon5_power2/sum(donjon5_power2)
    donjon_power3 = donjon_power3/sum(donjon_power3)
    
    return donjon5_power, donjon5_power2, donjon_power3

def parse_DONJON_fluxes():
    path_to_DISTR = "/home/p117902/working_dir/PolyBWR_project/Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/DISTR_res"
    # retrieve DONJON5 flux distribution
    donjon5_flux_G1 = []
    donjon5_flux_G2 = []
    with open(f"{path_to_DISTR}/Flux01.res", "r") as f:
        donjon5_flux_file1 = f.readlines()
        for line in donjon5_flux_file1:
            if "PLANE-Z" in line:
                plane_number = int(line.split("#")[-1])
                line_index = donjon5_flux_file1.index(line)
                donjon5_flux_G1.append(float(donjon5_flux_file1[line_index+2].strip()))
    with open(f"{path_to_DISTR}/Flux02.res", "r") as f:
        donjon5_flux_file2 = f.readlines()
        for line in donjon5_flux_file2:
            if "PLANE-Z" in line:
                plane_number = int(line.split("#")[-1])
                line_index = donjon5_flux_file2.index(line)
                donjon5_flux_G2.append(float(donjon5_flux_file2[line_index+2]))
    return donjon5_flux_G1, donjon5_flux_G2

def normalize_fluxes(flux_G1, flux_G2):
    # normalize the fluxes
    flux_G1 = np.array(flux_G1)/sum(np.array(flux_G1))
    flux_G2 = np.array(flux_G2)/sum(np.array(flux_G2))
    return flux_G1, flux_G2                

def compare_results(donjon5_powers, sum_fiss_rates):
    # compare the results
    fig3, ax3 = plt.subplots()
    ax3.plot(range(number_axial_slices), sum_fiss_rates, label="Fission rates")
    ax3.plot(range(number_axial_slices), donjon5_powers[0], label="Donjon5 power1")
    ax3.plot(range(number_axial_slices), donjon5_powers[1], label="Donjon5 power2")
    ax3.plot(range(number_axial_slices), donjon5_powers[2], label="Donjon power3")
    fig3.savefig(f"{case_name}_comparison.png")

def compute_relative_error(donjon5_rates, S2_rates):
    # calculate the relative error
    relative_error = np.abs(donjon5_rates-S2_rates)*100/S2_rates
    return relative_error

def plot_relative_error(relative_error, relative_error2, relative_error3):
    fig4, ax4 = plt.subplots()
    ax4.plot(range(number_axial_slices), relative_error, label="Relative error 0.1 relax")
    ax4.plot(range(number_axial_slices), relative_error2, label="Relative error 0.5 relax")
    ax4.plot(range(number_axial_slices), relative_error3, label="Relative error 0.9 relax")
    ax4.set_ylabel("Relative error (%)")
    ax4.set_xlabel("Axial slice")
    fig4.legend()
    fig4.savefig(f"{case_name}_relative_error.png")

    
case_name = "AT10_24UOX_3D_MPHYS_mc" # name of the serpent2 case
path_to_S2results = "/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/" # path to the serpent2 results
# Read the serpent output file
det = st.read(f'{path_to_S2results}/{case_name}_det0.m')
number_axial_slices = 20
n_mat_per_slice = 4
print(det.detectors)
#detectors = []
sum_fiss_rates, sum_n_gamma_rates, sum_heat_prod = parse_detector_rates(number_axial_slices, det)
donjon5_powers = parse_DONJON_power()
compare_results(donjon5_powers, sum_fiss_rates)
#relative_error, relative_error2, relative_error3 = compute_relative_error(donjon5_powers[0], donjon5_powers[1], donjon5_powers[2], sum_fiss_rates)
relative_error = compute_relative_error(donjon5_powers[0], sum_heat_prod)
relative_error2 = compute_relative_error(donjon5_powers[1], sum_n_gamma_rates)
relative_error3 = compute_relative_error(donjon5_powers[2], sum_fiss_rates)
# Plot the results
plot_relative_error(relative_error, relative_error2, relative_error3)


# begin fluxes post treatment
S2_flux_G1, S2_flux_G2 = parse_detector_fluxes(number_axial_slices, det)
donjon5_flux_G1, donjon5_flux_G2 = parse_DONJON_fluxes()

S2_flux_G1, S2_flux_G2 = normalize_fluxes(S2_flux_G1, S2_flux_G2)
donjon5_flux_G1, donjon5_flux_G2 = normalize_fluxes(donjon5_flux_G1, donjon5_flux_G2)
print(len(S2_flux_G1))
print(len(donjon5_flux_G1))

fig1, ax1 = plt.subplots()
ax1.plot(range(number_axial_slices), S2_flux_G1, label="S2 flux G1")
ax1.plot(range(number_axial_slices), donjon5_flux_G1, label="Donjon5 flux G1")
fig1.legend()
fig1.savefig(f"{case_name}_flux_G1.png")

fig2, ax2 = plt.subplots()
ax2.plot(range(number_axial_slices), S2_flux_G2, label="S2 flux G2")
ax2.plot(range(number_axial_slices), donjon5_flux_G2, label="Donjon5 flux G2")
fig2.legend()
fig2.savefig(f"{case_name}_flux_G2.png")

error_flux_G1 = compute_relative_error(donjon5_flux_G1, S2_flux_G1)
error_flux_G2 = compute_relative_error(donjon5_flux_G2, S2_flux_G2)

fig5, ax5 = plt.subplots()
ax5.plot(range(number_axial_slices), error_flux_G1, label="Relative error flux G1")
ax5.plot(range(number_axial_slices), error_flux_G2, label="Relative error flux G2")
fig5.legend()
fig5.savefig(f"{case_name}_relative_error_fluxes.png")
# end flux

error_G1G2 = compute_relative_error(donjon5_flux_G1, S2_flux_G2)
error_G2G1 = compute_relative_error(donjon5_flux_G2, S2_flux_G1)

fig6, ax6 = plt.subplots()
ax6.plot(range(number_axial_slices), error_G1G2, label="Relative error G1G2")
ax6.plot(range(number_axial_slices), error_G2G1, label="Relative error G2G1")
fig6.legend()
fig6.savefig(f"{case_name}_relative_error_G1G2_G2G1.png")

