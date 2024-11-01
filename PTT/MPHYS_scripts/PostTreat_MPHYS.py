# Python3 script to post treat serpent2 rates output
# Author : R. Guasch
# Date : 11/10/2024
# Purpose : Post treat serpent2 rates output to validate power distribution in Multi-physics simulations

import numpy as np
import matplotlib.pyplot as plt
import serpentTools as st
import os
import sys

### Parsing functions
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
    #print(fiss_rates["U235"])           

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

def parse_DONJON_power(height, mesh_size, pow_relax, power_scaling_factor, th_relax, voidFractionCorrel, frfaccorel, P2Pcorel, Ptot=40000):
    if height == 3.8:
        read_id = "h380"
    elif height == 1.555:
        read_id = "h1555"

    if pow_relax == False:
        power_relax_id = "non_relaxedPOW"
    else:
        power_relax_id = f"relaxedPow_{pow_relax}"
    if th_relax == False:
        th_relax_id = "non_relaxedTH"
    else:
        th_relax_id = f"relaxedTh_{th_relax}"
    # retrieve DONJON5 power distribution
    donjon5_power = np.loadtxt(f"/home/p117902/working_dir/PolyBWR_project/Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/{voidFractionCorrel}_{frfaccorel}_{P2Pcorel}/{read_id}/mesh{mesh_size}_{power_scaling_factor}/Data/Power_Distrib_24UOX_mesh{mesh_size}_BiCG_EPRIvoidModel_{power_relax_id}_{th_relax_id}.txt")
    # normalize the donjon5 power distribution
    donjon5_power = np.array(donjon5_power)/np.sum(donjon5_power)*Ptot # Renormalize to Ptot
    
    return donjon5_power

def parse_MPHYS_output_field(height, mesh_size, power_scaling_factor, field, pow_relax, th_relax, voidFractionCorrel,  frfaccorel, P2Pcorel):
    if height == 3.8:
        read_id = "h380"
    elif height == 1.555:
        read_id = "h1555"
    if pow_relax == False:
        power_relax_id = "non_relaxedPOW"
    else:
        power_relax_id = f"relaxedPOW_{pow_relax}"
    if th_relax == False:
        th_relax_id = "non_relaxedTH"
    else:
        th_relax_id = f"relaxedTH_{th_relax}"
    # retrieve DONJON5 power distribution
    field_values = np.loadtxt(f"/home/p117902/working_dir/PolyBWR_project/Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/{voidFractionCorrel}_{frfaccorel}_{P2Pcorel}/{read_id}/mesh{mesh_size}_{power_scaling_factor}/Data/{field}_24UOX_mesh{mesh_size}_BiCG_EPRIvoidModel_{power_relax_id}_{th_relax_id}.txt")
    return field_values

def parse_GenFoam_output_field(field, power_scaling_factor, voidFractionCorrel,  frfaccorel, P2Pcorel):
    """
    field = voidFraction, Pressure, T_liquid, T_vapour
    corresponding GeNFoam field : alpha.vapour / p / T.liquid / T.vapour
    voidFractionCorrel = "EPRIvoidModel"
    frfaccorel = "Churchill" / "blasius"
    P2Pcorel = "HEM1" / "lockhartMartinelli"
    """
    internal_field_values = []
    if field == "voidFraction":
        field = "alpha.vapour"
    elif field == "Pressure":
        field = "p"
    elif field == "T_liquid":
        field = "T.liquid"
    elif field == "T_vapour":
        field = "T.vapour"
    elif field == "U_liquid":
        field = "U.liquid"
    elif field == "U_vapour":
        field = "U.vapour"
    file_path = f"/home/p117902/working_dir/PolyBWR_project/PTT/MPHYS_scripts/GeNFoam_data/{voidFractionCorrel}_{frfaccorel}_{P2Pcorel}/fluidRegion_{power_scaling_factor}/{field}"
    with open(file_path, 'r') as file:
        lines = file.readlines()
    in_internal_field = False
    if field == "U.liquid" or field == "U.vapour":
        for i, line in enumerate(lines):
            # Locate the start of the internalField block
            
            if "internalField" in line and "nonuniform List<vector>" in line:
                text_to_parse = lines[i+3:i+78]
                #print(text_to_parse)
                #print(len(text_to_parse))
                for line in text_to_parse:
                    # Split line into values and strip any whitespace
                    value = line.strip().split()[0]
                    value = value.replace("(", "")
                    #print(value)
                    #for value in values:
                    internal_field_values.append(float(value))

    else:
        for i, line in enumerate(lines):
            # Locate the start of the internalField block
            if "internalField" in line and "nonuniform List<scalar>" in line:
                in_internal_field = True
                # The next line contains the number of values (skip it)
                i += 2  # Skip to the line with values in parentheses
                # Continue reading values in parentheses
                while '(' not in lines[i]:
                    i += 1
                i += 1  # Move to the first value line within parentheses
                
                # Collect values until closing parenthesis
                while ')' not in lines[i]:
                    # Split line into values and strip any whitespace
                    values = lines[i].strip().split()
                    if field == "U.liquid" or field == "U.vapour":
                        print(line)
                    internal_field_values.extend([float(value) for value in values])
                    #print(internal_field_values)
                    i += 1
                break

    return np.array(internal_field_values)

def parse_DONJON_Keffs(height, mesh_size, power_scaling_factor, pow_relax,th_relax, voidFractionCorrel,  frfaccorel, P2Pcorel):
    if height == 3.8:
        read_id = "h380"
    elif height == 1.555:
        read_id = "h1555"
    if pow_relax == False:
        power_relax_id = "non_relaxedPOW"
    else:
        power_relax_id = f"relaxedPow_{pow_relax}"
    if th_relax == False:
        th_relax_id = "non_relaxedTH"
    else:
        th_relax_id = f"relaxedTh_{th_relax}"
    # retrieve DONJON5 Keffs
    donjon5_keffs = np.loadtxt(f"/home/p117902/working_dir/PolyBWR_project/Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/{voidFractionCorrel}_{frfaccorel}_{P2Pcorel}/{read_id}/mesh{mesh_size}_{power_scaling_factor}/Data/Keffs_24UOX_mesh{mesh_size}_BiCG_EPRIvoidModel_{power_relax_id}_{th_relax_id}.txt")
    
    return donjon5_keffs

def parse_DONJON_fluxes(height, mesh_size, power_scaling_factor, voidFractionCorrel,  frfaccorel, P2Pcorel):
    if height == 3.8:
        read_id = "h380"
    elif height == 1.555:
        read_id = "h1555"
    path_to_DISTR = f"/home/p117902/working_dir/PolyBWR_project/Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/{voidFractionCorrel}_{frfaccorel}_{P2Pcorel}/{read_id}/mesh{mesh_size}_{power_scaling_factor}/Data"
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


### Normalization functions
def normalize_fluxes(flux_G1, flux_G2):
    # normalize the fluxes to Ptot
    flux_G1 = np.array(flux_G1)/sum(np.array(flux_G1))
    flux_G2 = np.array(flux_G2)/sum(np.array(flux_G2))
    return flux_G1, flux_G2             

def normalize_S2_power(Serpent2_power, Ptot):
    Serpent2_power = np.array(Serpent2_power)/sum(np.array(Serpent2_power))
    Serpent2_power = Serpent2_power*Ptot
    return Serpent2_power   


### Homogenization/interpolation functions

def average_neighbors(arr):
    # Ensure the input array has an even length
    if len(arr) % 2 != 0:
        raise ValueError("Input array length must be even.")
    
    # Sum every two neighboring elements
    result = np.add(arr[::2], arr[1::2])/2
    
    return result

def homogenize_field(fields_list):
    """
    Homogenize a field by averaging every two neighboring elements
    Comapre 160->80 to 80, 80->40 to 40, 40->20 to 20, 20->10 to 10
    """
    homogenized_fields = []
    most_converged = fields_list[-1]
    for i in range(len(fields_list)):
        homogenized_fields.append(average_neighbors(fields_list[i]))
    #homogenized_field = average_neighbors(field)
    return homogenized_fields

def condense_most_converged_field(most_converged_field):
    """
    condense the most converged field to compare to lower mesh sizes
    """
    most_conv_80 = average_neighbors(most_converged_field)
    most_conv_40 = average_neighbors(most_conv_80)
    most_conv_20 = average_neighbors(most_conv_40)
    most_conv_10 = average_neighbors(most_conv_20)
    return most_conv_10, most_conv_20, most_conv_40, most_conv_80 


# write interpolation function : take a list from larger mesh size and interpolate to get the values corresponding to a smaller mesh.
def interpolate_in_larger_mesh(height, smaller_mesh, larger_mesh_distr, interp_type="linear"):
    """
    length of a mesh give number of axial elements associated, so we can interpolate the values in the larger mesh at the smaller mesh's points
    """
    #print(f"smaller_mesh : {smaller_mesh}")
    n_larger_mesh = len(larger_mesh_distr)
    z_boundaries_larger_mesh = np.linspace(0, height, n_larger_mesh + 1)
    z_values = (z_boundaries_larger_mesh[:-1] + z_boundaries_larger_mesh[1:]) / 2  # Midpoints of control volumes
    # interpolate the values
    if interp_type == "linear":
        interpolated_values = np.interp(smaller_mesh, z_values, np.array(larger_mesh_distr))
    elif interp_type == "quadratic":
        poly_coeffs = np.polyfit(z_values, np.array(larger_mesh_distr), deg = 2)
        poly_interp = np.poly1d(poly_coeffs)
        interpolated_values = poly_interp(smaller_mesh)
    elif interp_type == "cubic":
        poly_coeffs = np.polyfit(z_values, np.array(larger_mesh_distr), deg = 3)
        poly_interp = np.poly1d(poly_coeffs)
        interpolated_values = poly_interp(smaller_mesh)
    return interpolated_values


### Error computing functions
def compute_relative_error(test_field, ref_field):
    # calculate the relative error
    relative_error = np.abs(test_field-ref_field)*100/ref_field
    return relative_error

def compute_RMS_error(field_test, field_ref, error_type):

    # calculate the quadratic error in units of field
    if error_type == "percent":
        RMS_error = np.sqrt(np.sum(((field_test-field_ref)*100/field_ref)**2)/len(field_ref))
    else:
        RMS_error = np.sqrt(np.sum(((field_test-field_ref))**2)/len(field_ref))
    return RMS_error

def error_Keffs_pcm(D5_keff, S2_keff):
    # compute delta_keff in pcm for converged Keff value
    delta_keff = (D5_keff - S2_keff)*1e5

    return delta_keff

def compute_nodal_error_to_most_converged(MPHYS_field_list, most_converged_field, interp_type, height):
    """
    compute the nodal error of each field to the most converged one
    """
    
    nodal_errors = []
    percent_nodal_errors = []
    for field in MPHYS_field_list:
        n = len(field)
        z_boundaries = np.linspace(0, height, n + 1)
        z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
        interpolated_field = interpolate_in_larger_mesh(height, z_values, most_converged_field, interp_type)
        nodal_error = np.abs(field - interpolated_field)
        nodal_errors.append(nodal_error)
        percent_nodal_error = nodal_error*100/interpolated_field
        percent_nodal_errors.append(percent_nodal_error)
        
    return nodal_errors, percent_nodal_errors

### Plotting functions
def compare_results(field, unit, donjon5_fields, S2_field, height, n_slices, name, save_name):
    if " " in field:
        save_field = field.replace(" ", "_")
    else:
        save_field = field
    if height == 3.8:
        save_id = "h380"
    elif height == 1.555:
        save_id = "h1555"
    # compute z_mesh
    z_boundaries = np.linspace(0, height, n_slices + 1)
    z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes

    # compare the results
    fig3, ax3 = plt.subplots()
    if S2_field is not None:
        ax3.plot(z_values, S2_field, label=f"S2 {field}", marker="x", linestyle="--", linewidth=0.5)
        S2_id = "S2"
    else:
        S2_id = ""
    for i in range(len(donjon5_fields)):
        ax3.plot(z_values, donjon5_fields[i], label=f"Donjon5 {field} {i}")
    
    ax3.set_title(f"Comparison of power distributions {name}")
    ax3.set_ylabel(f"{field} ({unit})")
    ax3.set_xlabel("Height (m)")
    ax3.grid()
    #fig3.legend(loc="best")
    fig3.savefig(f"Figures_{save_id}/AT10_24UOX_MPHYS_D5{S2_id}_{save_field}_comparison_{save_name}.png")

def plot_relative_error(relative_error_dict, height, type, power_scaling_factors_to_plot):
    if height == 3.8:
        save_id = "h380"
    elif height == 1.555:
        save_id = "h1555"
    fig, ax = plt.subplots()
    for key in relative_error_dict.keys():
        if "10" in key:
            # Compute the boundaries and midpoints of each control volume
            z_boundaries = np.linspace(0, height, 10 + 1)
        elif "20" in key:
            z_boundaries = np.linspace(0, height, 20 + 1)
        elif "40" in key:
            z_boundaries = np.linspace(0, height, 40 + 1)
        z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes in z
        power_scaling_factor = int(key.split(" ")[5])
        if power_scaling_factor in power_scaling_factors_to_plot:
            ax.plot(z_values, relative_error_dict[key], label=f"Relative error {key}", marker="x", linestyle="--", linewidth=0.5)
    ax.plot(z_values, 3.00*np.ones(len(z_values)), color="red", linestyle="--")
    ax.set_ylabel("Relative error (%)")
    ax.set_xlabel("Height m")
    fig.legend(loc="best")
    fig.savefig(f"Figures_{save_id}/AT10_24UOX_MPHYS_relative_errors_{type}.png")
    plt.close(fig)
    return

def plot_condensed_relative_error(relative_error_dict, height, type, power_scaling_factors_to_plot):
    if height == 3.8:
        save_id = "h380"
    elif height == 1.555:
        save_id = "h1555"
    fig, ax = plt.subplots()
    for key in relative_error_dict.keys():
        if "10" in key:
            # Compute the boundaries and midpoints of each control volume
            z_boundaries = np.linspace(0, height, 10 + 1)
        elif "20" in key:
            z_boundaries = np.linspace(0, height, 20 + 1)
        elif "40" in key:
            z_boundaries = np.linspace(0, height, 40 + 1)
        z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes in z
        power_scaling_factor = int(key.split(" ")[5])
        if power_scaling_factor in power_scaling_factors_to_plot:
            ax.plot(z_values, relative_error_dict[key], label=f"Relative error {key}", marker="x", linestyle="--", linewidth=0.5)
    ax.plot(z_values, 3.00*np.ones(len(z_values)), color="red", linestyle="--")
    ax.set_ylabel("Relative error (%)")
    ax.set_xlabel("Height m")
    fig.legend(loc="best")
    fig.savefig(f"Figures_{save_id}/AT10_24UOX_MPHYS_relative_errors_{type}.png")
    plt.close(fig)
    return

def plot_quadratic_errors(quadratic_error_dict, type, max_slices):
    """
    quadratic_error_dict : dictionary containing the quadratic errors for each mesh size
    type : string, type of the plot, flux or power
    max_slices : int, maximum number of axial slices
    """
    fig, ax = plt.subplots()
    z_slices = np.linspace(0, max_slices, 10)
    #print(f"z_slices : {z_slices}")
    for key in quadratic_error_dict.keys():
        
        num_axial_slices = key.split(" ")[3].split(",")[0]
        #print(f"num_axial_slices : {num_axial_slices}")
        print(f"key = {key}, quadratic_error_dict[key] : {quadratic_error_dict[key]}")
        ax.plot(float(num_axial_slices), quadratic_error_dict[key], label=f"{key}", marker="x", linestyle="--", linewidth=0.5)
    #ax5.plot(z_slices, 3.00*np.ones(len(z_slices)), color="red", linestyle="--")
    ax.set_ylabel("Quadratic error")
    ax.set_xlabel("Number of axial slices")
    #fig5.legend()
    fig.savefig(f"Figures_PT/AT10_24UOX_MPHYS_quadratic_errors_{type}.png")
    plt.close(fig)
    return

def plot_single_distr(title, field, unit, height,  n_slices, power_distr, name, save_name):
    if height == 3.8:
        save_id = "h380"
    elif height == 1.555:
        save_id = "h1555"
    # compute z_mesh
    z_boundaries = np.linspace(0, height, n_slices + 1)
    z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
    if " " in field:
        save_field = field.replace(" ", "_")
    else:
        save_field = field
    # plot distribution
    fig, ax = plt.subplots()
    ax.plot(z_values, power_distr, label=f"{field} distribution {name}", marker="x", linestyle="--", linewidth=0.5)
    ax.set_ylabel(f"{field} ({unit})")
    ax.set_xlabel("Height (m)")
    ax.set_title(title)
    ax.legend(loc="best")
    ax.grid()
   
    fig.savefig(f"Figures_{save_id}/AT10_24UOX_MPHYS_{save_field}_distribution_{save_name}.png")
    plt.close(fig)
    return

def plot_Keff_convergence(keffs, height, name, save_name):
    if height == 3.8:
        save_id = "h380"
    elif height == 1.555:
        save_id = "h1555"
    # plot keffs
    fig, ax = plt.subplots()
    ax.plot(keffs, label=f"Keffs {name}", marker="x", linestyle="--", linewidth=0.5)
    ax.set_ylabel("Keff")
    ax.set_xlabel("Number of iterations")
    ax.grid()
    fig.legend()
    fig.savefig(f"Figures_{save_id}/AT10_24UOX_MPHYS_keffs_conv_{save_name}.png")
    plt.close(fig)
    return

def plot_keff_error_spatial_convergence(delta_keffs_dict, height, name, save_name):
    """
    idea is to see how delta_Keff evolves with the number of axial slices increasing
    """
    if height == 3.8:
        save_id = "h380"
    elif height == 1.555:
        save_id = "h1555"
    fig, ax = plt.subplots()
    for key in delta_keffs_dict.keys():
        if "10" in key:
            n = 10
        elif "20" in key:
            n = 20
        elif "40" in key:
            n = 40
        elif "50" in key:
            n = 50
        elif "70" in key:
            n = 70
        elif "80" in key:
            n = 80
        elif "160" in key:
            n = 160
        ax.plot(n, delta_keffs_dict[key], marker="x", linestyle="--", linewidth=0.5)
    ax.set_ylabel("Delta Keff (pcm)")
    ax.set_xlabel("Number of axial slices")
    ax.set_title(f"Delta Keffs {name}")
    ax.grid()
    ax.legend(loc="best")
    fig.savefig(f"Figures_{save_id}/AT10_24UOX_MPHYS_delta_keffs_{save_name}.png")
    plt.close(fig)

    return

def plot_keff_D5S2_spatial_convergence(height, n_slices, D5_keffs, S2_keffs, name, save_name):
    """
    idea is to see how Keff evolves with the number of axial slices increasing : should reach spatial convergence !
    """
    if height == 3.8:
        save_id = "h380"
    elif height == 1.555:
        save_id = "h1555"
    fig, ax = plt.subplots()
    ax.plot(n_slices, D5_keffs, label=f"Keffs D5 {name}", marker="x", linestyle="--", linewidth=0.5)
    ax.plot(n_slices, S2_keffs, label=f"Keffs S2 {name}", marker="*", linestyle="--", linewidth=0.5)
    ax.set_ylabel("Keff")
    ax.set_xlabel("Number of axial slices")
    ax.set_title(f"Keffs {name}")
    ax.grid()
    ax.legend(loc="best")
    fig.savefig(f"Figures_{save_id}/AT10_24UOX_MPHYS_D5S2_keffs_{save_name}.png")
    plt.close(fig)

    return

def plot_nodal_errors_with_S2(height, n_slices_list, power_list, nodal_errors_dict, field, unit):
    """
    plot the nodal errors for each number of axial slices, error dict built from comparison with S2
    """
    if height == 3.8:
        save_id = "h380"
    elif height == 1.555:
        save_id = "h1555"
    
    fig, ax = plt.subplots()
    for n in n_slices_list:
        print(f"n = {n}")
        #z_meshes = []  
        z_boundaries = np.linspace(0, height, n + 1)
        z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
        #z_meshes.append(z_values)
        for power in power_list:
            #for i in range(len(n_slices_list)):
            key = f"Nodal error on power : {n} axial slices, {power}"
            print(f"key : {key}")
            print(f"nodal_errors_dict[key] : {nodal_errors_dict[key]}")
            ax.plot(z_values, nodal_errors_dict[key], label=f"{key}", marker="x", linestyle="--", linewidth=0.5)
    ax.plot(z_values, 5.00*np.ones(len(z_values)), color="red", linestyle="--")
    ax.plot(z_values, -5.00*np.ones(len(z_values)), color="red", linestyle="--")
    ax.set_ylabel(f"Nodal error on {field} {unit}")
    ax.set_xlabel("Height (m)")
    ax.grid()
    ax.legend(loc="best")
    fig.savefig(f"Figures_{save_id}/AT10_24UOX_MPHYS_nodal_errors_{field}.png")
    plt.close(fig)
    return

def plot_nodal_spatial_errors(height, field, unit, nodal_spatial_errors, interp_type):
    """
    plot the nodal errors comapring spatially converged fields
    """
    if height == 3.8:
        save_id = "h380"
    elif height == 1.555:
        save_id = "h1555"

    if interp_type == "linear":
        directory = "linear_interp"
    elif interp_type == "quadratic":
        directory = "quadratic_interp"
    elif interp_type == "cubic":
        directory = "cubic_interp"
    elif interp_type == "condensed":
        directory = "condensed"
    if " " in field:
        save_field = field.replace(" ", "_")
    else:
        save_field = field
    fig, ax = plt.subplots()
    for nodal_errors in nodal_spatial_errors:
        n=len(nodal_errors)
        # compute z_mesh
        z_boundaries = np.linspace(0, height, n + 1)
        z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2
        ax.plot(z_values, nodal_errors, label=f"nodal error on {field}, n={n}", marker="x", linestyle="--", linewidth=0.5)
    if unit != "":
        ax.set_ylabel(f"Nodal error on {field} ({unit})")
    else:
        ax.set_ylabel(f"Nodal error on {field}")
    ax.set_xlabel("Height (m)")
    ax.set_title(f"Nodal errors on {field}, compared to n=160 axial mesh")
    ax.grid()
    ax.legend(loc="best")
    if unit == "%":
        fig.savefig(f"Figures_{save_id}/{directory}/MPHYS_{save_field}_nodal_spatial_percent_errors_{interp_type}_on_mesh160.png")
    else:
        fig.savefig(f"Figures_{save_id}/{directory}/MPHYS_{save_field}_nodal_spatial_errors_{interp_type}_on_mesh160.png")

    plt.close(fig)
    return

def plot_RMS_errors(field, unit, height, n_slices_list, quadratic_error_dict):
    """
    plot the quadratic errors for each number of axial slices
    """
    if height == 3.8:
        save_id = "h380"
    elif height == 1.555:
        save_id = "h1555"
    if " " in field:
        save_field = field.replace(" ", "_")
    else:
        save_field = field
    RMS_errors_pow1 = []
    RMS_errors_pow2 = []
    fig, ax = plt.subplots()
    for n in n_slices_list:
       
        key_40 = f"RMS error on power : {n} axial slices, 40kW"
        key_20 = f"RMS error on power : {n} axial slices, 20kW"
        RMS_errors_pow1.append(quadratic_error_dict[key_40])
        RMS_errors_pow2.append(quadratic_error_dict[key_20])
    
    ax.plot(n_slices_list, np.array(RMS_errors_pow1), label=f"RMS error on power density, Ptot = 40kW", marker="x", linestyle="--", linewidth=0.8)
    ax.plot(n_slices_list, np.array(RMS_errors_pow2), label=f"RMS error on power density, Ptot = 20kW", marker="D", linestyle="--", linewidth=0.8)
    ax.plot(n_slices_list, 3.00*np.ones(len(n_slices_list)), color="red", linestyle="--")
    ax.set_ylabel(f"RMS error on {field} ({unit})")
    ax.set_xlabel("Number of axial nodes")
    ax.grid()
    ax.legend(loc="best")
    if unit == "%":
        fig.savefig(f"Figures_{save_id}/MPHYS_D5S2_{save_field}_quadratic_errors_percent.png")
    else:
        fig.savefig(f"Figures_{save_id}/MPHYS_D5S2_{save_field}_quadratic_errors.png")
    plt.close(fig)
    return

def compare_all_numbers_slices(field, unit, donjon_field_list, serpent_field_list, n_slices_list):
    """
    Compare the distributions for all numbers of axial slices
    """
    if height == 3.8:
        save_id = "h380"
    elif height == 1.555:
        save_id = "h1555"   

    if " " in field:
        save_field = field.replace(" ", "_")
    else:
        save_field = field
    fig, ax = plt.subplots()
    for i in range(len(n_slices_list)):
        ax.plot(donjon_field_list[i], label=f"Donjon5 {field} {n_slices_list[i]} axial slices")
        if serpent_field_list is not None:
            ax.plot(serpent_field_list[i], label=f"S2 {field} {n_slices_list[i]} axial slices")
    ax.set_ylabel(f"{field} ({unit})")
    ax.set_xlabel("Height (m)")
    ax.grid()
    ax.legend(loc="best")
    fig.savefig(f"Figures_{save_id}/AT10_24UOX_MPHYS_D5S2_comparison_{save_field}_all_slices.png")
    plt.close(fig)
    return

def plot_spatial_convergence_field(field, unit, fields_list, height, n_slices_list):
    """
    plot the spatial convergence of a field
    field = string, name of the field
    unit = string, unit of the field
    field_list = list of fields for each number of axial slices
    n_slices_list = list of number of axial slices
    """
    if height == 3.8:
        save_id = "h380"
    elif height == 1.555:
        save_id = "h1555"

    if " " in field:
        save_field = field.replace(" ", "_")
    else:
        save_field = field
    fig, ax = plt.subplots()
    for n in range(len(n_slices_list)):
        # compute z_mesh
        z_boundaries = np.linspace(0, height, n_slices_list[n] + 1)
        z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
        if n == len(n_slices_list)-1:
            ax.plot(z_values, fields_list[n], label=f"{field} {n_slices_list[n]} axial nodes", marker=".", linestyle="--", linewidth=0.2)
        else:
            ax.plot(z_values, fields_list[n], label=f"{field} {n_slices_list[n]} axial nodes", marker="x", linestyle="-.", linewidth=0.3)
    ax.set_ylabel(f"{field} ({unit})")
    ax.set_xlabel("Height (m)")
    ax.set_title(f"Spatial convergence of {field}")
    ax.grid()
    ax.legend(loc="best")
    fig.savefig(f"Figures_{save_id}/MPHYS_{save_field}_spatial_convergence.png")
    plt.close(fig)
        
def plot_piecewise_constant_field(field, unit, field_list, height, n_slices_list):
    """
    plot the piecewise constant field
    field = string, name of the field
    unit = string, unit of the field
    field_list = list of fields for each number of axial slices
    n_slices_list = list of number of axial slices
    """
    if height == 3.8:
        save_id = "h380"
    elif height == 1.555:
        save_id = "h1555"

    fig, ax = plt.subplots()
    for n in range(len(n_slices_list)):
        # compute z_mesh
        z_boundaries = np.linspace(0, height, n_slices_list[n] + 1)
        #z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
        # Generate the x and y values for plotting
        z_vals = []
        y_vals = []
        intervals = [(z_boundaries[i], z_boundaries[i+1]) for i in range(len(z_boundaries)-1)]
        for (start, end), value in zip(intervals, field_list[n]):
            z_vals.extend([start, end])  # Extend x values to cover the interval start and end
            y_vals.extend([value, value]) # Extend y values to cover the interval value
        
        ax.step(z_vals, y_vals, label=f"{field} : {n_slices_list[n]} axial nodes")
    if unit != "":
        ax.set_ylabel(f"{field} ({unit}) ")
    else:
        ax.set_ylabel(f"{field}")
    ax.set_xlabel("Height (m)")
    ax.set_title(f"{field} evolution for different axial meshes")
    ax.grid()
    ax.legend(loc="best")
    fig.savefig(f"Figures_{save_id}/MPHYS_{field}_piecewise_constant.png")
    plt.close(fig)

def plot_GeNFoam_vs_MPHYS(field, unit, GF_field, MPHYS_field, height, n_axial_nodes):
    """
    plot the piecewise constant field
    field = string, name of the field
    unit = string, unit of the field
    field_list = list of fields for each number of axial slices
    n_slices_list = list of number of axial slices
    """
    #fields_list = [GF_field, MPHYS_field]
    fig, ax = plt.subplots()
    # compute z_mesh
    z_boundaries = np.linspace(0, height, n_axial_nodes + 1)
    #z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
    # Generate the x and y values for plotting
    z_vals = []
    y_vals = []
    intervals = [(z_boundaries[i], z_boundaries[i+1]) for i in range(len(z_boundaries)-1)]
    #for n in range(len(fields_list)):
    for (start, end), value in zip(intervals, GF_field):
        z_vals.extend([start, end])  # Extend x values to cover the interval start and end
        y_vals.extend([value, value]) # Extend y values to cover the interval value
    ax.step(z_vals, y_vals, label = f"{field} : GeN-Foam")
    z_vals = []
    y_vals = []
    for (start, end), value in zip(intervals, MPHYS_field):
        z_vals.extend([start, end])  # Extend x values to cover the interval start and end
        y_vals.extend([value, value]) # Extend y values to cover the interval value
    ax.step(z_vals, y_vals, label = f"{field} : THM prototype")
    if unit != "":
        ax.set_ylabel(f"{field} ({unit}) ")
    else:
        ax.set_ylabel(f"{field}")
    ax.set_xlabel("Height (m)")
    ax.set_title(f"{field} : GeN-Foam vs THM prototype")
    ax.grid()
    ax.legend(loc="best")
    fig.savefig(f"Figures_GeNFoam_THM/{field}_piecewise_constant.png")
    plt.close(fig)

def plot_GeNFoam_vs_MPHYS_vl(field, unit, GF_L_field, GF_V_field, MPHYS_field, height, n_axial_nodes):
    """
    Variant to plot GeN-Foam vs THM prototype fields with GeNFoam vapour and liquid phases fields
    plot the piecewise constant field
    field = string, name of the field
    unit = string, unit of the field
    field_list = list of fields for each number of axial slices
    n_slices_list = list of number of axial slices
    """
    #fields_list = [GF_field, MPHYS_field]
    fig, ax = plt.subplots()
    # compute z_mesh
    z_boundaries = np.linspace(0, height, n_axial_nodes + 1)
    z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
    # Generate the x and y values for plotting

    ax.plot(z_values, MPHYS_field, label = f"{field} : THM prototype, mixture", marker="x", linestyle="--", linewidth=0.5)
    ax.plot(z_values, GF_L_field, label = f"{field} : GeN-Foam, liquid phase", marker="x", linestyle="-.", linewidth=0.5)
    ax.plot(z_values, GF_V_field, label = f"{field} : GeN-Foam, vapour phase", marker="x", linestyle=":", linewidth=0.5)
    if unit != "":
        ax.set_ylabel(f"{field} ({unit}) ")
    else:
        ax.set_ylabel(f"{field}")
    ax.set_xlabel("Height (m)")
    ax.set_title(f"{field} : GeN-Foam vs THM prototype")
    ax.grid()
    ax.legend(loc="best")
    fig.savefig(f"Figures_GeNFoam_THM/{field}_dots.png")
    plt.close(fig)

    # plot the piecewise constant field
    fig, ax = plt.subplots()
    # compute z_mesh
    z_boundaries = np.linspace(0, height, n_axial_nodes + 1)
    z_vals = []
    y_vals = []
    intervals = [(z_boundaries[i], z_boundaries[i+1]) for i in range(len(z_boundaries)-1)]
    #for n in range(len(fields_list)):
    for (start, end), value in zip(intervals, GF_L_field):
        z_vals.extend([start, end])  # Extend x values to cover the interval start and end
        y_vals.extend([value, value]) # Extend y values to cover the interval value
    ax.step(z_vals, y_vals, label = f"{field} : GeN-Foam, liquid phase")
    z_vals = []
    y_vals = []
    for (start, end), value in zip(intervals, GF_V_field):
        z_vals.extend([start, end])  # Extend x values to cover the interval start and end
        y_vals.extend([value, value]) # Extend y values to cover the interval value
    ax.step(z_vals, y_vals, label = f"{field} : GeN-Foam, vapor phase")
    z_vals = []
    y_vals = []
    for (start, end), value in zip(intervals, MPHYS_field):
        z_vals.extend([start, end])  # Extend x values to cover the interval start and end
        y_vals.extend([value, value]) # Extend y values to cover the interval value
    ax.step(z_vals, y_vals, label = f"{field} : THM prototype")
    if unit != "":
        ax.set_ylabel(f"{field} ({unit}) ")
    else:
        ax.set_ylabel(f"{field}")
    ax.set_xlabel("Height (m)")
    ax.set_title(f"{field} : GeN-Foam vs THM prototype")
    ax.grid()
    ax.legend(loc="best")
    fig.savefig(f"Figures_GeNFoam_THM/{field}_piecewise_constant.png")
    plt.close(fig)

def plot_nodal_errors_with_GeNFoam(height, nodal_errors_P, nodal_errors_VF, nodal_errors_TL, nodal_errors_TV, nodal_errors_UL, nodal_errors_UV):
    """
    plot the nodal errors for each field
    """
    fig, ax = plt.subplots()
    # compute z_mesh
    z_boundaries = np.linspace(0, height, len(nodal_errors_P) + 1)
    z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
    ax.plot(z_values, nodal_errors_P, label="Nodal error on Pressure", marker="x", linestyle="--", linewidth=0.5)
    ax.plot(z_values, nodal_errors_VF, label="Nodal error on Void Fraction", marker="x", linestyle="--", linewidth=0.5)
    ax.plot(z_values, nodal_errors_TL, label="Nodal error on T.liquid", marker="x", linestyle="--", linewidth=0.5)
    ax.plot(z_values, nodal_errors_TV, label="Nodal error on T.vapor", marker="x", linestyle="--", linewidth=0.5)
    #ax.plot(z_values, nodal_errors_UL, label="Nodal error on U.liquid", marker="x", linestyle="--", linewidth=0.5)
    #ax.plot(z_values, nodal_errors_UV, label="Nodal error on U.vapor", marker="x", linestyle="--", linewidth=0.5)
    ax.set_ylabel("Nodal error (%)")
    ax.set_xlabel("Height (m)")
    ax.grid()
    ax.legend(loc="best")
    fig.savefig(f"Figures_GeNFoam_THM/GeNFoam_nodal_errors.png")
    plt.close(fig)
    return


### LaTeX table printing functions time !!!
def print_latex_table_from_RMS_errors_dict(n_slices_list, quadratic_error_dict):
    """
    print a latex table with the quadratic errors for each number of axial slices
    """
    print("\\begin{table}[H]")
    print("\\centering")
    print("\\begin{tabular}{|c|c|}")
    print("\\hline")
    print("Number of axial slices & RMS error on power \\\\")
    print("\\hline")
    for i in range(len(n_slices_list)):
        key = f"RMS error on power : {n_slices_list[i]} axial slices"
        print(f"{n_slices_list[i]} & {quadratic_error_dict[key]} \\\\")
    print("\\hline")
    print("\\end{tabular}")
    print("\\end{table}")
    return

def print_latex_table_from_RMS_errors_list(field, unit_of_error, n_slices_list, RMS_error_list):
    """
    print a latex table with the quadratic errors for each number of axial slices
    """
    print("\\begin{table}[H]")
    print("\\centering")
    print("\\begin{tabular}{|c|c|}")
    print("\\hline")
    print(f"Number of axial slices & RMS error on {field} ({unit_of_error}) \\\\")
    print("\\hline")
    for i in range(len(n_slices_list)):
        print(f"{n_slices_list[i]} & {RMS_error_list[i]} \\\\")
    print("\\hline")
    print("\\end{tabular}")
    print("\\end{table}")
    return

# Post treatment and checker functions

def post_treat_Fields_SpatialConvergence_interp(height, power_scaling_factor, number_axial_slices, pow_relax, th_relax, voidFractionCorrel,  frfaccorel, P2Pcorel):
    """
    power scaling factor : int, power scaling factor for the DONJON5 simulations
    number_axial_slices : list, number of axial slices for the simulations
    pow_relax : float/bool, relaxation parameter for the power distribution
    th_relax : float/bool, relaxation parameter for the thermal-hydraulic parameters
    check spatial convergence of all fields from the MPHYS output
    1) Power density distribution for all meshes
    2) Fuel Temperature distribution for all meshes
    3) Water Temperature distribution for all meshes
    4) Water density distribution for all meshes
    5) Void fraction distribution for all meshes
    """
    print(f"$$$ - BEGIN post treatment of fields spatial convergence for {height} m height, {power_scaling_factor} power scaling factor")
    print(f"The number of axial slices are : {number_axial_slices}")
    print(f"The power and TH relaxation parameters are : {pow_relax}, {th_relax}")
    print(f"The void fraction, frfaccorel and P2Pcorel correlations are : {voidFractionCorrel}, {frfaccorel}, {P2Pcorel}")

    Qfiss_all_meshes = [] 
    TF_all_meshes = []
    TC_all_meshes = []
    DC_all_meshes = []
    VF_all_meshes = []
    for n in number_axial_slices:
        Qfiss_all_meshes.append(parse_MPHYS_output_field(height, n, power_scaling_factor, "Qfiss", pow_relax, th_relax, voidFractionCorrel,  frfaccorel, P2Pcorel))
        TF_all_meshes.append(parse_MPHYS_output_field(height, n, power_scaling_factor, "TeffFuel", pow_relax, th_relax, voidFractionCorrel,  frfaccorel, P2Pcorel))
        TC_all_meshes.append(parse_MPHYS_output_field(height, n, power_scaling_factor, "Twater", pow_relax, th_relax, voidFractionCorrel,  frfaccorel, P2Pcorel))
        DC_all_meshes.append(parse_MPHYS_output_field(height, n, power_scaling_factor, "rho", pow_relax, th_relax, voidFractionCorrel,  frfaccorel, P2Pcorel))
        VF_all_meshes.append(parse_MPHYS_output_field(height, n, power_scaling_factor, "voidFraction", pow_relax, th_relax, voidFractionCorrel,  frfaccorel, P2Pcorel))

    # plot the spatial convergence of the fields
    plot_spatial_convergence_field("power density distribution", "W/$m^3$", Qfiss_all_meshes, height, number_axial_slices)
    plot_piecewise_constant_field("Power density", "W/$m^3$", Qfiss_all_meshes, height, number_axial_slices)
    plot_spatial_convergence_field("fuel temperature", "K", TF_all_meshes, height, number_axial_slices)
    plot_piecewise_constant_field("Fuel temperature", "K", TF_all_meshes, height, number_axial_slices)
    plot_spatial_convergence_field("coolant temperature", "K", TC_all_meshes, height, number_axial_slices)
    plot_piecewise_constant_field("Coolant temperature", "K", TC_all_meshes, height, number_axial_slices)
    plot_spatial_convergence_field("coolant density", "kg/$m^3$", DC_all_meshes, height, number_axial_slices)
    plot_piecewise_constant_field("Coolant density", "kg/$m^3$", DC_all_meshes, height, number_axial_slices)
    plot_spatial_convergence_field("void fraction", "", VF_all_meshes, height, number_axial_slices)
    plot_piecewise_constant_field("Void fraction", "", VF_all_meshes, height, number_axial_slices)

    # comapre all distribution to the most spatially converged one, consider different interpolation types
    interp_type ="linear"
    spatial_errors_Qfiss, percent_spatial_errors_Qfiss = compute_nodal_error_to_most_converged(Qfiss_all_meshes[:-1], Qfiss_all_meshes[-1], interp_type, height)
    #print(f"spatial_errors_Qfiss : {spatial_errors_Qfiss}")
    spatial_errors_TF, percent_spatial_errors_TF = compute_nodal_error_to_most_converged(TF_all_meshes[:-1], TF_all_meshes[-1], interp_type, height)
    spatial_errors_TC, percent_spatial_errors_TC = compute_nodal_error_to_most_converged(TC_all_meshes[:-1], TC_all_meshes[-1], interp_type, height)
    spatial_errors_DC, percent_spatial_errors_DC = compute_nodal_error_to_most_converged(DC_all_meshes[:-1], DC_all_meshes[-1], interp_type, height)
    spatial_errors_VF, percent_spatial_errors_VF = compute_nodal_error_to_most_converged(VF_all_meshes[:-1], VF_all_meshes[-1], interp_type, height)

    # plot the nodal spatial errors for each field
    plot_nodal_spatial_errors(height, "power density distribution", "W/$m^3$", spatial_errors_Qfiss, interp_type) 
    plot_nodal_spatial_errors(height, "power density distribution", "%", percent_spatial_errors_Qfiss, interp_type)
    plot_nodal_spatial_errors(height, "fuel temperature", "K", spatial_errors_TF, interp_type)
    plot_nodal_spatial_errors(height, "fuel temperature", "%", percent_spatial_errors_TF, interp_type)
    plot_nodal_spatial_errors(height, "coolant temperature", "K", spatial_errors_TC, interp_type)
    plot_nodal_spatial_errors(height, "coolant temperature", "%", percent_spatial_errors_TC, interp_type)
    plot_nodal_spatial_errors(height, "coolant density", "kg/$m^3$", spatial_errors_DC, interp_type)
    plot_nodal_spatial_errors(height, "coolant density", "%", percent_spatial_errors_DC, interp_type)
    plot_nodal_spatial_errors(height, "void fraction", "", spatial_errors_VF, interp_type)
    plot_nodal_spatial_errors(height, "void fraction", "%", percent_spatial_errors_VF, interp_type)

    # compute RMS errors (%) on the fields vs the most spatially converged one

    RMS_errors_Qfiss = []
    RMS_errors_TF = []
    RMS_errors_TC = []
    RMS_errors_DC = []
    RMS_errors_VF = []

    for error in percent_spatial_errors_Qfiss:
        RMS_errors_Qfiss.append(np.sqrt(np.sum(error**2)/len(error)))
    for error in spatial_errors_TF:
        RMS_errors_TF.append(np.sqrt(np.sum(error**2)/len(error)))
    for error in spatial_errors_TC:
        RMS_errors_TC.append(np.sqrt(np.sum(error**2)/len(error)))
    for error in spatial_errors_DC:
        RMS_errors_DC.append(np.sqrt(np.sum(error**2)/len(error)))
    for error in spatial_errors_VF:
        RMS_errors_VF.append(np.sqrt(np.sum(error**2)/len(error)))


    # plot the RMS errors for each field
        
    print_latex_table_from_RMS_errors_list("power density distribution", "$%$", number_axial_slices[:-1], RMS_errors_Qfiss)
    print_latex_table_from_RMS_errors_list("fuel temperature", "K", number_axial_slices[:-1], RMS_errors_TF)
    print_latex_table_from_RMS_errors_list("coolant temperature", "K", number_axial_slices[:-1], RMS_errors_TC)
    print_latex_table_from_RMS_errors_list("coolant density", "$kg/m^3%", number_axial_slices[:-1], RMS_errors_DC)
    print_latex_table_from_RMS_errors_list("void fraction", "%", number_axial_slices[:-1], RMS_errors_VF)
    
    print(f"$$$ - END post treatment of fields spatial convergence for {height} m height, {power_scaling_factor} power scaling factor")

    return

def post_treat_Fields_spatial_convergence_condense(height, power_scaling_factor, number_axial_slices, pow_relax, th_relax, voidFractionCorrel,  frfaccorel, P2Pcorel):
    # Interpolation doesnt really work all that well due to the sharp gradients in the fields
    Qfiss_all_meshes = [] 
    TF_all_meshes = []
    TC_all_meshes = []
    DC_all_meshes = []
    VF_all_meshes = []
    for n in number_axial_slices:
        Qfiss_all_meshes.append(parse_MPHYS_output_field(height, n, power_scaling_factor, "Qfiss", pow_relax, th_relax, voidFractionCorrel,  frfaccorel, P2Pcorel))
        TF_all_meshes.append(parse_MPHYS_output_field(height, n, power_scaling_factor, "TeffFuel", pow_relax, th_relax, voidFractionCorrel,  frfaccorel, P2Pcorel))
        TC_all_meshes.append(parse_MPHYS_output_field(height, n, power_scaling_factor, "Twater", pow_relax, th_relax, voidFractionCorrel,  frfaccorel, P2Pcorel))
        DC_all_meshes.append(parse_MPHYS_output_field(height, n, power_scaling_factor, "rho", pow_relax, th_relax, voidFractionCorrel,  frfaccorel, P2Pcorel))
        VF_all_meshes.append(parse_MPHYS_output_field(height, n, power_scaling_factor, "voidFraction", pow_relax, th_relax, voidFractionCorrel,  frfaccorel, P2Pcorel))

    n_slices_to_plot = [40, 80, 160]
    Qfiss_to_plot = [Qfiss_all_meshes[2], Qfiss_all_meshes[3], Qfiss_all_meshes[4]]
    #plot_spatial_convergence_field("power density distribution special plot", "W/m3", Qfiss_to_plot, n_slices_to_plot)
    #plot_piecewise_constant_field("power density distribution", "W/m3", Qfiss_to_plot, n_slices_to_plot)
    

    # attempt averaging most converged fields to compare with lower mesh sizes

    most_converged_Qfiss = Qfiss_all_meshes[-1]
    most_converged_TF = TF_all_meshes[-1]
    most_converged_TC = TC_all_meshes[-1]
    most_converged_DC = DC_all_meshes[-1]
    most_converged_VF = VF_all_meshes[-1]

    most_conv_Qfiss_condensed = condense_most_converged_field(most_converged_Qfiss)
    most_conv_TF_condensed = condense_most_converged_field(most_converged_TF)
    most_conv_TC_condensed = condense_most_converged_field(most_converged_TC)
    most_conv_DC_condensed = condense_most_converged_field(most_converged_DC)
    most_conv_VF_condensed = condense_most_converged_field(most_converged_VF)

    # compute the relative errors for each field

    rel_errors_Qfiss = []
    abs_errors_Qfiss = []
    rel_errors_TF = []
    abs_errors_TF = []
    rel_errors_TC = []
    abs_errors_TC = []
    rel_errors_DC = []
    abs_errors_DC = []
    rel_errors_VF = []
    abs_errors_VF = []

    for i in range(len(most_conv_Qfiss_condensed)):
        if len(Qfiss_all_meshes[i]) == 70:
            print("skipping 70 axial nodes mesh")
        else:
            print(f"Computing relative/absolute errors for raw data on mesh {number_axial_slices[i]} to most converged mesh, condensed")
            rel_errors_Qfiss.append(compute_relative_error(Qfiss_all_meshes[i], most_conv_Qfiss_condensed[i]))
            abs_errors_Qfiss.append(np.abs(Qfiss_all_meshes[i] - most_conv_Qfiss_condensed[i]))
            
            rel_errors_TF.append(compute_relative_error(TF_all_meshes[i], most_conv_TF_condensed[i]))
            abs_errors_TF.append(np.abs(TF_all_meshes[i] - most_conv_TF_condensed[i]))
            
            rel_errors_TC.append(compute_relative_error(TC_all_meshes[i], most_conv_TC_condensed[i]))
            abs_errors_TC.append(np.abs(TC_all_meshes[i] - most_conv_TC_condensed[i]))
            
            rel_errors_DC.append(compute_relative_error(DC_all_meshes[i], most_conv_DC_condensed[i]))
            abs_errors_DC.append(np.abs(DC_all_meshes[i] - most_conv_DC_condensed[i]))

            rel_errors_VF.append(compute_relative_error(VF_all_meshes[i], most_conv_VF_condensed[i]))
            abs_errors_VF.append(np.abs(VF_all_meshes[i] - most_conv_VF_condensed[i]))


    RMS_Qfiss = []
    RMS_TF = []
    RMS_TC = []
    RMS_DC = []
    RMS_VF = []

    for i in range(len(rel_errors_Qfiss)):
        RMS_Qfiss.append(compute_RMS_error(Qfiss_all_meshes[i], most_conv_Qfiss_condensed[i], "percent"))
        RMS_TF.append(compute_RMS_error(TF_all_meshes[i], most_conv_TF_condensed[i], "absolute"))
        RMS_TC.append(compute_RMS_error(TC_all_meshes[i], most_conv_TC_condensed[i], "absolute"))
        RMS_DC.append(compute_RMS_error(DC_all_meshes[i], most_conv_DC_condensed[i], "absolute"))
        RMS_VF.append(compute_RMS_error(VF_all_meshes[i], most_conv_VF_condensed[i], "absolute"))

    # plot the relative errors for each field
    plot_nodal_spatial_errors(height, "power density", "%", rel_errors_Qfiss, "condensed")
    plot_nodal_spatial_errors(height, "fuel temperature", "%", rel_errors_TF, "condensed")
    plot_nodal_spatial_errors(height, "coolant temperature", "%", rel_errors_TC, "condensed")
    plot_nodal_spatial_errors(height, "coolant density", "%", rel_errors_DC, "condensed")
    plot_nodal_spatial_errors(height, "void fraction", "%", rel_errors_VF, "condensed")

    # plot absolute errors for each field

    plot_nodal_spatial_errors(height, "power density", "W/$m^3$", abs_errors_Qfiss, "condensed")
    plot_nodal_spatial_errors(height, "fuel temperature", "K", abs_errors_TF, "condensed")
    plot_nodal_spatial_errors(height, "coolant temperature", "K", abs_errors_TC, "condensed")
    plot_nodal_spatial_errors(height, "coolant density", "kg/$m^3$", abs_errors_DC, "condensed")
    plot_nodal_spatial_errors(height, "void fraction", "", abs_errors_VF, "condensed")


    # print the RMS errors for each field
    print_latex_table_from_RMS_errors_list("power density", "W/$m^3$", number_axial_slices[:-1], RMS_Qfiss)
    print_latex_table_from_RMS_errors_list("fuel temperature", "K", number_axial_slices[:-1], RMS_TF)
    print_latex_table_from_RMS_errors_list("coolant temperature", "K", number_axial_slices[:-1], RMS_TC)
    print_latex_table_from_RMS_errors_list("coolant density", "kg/$m^3$", number_axial_slices[:-1], RMS_DC)
    print_latex_table_from_RMS_errors_list("void fraction", "", number_axial_slices[:-1], RMS_VF)

def post_treat_D5vsS2(height, power_scaling_factor, number_axial_slices, pow_relax, th_relax, fuel_volume, voidFractionCorrel,  frfaccorel, P2Pcorel):
    """
    power scaling factor : list of int, power scaling factor for the DONJON5 simulations, used to study different power normalizations
    number_axial_slices : list of int, number of axial slices used to discretize the geometry
    pow_relax : float/bool, relaxation parameter for the power distribution
    th_relax : float/bool, relaxation parameter for the thermal-hydraulic parameters
    fuel_volume : float, volume of the fuel in the geometry, used to normalize the power distributions to W/m3
    """
    if height == 3.8:
        read_id = ""
    elif height == 1.555:
        read_id = "_h155"
    NODAL_errors = {}
    RMS_errors_Power = {}

    delta_keffs = {}

    donjon_powers = []
    serpent_powers = []

    D5_keffs = []
    #D5_keffs_dict = {}
    S2_keffs = []

    for p in power_scaling_factor:
        for n in number_axial_slices:
            bundle_volume = fuel_volume / n

            if p == 1:
                power = "40kW"
                ptot = 40000
            elif p == 2:
                power = "20kW"
                ptot = 20000
            elif p == 4:
                power = "10kW"
                ptot = 10000
            elif p == 8:
                power = "5kW"
                ptot = 5000
            # Recover S2 data

            path_to_S2results = f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64" #/MPHYS/mesh{n}" # path to the serpent2 results
            S2_res = f"AT10_24UOX_MPHYS{read_id}_mesh{n}_{p}_mc"
            # Recover S2 Keff
            res = st.read(f'{path_to_S2results}/{S2_res}_res.m') # read the serpent2 output file

            serpent_keff=res.resdata["absKeff"]
            S2_keffs.append(serpent_keff[0])
            #print(f"S2 Keff : {serpent_keff}")

            # Recover S2 detector data
            det = st.read(f'{path_to_S2results}/{S2_res}_det0.m') # read the serpent2 output file
            # Parse Reaction rates (fission, n_gamma, heat production) for S2 and power for DONJON5
            sum_fiss_rates, sum_n_gamma_rates, sum_heat_prod = parse_detector_rates(n, det)
            # Normalize the S2 power distribution, to compare with DONJON5, turn it into W/m3 to facilitate comparisons between meshes
            Serpent2_power_dens = normalize_S2_power(sum_heat_prod, ptot/bundle_volume)
            serpent_powers.append(Serpent2_power_dens)
            # plot the power distribution
            plot_single_distr(f"Axial Serpent2 power density, normalized to {power}", "Power density", "$W/m^3$", height, n, Serpent2_power_dens, f"S2 : {n} axial slices", f"S2_{n}_axial_slices_{p}")

            # Recover DONJON5 data
            # Parse DONJON5 Keffs
            DONJON_keffs = parse_DONJON_Keffs(height, n,  p, pow_relax, th_relax, voidFractionCorrel,  frfaccorel, P2Pcorel)
            D5_keffs.append(DONJON_keffs[-1])

            # Compare Keffs
            delta_keffs[f"Delta Keffs : {n} axial slices"] = error_Keffs_pcm(DONJON_keffs[-1], serpent_keff[0])
            print(f"Delta Keffs : {n} axial slices at power = {power} : {delta_keffs[f'Delta Keffs : {n} axial slices']} pcm")

            # multiphysics convergence of Keffs :
            plot_Keff_convergence(DONJON_keffs, height, f"D5 : {n} axial slices", f"D5_{n}_axial_slices")

            # Parse DONJON5 power distribution
            donjon5_power_dens = parse_DONJON_power(height, n, pow_relax, p, th_relax, voidFractionCorrel,  frfaccorel, P2Pcorel, ptot/bundle_volume) # normalize the power distribution to W/m3
            donjon_powers.append(donjon5_power_dens)
            # plot the power distribution
            plot_single_distr(f"Axial DONJON5 power density, normalized to {power}", "Power density", "$W/m^3$", height, n, donjon5_power_dens, f"D5 : {n} axial slices", f"D5_{n}_axial_slices_{p}")

            # Compare D5-S2 results        
            #relative_error = compute_relative_error(donjon5_power, sum_heat_prod)
            compare_results("Power", "$W/m^3$", [donjon5_power_dens], Serpent2_power_dens, height, n, f"{n} axial slices", f"{n}_axial_slices")
            nodal_error = (donjon5_power_dens - Serpent2_power_dens)*100/Serpent2_power_dens
            NODAL_errors[f"Nodal error on power : {n} axial slices, {power}"] = nodal_error
            RMS_error = compute_RMS_error(donjon5_power_dens, Serpent2_power_dens, "percent")
            RMS_errors_Power[f"RMS error on power : {n} axial slices, {power}"] = RMS_error

            # Compare axial power distributions
            #ERRORS_r[f"Power : Mesh {n}, power {p}"] = relative_error
            #quadratic_error_r[f"Power : Mesh {n}, power {p}"] = compute_quadratic_error(donjon5_power, sum_heat_prod)

        #plot_keff_D5S2_spatial_convergence(height, number_axial_slices, D5_keffs, S2_keffs, f"spatial convergence, p = {p}", f"spatial_convergence_p{p}")
    # plot keff spatial convergence
    plot_keff_error_spatial_convergence(delta_keffs, height, "spatial convergence", "spatial_convergence")

    # plot nodal errors on power for each node
    plot_nodal_errors_with_S2(height, [20, 70, 160], ["40kW", "20kW"], NODAL_errors, field="power density", unit="%")

    # plot quadratic errors on power, see how it evolves with the number of axial slices
    plot_RMS_errors("power density", "%", height, [10,20,40,50,70,80,160], RMS_errors_Power)

    
    return

def consistency_check_Qfiss(n,p,height):
    if height == 3.8:
        save_id = "h380"
    elif height == 1.555:
        save_id = "h1555"
    Qfiss_distribs = []
    # Recover MPHYS outputs for n=10, p=1
    relax_options = [(False,False), (False, 0.9), (0.5, 0.9), (0.5, False), (0.9, 0.9), (0.5, 0.8), (0.5,0.7), (0.5,0.5), (0.5,0.2), (0.2,0.2)]
    for relax in relax_options:
        Qfiss = parse_MPHYS_output_field(n, p, "Qfiss", relax[0], relax[1])
        Qfiss_distribs.append(Qfiss)
        # plot the power density distribution
        #plot_single_distr("Power density", "W/m3", n, Qfiss, f"Qfiss : {n} axial slices", f"Qfiss_{n}_axial_slices_{relax[0]}_{relax[1]}")

    compare_results("Power density", "W/m3", Qfiss_distribs, None, height, n, f"{n} axial slices", f"{n}_axial_slices")

    diffQfiss = []
    for i in range(1,len(Qfiss_distribs)):
        diffQfiss.append((Qfiss_distribs[i] - Qfiss_distribs[0])*100/Qfiss_distribs[0])

    fig,ax = plt.subplots()
    for i in range(len(diffQfiss)):
        ax.plot(diffQfiss[i], label=f"diffQfiss_{i}")
    ax.set_ylabel("Relative difference (%)")
    fig.legend()
    fig.savefig(f"Figures_{save_id}/MPHYS_diffs_Qdiff_mesh{n}_{p}_RELAXStudy.png")
    plt.close(fig)

def consistency_check_TH_params(n,p):
    TF_list = []
    TC_list = []
    DC_list = []

    # Recover MPHYS outputs for n=10, p=1
    relax_options = [(False,False), (False, 0.9), (0.9, False), (0.9, 0.9), (0.5, 0.8), (0.5,0.7), (0.5,0.5), (0.5,0.2), (0.2,0.2)]
    for relax in relax_options:
        TF = parse_MPHYS_output_field(n, p, "TeffFuel", relax[0], relax[1])
        TC = parse_MPHYS_output_field(n, p, "Twater", relax[0], relax[1])
        DC = parse_MPHYS_output_field(n, p, "rho", relax[0], relax[1])
        TF_list.append(TF)
        TC_list.append(TC)
        DC_list.append(DC)
        # plot the power density distribution
        #plot_single_distr("TFuel", "K", n, TF, f"TF : {n} axial slices", f"TF_{n}_axial_slices_{relax[0]}_{relax[1]}")
        #plot_single_distr("TCoolant", "K", n, TC, f"TC : {n} axial slices", f"TC_{n}_axial_slices_{relax[0]}_{relax[1]}")
        #plot_single_distr("DCoolant", "kg/m3", n, DC, f"DC : {n} axial slices", f"DC_{n}_axial_slices_{relax[0]}_{relax[1]}")

    compare_results("TFuel", "K", TF_list, None, n, f"{n} axial slices", f"{n}_axial_slices")
    compare_results("TCoolant", "K", TC_list, None, n, f"{n} axial slices", f"{n}_axial_slices")
    compare_results("DCoolant", "kg/m3", DC_list, None, n, f"{n} axial slices", f"{n}_axial_slices")

    diffTF = []
    diffTC = []
    diffDC = []
    for i in range(1,len(TF_list)):
        
        diffTF.append((TF_list[i] - TF_list[0])*100/TF_list[0])
        diffTC.append((TC_list[i] - TC_list[0])*100/TC_list[0])
        diffDC.append((DC_list[i] - DC_list[0])*100/DC_list[0])

    fig,ax = plt.subplots()
    for i in range(len(diffTF)):
        ax.plot(diffTF[i], label=f"diffTF_{i}")
    ax.set_ylabel("Relative difference (%)")
    fig.legend()
    fig.savefig(f"Figures_PT/MPHYS_diffs_TF_mesh{n}_{p}_RELAXStudy.png")
    plt.close(fig)

    fig,ax = plt.subplots()
    for i in range(len(diffTC)):
        ax.plot(diffTC[i], label=f"diffTC_{i+1}")
    ax.set_ylabel("Relative difference (%)")
    fig.legend()
    fig.savefig(f"Figures_PT/MPHYS_diffs_TC_mesh{n}_{p}_RELAXStudy.png")
    plt.close(fig)

    fig,ax = plt.subplots()
    for i in range(len(diffDC)):
        ax.plot(diffDC[i], label=f"diffDC_{i+1}")
    ax.set_ylabel("Relative difference (%)")
    fig.legend()
    fig.savefig(f"Figures_PT/MPHYS_diffs_DC_mesh{n}_{p}_RELAXStudy.png")
    plt.close(fig)

    return

def post_treat_GeNFoam_fields(pow, voidFractionCorrel, frfaccorel, P2Pcorel, alpha_S):

    
    # Recover GenFoam outputs
    voidFraction_GF = parse_GenFoam_output_field("voidFraction", pow, voidFractionCorrel,  frfaccorel, P2Pcorel)
    Pressure_GF = parse_GenFoam_output_field("Pressure", pow, voidFractionCorrel,  frfaccorel, P2Pcorel)
    Temperature_L_GF = parse_GenFoam_output_field("T_liquid", pow, voidFractionCorrel,  frfaccorel, P2Pcorel)
    Temperature_V_GF = parse_GenFoam_output_field("T_vapour", pow, voidFractionCorrel,  frfaccorel, P2Pcorel)
    U_L_GF = parse_GenFoam_output_field("U_liquid", pow, voidFractionCorrel,  frfaccorel, P2Pcorel)
    U_V_GF = parse_GenFoam_output_field("U_vapour", pow, voidFractionCorrel,  frfaccorel, P2Pcorel)

    # Recover MPHYS outputs
    TC = parse_MPHYS_output_field(1.555, 70, pow, "Twater", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)
    DC = parse_MPHYS_output_field(1.555, 70, pow, "rho", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)
    VF = parse_MPHYS_output_field(1.555, 70, pow, "voidFraction", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)
    P = parse_MPHYS_output_field(1.555, 70, pow, "Pressure", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)
    U = parse_MPHYS_output_field(1.555, 70, pow, "Velocity", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)
    H = parse_MPHYS_output_field(1.555, 70, pow, "Enthalpy", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)

    # exclude the first 5 points
    voidFraction_GF = voidFraction_GF[5:]
    Pressure_GF = Pressure_GF[5:]
    Temperature_L_GF = Temperature_L_GF[5:]
    Temperature_V_GF = Temperature_V_GF[5:]
    U_L_GF = U_L_GF[5:]
    U_V_GF = U_V_GF[5:]


    # renormalize GeN-Foam void fraction so that alpha_vapour + alpha_liquid = 1
    voidFraction_GF = voidFraction_GF/(1-alpha_S)
    print(f"voidFraction_GF : {voidFraction_GF}")
    print(f"VF : {VF}")

    plot_GeNFoam_vs_MPHYS("Void fraction", "", voidFraction_GF, VF, height, 70)
    plot_GeNFoam_vs_MPHYS("Pressure", "Pa", Pressure_GF, P, height, 70)
    plot_GeNFoam_vs_MPHYS_vl("Coolant temperature", "K", Temperature_L_GF, Temperature_V_GF, TC, height, 70)
    plot_GeNFoam_vs_MPHYS_vl("Coolant velocity", "m/s", U_L_GF, U_V_GF, U, height, 70)
    print("Liquid phase velocity : ", U_L_GF)
    print("Vapour phase velocity : ", U_V_GF)
    print("Mixture velocity : ", U)

    nodal_errors_P = P - Pressure_GF
    nodal_errors_VF =  VF - voidFraction_GF
    nodal_errors_TL =  TC - Temperature_L_GF
    nodal_errors_TV =  TC - Temperature_V_GF
    nodal_errors_UL = U - U_L_GF
    nodal_errors_UV =  U - U_V_GF
    
    RMS_error_P = compute_RMS_error(Pressure_GF, P, "absolute")
    RMS_error_VF = compute_RMS_error(voidFraction_GF, VF, "absolute")
    RMS_error_TL = compute_RMS_error(Temperature_L_GF, TC, "absolute")
    RMS_error_TV = compute_RMS_error(Temperature_V_GF, TC, "absolute")
    RMS_error_UL = compute_RMS_error(U_L_GF, U, "absolute")
    RMS_error_UV = compute_RMS_error(U_V_GF, U, "absolute")

    MAX_error_P = np.max(np.abs(nodal_errors_P))
    MAX_error_VF = np.max(np.abs(nodal_errors_VF))
    MAX_error_TL = np.max(np.abs(nodal_errors_TL))
    MAX_error_TV = np.max(np.abs(nodal_errors_TV))
    MAX_error_UL = np.max(np.abs(nodal_errors_UL))
    MAX_error_UV = np.max(np.abs(nodal_errors_UV))

    AVG_error_P = np.mean(np.abs(nodal_errors_P))
    AVG_error_VF = np.mean(np.abs(nodal_errors_VF))
    AVG_error_TL = np.mean(np.abs(nodal_errors_TL))
    AVG_error_TV = np.mean(np.abs(nodal_errors_TV))
    AVG_error_UL = np.mean(np.abs(nodal_errors_UL))
    AVG_error_UV = np.mean(np.abs(nodal_errors_UV))

    print(f"RMS error on Pressure : {RMS_error_P} Pa")
    print(f"RMS error on void fraction : {RMS_error_VF}")
    print(f"RMS error on Liquid temperature : {RMS_error_TL} K")
    print(f"RMS error on Vapour temperature : {RMS_error_TV} K")
    print(f"RMS error on Liquid velocity : {RMS_error_UL} m/s")
    print(f"RMS error on Vapour velocity : {RMS_error_UV} m/s")

    print(f"MAX error on Pressure : {MAX_error_P} Pa")
    print(f"MAX error on void fraction : {MAX_error_VF}")
    print(f"MAX error on Liquid temperature : {MAX_error_TL} K")
    print(f"MAX error on Vapour temperature : {MAX_error_TV} K")
    print(f"MAX error on Liquid velocity : {MAX_error_UL} m/s")
    print(f"MAX error on Vapour velocity : {MAX_error_UV} m/s")

    print(f"AVG error on Pressure : {AVG_error_P} Pa")
    print(f"AVG error on void fraction : {AVG_error_VF}")
    print(f"AVG error on Liquid temperature : {AVG_error_TL} K")
    print(f"AVG error on Vapour temperature : {AVG_error_TV} K")
    print(f"AVG error on Liquid velocity : {AVG_error_UL} m/s")
    print(f"AVG error on Vapour velocity : {AVG_error_UV} m/s")

    pressure_drop_GF = np.max(Pressure_GF) - np.min(Pressure_GF)
    pressure_drop_THM = np.max(P) - np.min(P)

    print(f"Pressure drop in GeN-Foam : {pressure_drop_GF} Pa")
    print(f"Pressure drop in THM : {pressure_drop_THM} Pa")
    error_pressure_drop = (pressure_drop_THM - pressure_drop_GF)
    print(f"Error on pressure drop : {error_pressure_drop} Pa")

    # Do the same but compute errors in %
    RMS_error_P = compute_RMS_error(Pressure_GF, P, "percent")
    RMS_error_VF = compute_RMS_error(voidFraction_GF, VF, "percent")
    RMS_error_TL = compute_RMS_error(Temperature_L_GF, TC, "percent")
    RMS_error_TV = compute_RMS_error(Temperature_V_GF, TC, "percent")
    RMS_error_UL = compute_RMS_error(U_L_GF, U, "percent")
    RMS_error_UV = compute_RMS_error(U_V_GF, U, "percent")
    
    nodal_percent_errors_P = (P - Pressure_GF)*100/Pressure_GF
    nodal_percent_errors_VF =  [(VF[i] - voidFraction_GF[i])*100/voidFraction_GF[i] if voidFraction_GF[i] != 0 else 0 for i in range(len(voidFraction_GF))]
    nodal_percent_errors_TL =  (TC - Temperature_L_GF)*100/Temperature_L_GF
    nodal_percent_errors_TV =  (TC - Temperature_V_GF)*100/Temperature_V_GF
    nodal_percent_errors_UL = (U - U_L_GF)*100/U_L_GF
    nodal_percent_errors_UV =  (U - U_V_GF)*100/U_V_GF

    MAX_percent_error_P = np.max(np.abs(nodal_percent_errors_P))
    MAX_percent_error_VF = np.max(np.abs(nodal_percent_errors_VF))
    MAX_percent_error_TL = np.max(np.abs(nodal_percent_errors_TL))
    MAX_percent_error_TV = np.max(np.abs(nodal_percent_errors_TV))
    MAX_percent_error_UL = np.max(np.abs(nodal_percent_errors_UL))
    MAX_percent_error_UV = np.max(np.abs(nodal_percent_errors_UV))

    AVG_percent_error_P = np.mean(np.abs(nodal_percent_errors_P))
    AVG_percent_error_VF = np.mean(np.abs(nodal_percent_errors_VF))
    AVG_percent_error_TL = np.mean(np.abs(nodal_percent_errors_TL))
    AVG_percent_error_TV = np.mean(np.abs(nodal_percent_errors_TV))
    AVG_percent_error_UL = np.mean(np.abs(nodal_percent_errors_UL))
    AVG_percent_error_UV = np.mean(np.abs(nodal_percent_errors_UV))

    print(f"RMS error on Pressure : {RMS_error_P} %")
    #print(f"RMS error on void fraction : {RMS_error_VF} %")
    print(f"RMS error on Liquid temperature : {RMS_error_TL} %")
    print(f"RMS error on Vapour temperature : {RMS_error_TV} %")
    print(f"RMS error on Liquid velocity : {RMS_error_UL} %")
    print(f"RMS error on Vapour velocity : {RMS_error_UV} %")

    print(f"MAX error on Pressure : {MAX_percent_error_P} %")
    #print(f"MAX error on void fraction : {MAX_percent_error_VF} %")
    print(f"MAX error on Liquid temperature : {MAX_percent_error_TL} %")
    print(f"MAX error on Vapour temperature : {MAX_percent_error_TV} %")
    print(f"MAX error on Liquid velocity : {MAX_percent_error_UL} %")
    print(f"MAX error on Vapour velocity : {MAX_percent_error_UV} %")

    print(f"AVG error on Pressure : {AVG_percent_error_P} %")
    print(f"AVG error on void fraction : {AVG_percent_error_VF} %")
    print(f"AVG error on Liquid temperature : {AVG_percent_error_TL} %")
    print(f"AVG error on Vapour temperature : {AVG_percent_error_TV} %")
    print(f"AVG error on Liquid velocity : {AVG_percent_error_UL} %")
    print(f"AVG error on Vapour velocity : {AVG_percent_error_UV} %")

    # plot the nodal errors
    plot_nodal_errors_with_GeNFoam(1.555, nodal_percent_errors_P, nodal_percent_errors_VF,
                                    nodal_percent_errors_TL, nodal_percent_errors_TV,
                                    nodal_percent_errors_UL, nodal_percent_errors_UV)


if __name__ == "__main__" :
    plot_1_2_pow_norms = True
    plot_THM_vs_GeNFoam_1_2 = True
    number_axial_slices = [10,20,40,50,70,80,160]
    power_scaling_factor = [1,2] #, 2, 4, 8] # run 4 and 8
    pitch = 1.295e-2
    height = 1.555
    clad_radius = 0.5140e-2  
    fuel_volume = np.pi*clad_radius**2*height
    alpha_S = (np.pi*clad_radius**2)/(pitch**2)
    
    pow_relax = False #pow_relax = [0.2,0.5,0.8, False]
    th_relax = False #th_relax = [0.2,0.5,0.8, False]

    voidFractionCorrel = 'EPRIvoidModel' # 'modBestion', 'HEM1', 'GEramp', 'EPRIvoidModel'
    frfaccorel = "blasius" # 'base', 'blasius', 'Churchill', 
    P2Pcorel = "lockhartMartinelli" # 'base', 'HEM1', 'HEM2', 'MNmodel', 'lockhartMartinelli'
    numericalMethod = "BiCG" # "FVM": Solves the system using matrix inversion with preconditioning.

    post_treat_D5vsS2(height, power_scaling_factor, number_axial_slices, pow_relax, th_relax, fuel_volume, voidFractionCorrel,  frfaccorel, P2Pcorel)
    #consistency_check_Qfiss(n=10, p=1)
    #consistency_check_TH_params(n=10,p=1)
    #consistency_check_TH_params(n=20,p=1)

    post_treat_Fields_SpatialConvergence_interp(height, 1, number_axial_slices, pow_relax, th_relax, voidFractionCorrel, frfaccorel, P2Pcorel)

    #number_axial_slices = [10, 20, 40, 80, 160]
    #post_treat_Fields_spatial_convergence_condense(height, 2, number_axial_slices, pow_relax, th_relax, voidFractionCorrel, frfaccorel, P2Pcorel)

    #post_treat_GeNFoam_fields(2, voidFractionCorrel, frfaccorel, P2Pcorel, alpha_S) # 

    Qfiss1 = parse_MPHYS_output_field(1.555, 160, 1, "Qfiss", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)
    Qfiss1_coarse = parse_MPHYS_output_field(1.555, 40, 1, "Qfiss", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)
    Qfiss2 = parse_MPHYS_output_field(1.555, 160, 2, "Qfiss", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)
    Qfiss2_coarse = parse_MPHYS_output_field(1.555, 40, 2, "Qfiss", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)

    VF1 = parse_MPHYS_output_field(1.555, 160, 1, "voidFraction", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)
    VF1_coarse = parse_MPHYS_output_field(1.555, 40, 1, "voidFraction", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)
    VF2 = parse_MPHYS_output_field(1.555, 160, 2, "voidFraction", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)
    VF2_coarse = parse_MPHYS_output_field(1.555, 40, 2, "voidFraction", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)


    if plot_1_2_pow_norms:
        fig, ax = plt.subplots()
        # compute z_mesh
        z_boundaries = np.linspace(0, height, 160 + 1)
        #z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
        # Generate the x and y values for plotting
        z_vals = []
        y_vals = []
        intervals = [(z_boundaries[i], z_boundaries[i+1]) for i in range(len(z_boundaries)-1)]
        for (start, end), value in zip(intervals, Qfiss1):
            z_vals.extend([start, end])  # Extend x values to cover the interval start and end
            y_vals.extend([value, value]) # Extend y values to cover the interval value
        
        ax.step(z_vals, y_vals, label=f"Power density, n=160 nodes : Ptot = 40kW")
        
        z_boundaries = np.linspace(0, height, 40 + 1)
        #z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
        # Generate the x and y values for plotting
        z_vals = []
        y_vals = []
        intervals = [(z_boundaries[i], z_boundaries[i+1]) for i in range(len(z_boundaries)-1)]
        for (start, end), value in zip(intervals, Qfiss1_coarse):
            z_vals.extend([start, end])
            y_vals.extend([value, value])
        ax.step(z_vals, y_vals, label=f"Power density n=40 nodes : Ptot = 40kW")

        z_boundaries = np.linspace(0, height, 160 + 1)
        z_vals = []
        y_vals = []
        intervals = [(z_boundaries[i], z_boundaries[i+1]) for i in range(len(z_boundaries)-1)]
        for (start, end), value in zip(intervals, Qfiss2):
            z_vals.extend([start, end])  # Extend x values to cover the interval start and end
            y_vals.extend([value, value]) # Extend y values to cover the interval value
        
        ax.step(z_vals, y_vals, label=f"Power density n=160 nodes : Ptot = 20kW")
        
        z_boundaries = np.linspace(0, height, 40 + 1)
        z_vals = []
        y_vals = []
        intervals = [(z_boundaries[i], z_boundaries[i+1]) for i in range(len(z_boundaries)-1)]
        for (start, end), value in zip(intervals, Qfiss2_coarse):
            z_vals.extend([start, end])
            y_vals.extend([value, value])
        ax.step(z_vals, y_vals, label=f"Power density n=40 nodes : Ptot = 20kW")
        
        ax.set_ylabel(f"Power density ($W/m^3$)")
        ax.set_xlabel("Height (m)")
        #ax.set_title(f"Power density distribution : Ptot = 40kW and 20kW")
        ax.grid()
        ax.legend(loc="best")
        fig.savefig(f"Figures_h1555/MPHYS_Qfiss_Ptot_1_2_piecewise_constant.png")
        plt.close(fig)


        # plotting Void fractions

        fig, ax = plt.subplots()
        # compute z_mesh
        z_boundaries = np.linspace(0, height, 160 + 1)
        #z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
        # Generate the x and y values for plotting
        z_vals = []
        y_vals = []
        intervals = [(z_boundaries[i], z_boundaries[i+1]) for i in range(len(z_boundaries)-1)]
        for (start, end), value in zip(intervals, VF1):
            z_vals.extend([start, end])  # Extend x values to cover the interval start and end
            y_vals.extend([value, value]) # Extend y values to cover the interval value
        ax.step(z_vals, y_vals, label=f"Void fraction, n=160 nodes : Ptot = 40kW")

        z_boundaries = np.linspace(0, height, 40 + 1)
        #z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
        # Generate the x and y values for plotting
        z_vals = []
        y_vals = []
        intervals = [(z_boundaries[i], z_boundaries[i+1]) for i in range(len(z_boundaries)-1)]
        for (start, end), value in zip(intervals, VF1_coarse):
            z_vals.extend([start, end])
            y_vals.extend([value, value])
        ax.step(z_vals, y_vals, label=f"Void fraction, n=40 nodes : Ptot = 40kW")

        z_boundaries = np.linspace(0, height, 160 + 1)
        z_vals = []
        y_vals = []
        intervals = [(z_boundaries[i], z_boundaries[i+1]) for i in range(len(z_boundaries)-1)]
        for (start, end), value in zip(intervals, VF2):
            z_vals.extend([start, end])  # Extend x values to cover the interval start and end
            y_vals.extend([value, value]) # Extend y values to cover the interval value
        
        ax.step(z_vals, y_vals, label=f"Void fraction, n=160 nodes : Ptot = 20kW")
        
        z_boundaries = np.linspace(0, height, 40 + 1)
        z_vals = []
        y_vals = []
        intervals = [(z_boundaries[i], z_boundaries[i+1]) for i in range(len(z_boundaries)-1)]
        for (start, end), value in zip(intervals, VF2_coarse):
            z_vals.extend([start, end])
            y_vals.extend([value, value])
        ax.step(z_vals, y_vals, label=f"Void fraction, n=40 nodes : Ptot = 20kW")

        
        ax.set_ylabel(f"Void fraction")
        ax.set_xlabel("Height (m)")
        #ax.set_title(f"Void fraction : Ptot = 40kW and 20kW")
        ax.grid()
        ax.legend(loc="best")
        fig.savefig(f"Figures_h1555/MPHYS_VF_Ptot_1_2_piecewise_constant.png")
        plt.close(fig)

    if plot_THM_vs_GeNFoam_1_2:
        # Parse GeNFoam fields for power = 40 kW and 20 kW
        voidFraction_GF_1 = parse_GenFoam_output_field("voidFraction", 1, voidFractionCorrel,  frfaccorel, P2Pcorel)
        voidFraction_GF_2 = parse_GenFoam_output_field("voidFraction", 2, voidFractionCorrel,  frfaccorel, P2Pcorel)

        Pressure_GF_1 = parse_GenFoam_output_field("Pressure", 1, voidFractionCorrel,  frfaccorel, P2Pcorel)
        Pressure_GF_2 = parse_GenFoam_output_field("Pressure", 2, voidFractionCorrel,  frfaccorel, P2Pcorel)

        Temperature_L_GF_1 = parse_GenFoam_output_field("T_liquid", 1, voidFractionCorrel,  frfaccorel, P2Pcorel)
        Temperature_L_GF_2 = parse_GenFoam_output_field("T_liquid", 2, voidFractionCorrel,  frfaccorel, P2Pcorel)

        Temperature_V_GF_1 = parse_GenFoam_output_field("T_vapour", 1, voidFractionCorrel,  frfaccorel, P2Pcorel)
        Temperature_V_GF_2 = parse_GenFoam_output_field("T_vapour", 2, voidFractionCorrel,  frfaccorel, P2Pcorel)

        # exclude the first 5 points
        voidFraction_GF_1 = voidFraction_GF_1[5:]
        voidFraction_GF_2 = voidFraction_GF_2[5:]
        Pressure_GF_1 = Pressure_GF_1[5:]
        Pressure_GF_2 = Pressure_GF_2[5:]
        Temperature_L_GF_1 = Temperature_L_GF_1[5:]
        Temperature_L_GF_2 = Temperature_L_GF_2[5:]
        Temperature_V_GF_1 = Temperature_V_GF_1[5:]
        Temperature_V_GF_2 = Temperature_V_GF_2[5:]

        # renormalize GeN-Foam void fraction so that alpha_vapour + alpha_liquid = 1
        voidFraction_GF_1 = voidFraction_GF_1/(1-alpha_S)
        voidFraction_GF_2 = voidFraction_GF_2/(1-alpha_S)
        
        # Parse THM fields for power = 40 kW and 20 kW

        TC_1 = parse_MPHYS_output_field(1.555, 70, 1, "Twater", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)
        TC_2 = parse_MPHYS_output_field(1.555, 70, 2, "Twater", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)

        #DC_1 = parse_MPHYS_output_field(1.555, 70, 1, "rho", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)
        #DC_2 = parse_MPHYS_output_field(1.555, 70, 2, "rho", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)

        VF_1 = parse_MPHYS_output_field(1.555, 70, 1, "voidFraction", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)
        VF_2 = parse_MPHYS_output_field(1.555, 70, 2, "voidFraction", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)

        P_1 = parse_MPHYS_output_field(1.555, 70, 1, "Pressure", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)
        P_2 = parse_MPHYS_output_field(1.555, 70, 2, "Pressure", False, False, voidFractionCorrel,  frfaccorel, P2Pcorel)


        # Plot pressures for powers 1 and 2, THM / GeNFoam 
        fig, ax = plt.subplots()
        # compute z_mesh
        z_boundaries = np.linspace(0, height, 70 + 1)
        #z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
        # Generate the x and y values for plotting
        z_vals = []
        y_vals = []
        intervals = [(z_boundaries[i], z_boundaries[i+1]) for i in range(len(z_boundaries)-1)]
        for (start, end), value in zip(intervals, Pressure_GF_1):
            z_vals.extend([start, end])
            y_vals.extend([value, value])
        ax.step(z_vals, y_vals, label=f"Pressure GeNFoam, Ptot = 40kW", linestyle="--", color = "red") 
        z_vals = []
        y_vals = []
        for (start, end), value in zip(intervals, P_1):
            z_vals.extend([start, end])
            y_vals.extend([value, value])
        ax.step(z_vals, y_vals, label=f"Pressure THM, Ptot = 40kW", linestyle="--", color = "green")
        z_vals = []
        y_vals = []
        for (start, end), value in zip(intervals, Pressure_GF_2):
            z_vals.extend([start, end])
            y_vals.extend([value, value])
        ax.step(z_vals, y_vals, label=f"Pressure GeNFoam, Ptot = 20kW", color = "red")
        z_vals = []
        y_vals = []
        for (start, end), value in zip(intervals, P_2):
            z_vals.extend([start, end])
            y_vals.extend([value, value])
        ax.step(z_vals, y_vals, label=f"Pressure THM, Ptot = 20kW", color = "green")
        ax.set_ylabel(f"Pressure (Pa)")
        ax.set_xlabel("Height (m)")
        #ax.set_title(f"Pressure : THM vs GeNFoam comparison")
        ax.grid()
        ax.legend(loc="best")
        fig.savefig(f"Figures_GeNFoam_THM/THM_vs_GeNFoam_Pressure_pow_1_2.png")
        plt.close(fig)

        # Plot void fractions for powers 1 and 2, THM / GeNFoam
        fig, ax = plt.subplots()
        # compute z_mesh
        z_boundaries = np.linspace(0, height, 70 + 1)
        #z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
        # Generate the x and y values for plotting
        z_vals = []
        y_vals = []
        intervals = [(z_boundaries[i], z_boundaries[i+1]) for i in range(len(z_boundaries)-1)]
        for (start, end), value in zip(intervals, voidFraction_GF_1):
            z_vals.extend([start, end])
            y_vals.extend([value, value])
        ax.step(z_vals, y_vals, label=f"Void fraction GeNFoam, Ptot = 40kW", linestyle="--", color="red")
        z_vals = []
        y_vals = []
        for (start, end), value in zip(intervals, VF_1):
            z_vals.extend([start, end])
            y_vals.extend([value, value])
        ax.step(z_vals, y_vals, label=f"Void fraction THM, Ptot = 40kW", linestyle="--", color="green")
        z_vals = []
        y_vals = []
        for (start, end), value in zip(intervals, voidFraction_GF_2):
            z_vals.extend([start, end])
            y_vals.extend([value, value])
        ax.step(z_vals, y_vals, label=f"Void fraction GeNFoam, Ptot = 20kW", color = "red")
        z_vals = []
        y_vals = []
        for (start, end), value in zip(intervals, VF_2):
            z_vals.extend([start, end])
            y_vals.extend([value, value])
        ax.step(z_vals, y_vals, label=f"Void fraction THM, Ptot = 20kW", color = "green")
        ax.set_ylabel(f"Void fraction")
        ax.set_xlabel("Height (m)")
        #ax.set_title(f"Void fraction : THM vs GeNFoam comparison")
        ax.grid()
        ax.legend(loc="best")
        fig.savefig(f"Figures_GeNFoam_THM/THM_vs_GeNFoam_Void_fraction_pow_1_2.png")
        plt.close(fig)


        # Plot liquid temperatures, vapour temperatures, mixture temperatures for powers 1 and 2, THM / GeNFoam
        fig, ax = plt.subplots()
        # compute z_mesh
        z_boundaries = np.linspace(0, height, 70 + 1)
        #z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
        # Generate the x and y values for plotting
        z_vals = []
        y_vals = []
        intervals = [(z_boundaries[i], z_boundaries[i+1]) for i in range(len(z_boundaries)-1)]
        for (start, end), value in zip(intervals, Temperature_L_GF_1):
            z_vals.extend([start, end])
            y_vals.extend([value, value])
        ax.step(z_vals, y_vals, label=f"T liquid GeNFoam, Ptot = 40kW", linestyle="--", color = "blue")
        z_vals = []
        y_vals = []
        for (start, end), value in zip(intervals, Temperature_L_GF_2):
            z_vals.extend([start, end])
            y_vals.extend([value, value])
        ax.step(z_vals, y_vals, label=f"T liquid GeNFoam, Ptot = 20kW", color = "blue", linestyle="-")
        z_vals = []
        y_vals = []
        for (start, end), value in zip(intervals, Temperature_V_GF_1):
            z_vals.extend([start, end])
            y_vals.extend([value, value])
        ax.step(z_vals, y_vals, label=f"T vapor GeNFoam, Ptot = 40kW", linestyle="--", color = "red")
        z_vals = []
        y_vals = []
        for (start, end), value in zip(intervals, Temperature_V_GF_2):
            z_vals.extend([start, end])
            y_vals.extend([value, value])
        ax.step(z_vals, y_vals, label=f"T vapor GeNFoam, Ptot = 20kW", linestyle="-", color = "red")
        z_vals = []
        y_vals = []
        for (start, end), value in zip(intervals, TC_1):
            z_vals.extend([start, end])
            y_vals.extend([value, value])
        ax.step(z_vals, y_vals, label=f"T mixture THM, Ptot = 40kW", linestyle="--", color = "green")
        z_vals = []
        y_vals = []
        for (start, end), value in zip(intervals, TC_2):
            z_vals.extend([start, end])
            y_vals.extend([value, value])
        ax.step(z_vals, y_vals, label=f"T mixture THM, Ptot = 20kW", linestyle="-", color = "green")
        ax.set_ylabel(f"Coolant temperature (K)")
        ax.set_xlabel("Height (m)")
        #ax.set_title(f"Coolant temperature : THM vs GeNFoam comparison")
        ax.grid()
        ax.legend(loc="best")
        fig.savefig(f"Figures_GeNFoam_THM/THM_vs_GeNFoam_Coolant_temperature_pow_1_2.png")
        plt.close(fig)








    




