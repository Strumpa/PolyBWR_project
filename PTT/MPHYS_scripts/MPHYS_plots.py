### Seperate File regrouping all the plotting functions for the MPHYS project
## Author: R. Guasch, re-organsied on 21/03/2025
# List of functions:
# - plot_nodal_errors_with_GeNFoam(height, nodal_errors_P, nodal_errors_VF, nodal_errors_TL, nodal_errors_TV, nodal_errors_UL, nodal_errors_UV)
# - plot_GeNFoam_vs_MPHYS_vl(field, unit, GF_L_field, GF_V_field, MPHYS_field, height, n_axial_nodes)
# - plot_GeNFoam_vs_MPHYS(field, unit, GF_field, MPHYS_field, height, n_axial_nodes)
# - plot_piecewise_constant_field(field, unit, field_list, height, n_slices_list)
# - plot_spatial_convergence_field(field, unit, fields_list, height, n_slices_list)
# - plot_RMS_errors(field, unit, height, n_slices_list, quadratic_error_dict)
# - plot_nodal_spatial_errors(height, field, unit, nodal_spatial_errors, interp_type)
# - plot_nodal_errors_with_S2(height, n_slices_list, power_list, nodal_errors_dict, field, unit)
# - plot_keff_D5S2_spatial_convergence(height, n_slices, D5_keffs, S2_keffs, name, save_name)
# - plot_keff_error_spatial_convergence(delta_keffs_dict, height, name, save_name)
# - compare_results(field, unit, donjon5_fields, S2_field, height, n_slices, name, save_name)
# - plot_relative_error(relative_error_dict, height, type, power_scaling_factors_to_plot)
# - plot_condensed_relative_error(relative_error_dict, height, type, power_scaling_factors_to_plot)
# - plot_quadratic_errors(quadratic_error_dict, type, max_slices)
# - plot_single_distr(title, field, unit, height,  n_slices, power_distr, name, save_name)
# - plot_Keff_convergence(keffs, height, name, save_name)

import numpy as np
import matplotlib.pyplot as plt

### Plotting functions
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
    
    return

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
    
    return

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



def compare_all_numbers_slices(height, field, unit, donjon_field_list, serpent_field_list, n_slices_list):
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





