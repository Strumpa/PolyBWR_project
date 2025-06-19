# Python3 script to retrive the absorption rates computed by the VDG: module
# Author : R. Guasch
# Date : 2024/08/06
# Purpose : analysis of absorption rates for U238, Gd155 and Gd157 for homogeneous cell study.
# Usage : python3 VDG_abs_rates.py


import numpy as np
import matplotlib.pyplot as plt
import re
import os

def parse_error_abs_rates(filename, input_case, output_case, grmin=52, grmax=206):
    """
    This function parses the HOM_C7_autop.result file and retrieves the relative errors on the absorption rates for U238, Gd155 and Gd157.
    The relative error is computed taking the AUTOSECOL (AUTO:) value as the reference. In HOM_C7_autop.x2m, it comapres the it to RSE, PT, and SUBG values.
    """
    # Open the file
    if input_case == "HOM_U5":
        isotopes = ["U235"]
        isotope_counter = {"U235": 0, "ALL": 0}
    elif input_case == "HOM_U5_U8":
        isotopes = ["U235", "U238"]
        isotopes_to_parse = isotopes
        isotope_counter = {"U235": 0, "U238": 0, "ALL": 0}
    elif input_case == "HOM_UOX":
        isotopes = ["U234", "U235", "U238"]
        isotopes_to_parse = isotopes
        isotope_counter = {"U234": 0, "U235": 0, "U238": 0, "ALL": 0}
    elif input_case == "HOM_UOX_Gd155":
        isotopes = ["U234", "U235", "U238", "Gd155"]
        isotopes_to_parse = isotopes
        isotope_counter = {"U234": 0, "U235": 0, "U238": 0, "Gd155": 0, "ALL": 0}
    elif input_case == "HOM_UOX_Gd157":
        isotopes = ["U234", "U235", "U238", "Gd157"]
        isotopes_to_parse = isotopes
        isotope_counter = {"U234": 0, "U235": 0, "U238": 0, "Gd157": 0, "ALL": 0}
    elif input_case == "HOM_UOXGd":
        isotopes = ["U234", "U235", "U238", "Gd155", "Gd157"]
        isotopes_to_parse = isotopes
        isotope_counter = {"U234": 0, "U235": 0, "U238": 0, "Gd155": 0, "Gd157": 0, "ALL": 0}
    length_VDG_table = grmax - grmin + 1
    with open(filename) as f:
        lines = f.readlines()

    print(f"length of isotopes = {len(isotopes)}")
    #print(isotopes_to_parse)
    if len(isotopes) == 1:
        for i in range(len(lines)):
            if f"SUM OF ALL ISOTOPES" in lines[i]:
                print(f"{isotopes[0]} line index is {i}")
                parsed_data = lines[i+3:i+length_VDG_table+3]
                groups, errors = get_groups_and_errors(parsed_data)
                plot_relative_error_histogram(groups, grmin, grmax,errors, isotopes[0], isotope_counter[isotopes[0]], input_case, output_case)
                parsed_MAX_INT_AVG_error = lines[i+length_VDG_table+5:i+length_VDG_table+9]
                print(f"Parsed MAX_INT_AVG_error for {isotopes[0]}")
                print(parsed_MAX_INT_AVG_error)
                generate_MAX_INT_AVG_error_table(parsed_MAX_INT_AVG_error, isotopes[0], isotope_counter[isotopes[0]])
                isotope_counter[isotopes[0]] += 1 
    else: 
        for iso in isotopes_to_parse:
            print(f"iso = {iso}")
            for i in range(len(lines)):
                if f" ISOTOPE ='{iso}" in lines[i]:
                    parsed_data = lines[i+3:i+length_VDG_table+3]
                    groups, errors = get_groups_and_errors(parsed_data)
                    print(f"iso = {iso}, counter = {isotope_counter[iso]}")
                    plot_relative_error_histogram(groups, grmin, grmax, errors, iso, isotope_counter[iso], input_case, output_case)
                    parsed_MAX_INT_AVG_error = lines[i+length_VDG_table+5:i+length_VDG_table+9]
                    print(f"Parsed MAX_INT_AVG_error for {iso}")
                    print(parsed_MAX_INT_AVG_error)
                    generate_MAX_INT_AVG_error_table(parsed_MAX_INT_AVG_error, iso, isotope_counter[iso])
                    isotope_counter[iso] += 1 
        for i in range(len(lines)):
            if "SUM OF ALL ISOTOPES" in lines[i]:
                print("in ALL ISOTOPES")
                parsed_data = lines[i+3:i+length_VDG_table+3]
                groups, errors = get_groups_and_errors(parsed_data)
                plot_relative_error_histogram(groups, grmin, grmax, errors, "ALL", isotope_counter["ALL"], input_case, output_case)
                parsed_MAX_INT_AVG_error = lines[i+length_VDG_table+5:i+length_VDG_table+9]
                generate_MAX_INT_AVG_error_table(parsed_MAX_INT_AVG_error, "ALL", isotope_counter["ALL"])
                isotope_counter["ALL"] += 1
    
    for i in range(len(lines)):
        if ">|Kinf AUTO =" in lines[i]:
            Kinf_AUTO = lines[i].strip().split("=")[1].strip().split(" ")[0]
        if ">|Kinf RSE =" in lines[i]:
            Kinf_RSE = lines[i].strip().split("=")[1].strip().split(" ")[0]
        if ">|Kinf PT =" in lines[i]:
            Kinf_PT = lines[i].strip().split("=")[1].strip().split(" ")[0]
        if ">|Kinf SUBG =" in lines[i]:
            Kinf_SUBG = lines[i].strip().split("=")[1].strip().split(" ")[0]

    compute_error_to_S2(input_case, Kinf_AUTO, Kinf_RSE, Kinf_PT, Kinf_SUBG)
    return

def get_groups_and_errors(lines):
    groups = []
    errors = []
    for line in lines:
        line = line.strip()
        if line:
            group, error = line.split(" ")[0], line.split(" ")[-2] 
            groups.append(int(group))
            errors.append(float(error))

    return groups, errors

def plot_relative_error_histogram(groups, grmin, grmax, relative_errors, isotope, counter, input_case, output_case):
    """
    Plots a histogram of relative errors in percentage with energy groups labeled from 1 to 295.

    Parameters:
    relative_errors (list or array-like): A list or array containing the relative errors in percentage.

    Returns:
    None
    """
     # Number of bins and energy groups
    num_bins = grmax - grmin + 1
    energy_groups = list(range(grmin, grmax + 1))
    tests = ["RSE", "PT", "SUBG"]

    if isotope == "ALL":
        isotope = "homogeneous mixture"
    
    # Ensure the length of relative_errors matches the number of energy groups
    if len(relative_errors) != num_bins:
        raise ValueError(f"The length of relative_errors must be {num_bins}.")

    # Plotting the data
    plt.figure(figsize=(12, 6))
    plt.bar(energy_groups, relative_errors, color='skyblue', edgecolor='black')

    # Adding titles and labels
    plt.title(f'Relative Error on absorption rate of {isotope} for {tests[counter]} vs AUTOSECOL')
    plt.xlabel('Energy Group')
    plt.ylabel('Relative Error (%)')
    plt.xticks(ticks=range(energy_groups[0], energy_groups[-1] + 1, 20))  # Set x-axis ticks every 20 groups
    path = os.getcwd()
    os.chdir(f"{path}/VDG_errors/{input_case}")
    plt.savefig(f"relative_errors_{isotope}_{tests[counter]}_{output_case}.png")
    # Display the plot
    plt.show()
    plt.close()
    os.chdir(path)
    
def generate_MAX_INT_AVG_error_table(parsed_MAX_INT_AVG_error, isotope, counter):
    """
    Generate a Latex table which presents the maximum error and its associated energy group, the average error, and the integrated error.
    """
    tests = ["RSE", "PT", "SUBG"]
    method = tests[counter]
    MAX = float(parsed_MAX_INT_AVG_error[1].strip().split("=")[-1].split("%")[0])
    print(f"MAX={MAX}")
    GROUP = int(parsed_MAX_INT_AVG_error[1].strip().split("IN GROUP")[-1])
    print(f"GROUP={GROUP}")
    AVG = float(parsed_MAX_INT_AVG_error[2].strip().split("=    ")[-1].split("%")[0])
    print(f"AVG={AVG}")
    INT = float(parsed_MAX_INT_AVG_error[3].strip().split("=")[-1].split("%")[0])
    # Begin the LaTeX table
    table = r"""
            \begin{table}[h!]
            \centering
            \begin{tabular}{|l|l|c|c|c|}
            \hline
            \textbf{Isotope} & \textbf{Method} & \textbf{Maximum Error (\%)} & \textbf{Average Error (\%)} & \textbf{Integrated Error (\%)} \\
            """
    table +=  f"{isotope} & {method} & {MAX:.3f} (Group {GROUP}) & {AVG:.3f} & {INT:.3f} \\\\"

        
        # End the LaTeX table
    table += r"""\hline
            \end{tabular}
            \caption{Maximum, average and integrated error for isotope : fill_this}
            \label{table:fill_this}
            \end{table}
            """
    print(table)

    return

def compute_error_to_S2(input_case, Kinf_AUTO, Kinf_RSE, Kinf_PT, Kinf_SUBG):
    """
    This is kinda stupid and hard coded for now. It computes the error to S2 for the Kinf values.   
    """
    print(f"input_case = {input_case}")
    if input_case == "HOM_U5":
        Kinf_S2_Pynjoy2016 = 1.88762
        Kinf_S2_sss_jeff311 = 1.88743
    elif input_case == "HOM_U5_U8":
        Kinf_S2_Pynjoy2016 = 1.28355
        Kinf_S2_sss_jeff311 = 1.28260
    elif input_case == "HOM_UOX":
        Kinf_S2_Pynjoy2016 = 1.27877
        Kinf_S2_sss_jeff311 = 1.27794
    elif input_case == "HOM_UOX_Gd155":
        Kinf_S2_Pynjoy2016 = 5.43432E-01
        Kinf_S2_sss_jeff311 = 5.42958E-01
    elif input_case == "HOM_UOX_Gd157":
        Kinf_S2_Pynjoy2016 = 4.65818E-01
        Kinf_S2_sss_jeff311 = 4.65498E-01
    elif input_case == "HOM_UOXGd":
        Kinf_S2_Pynjoy2016 = 4.20426E-01
        Kinf_S2_sss_jeff311 = 4.20228E-01

    error_AUTO_PyNjoy2016 = (float(Kinf_AUTO) - Kinf_S2_Pynjoy2016)*1e5
    error_RSE_PyNjoy2016 = (float(Kinf_RSE) - Kinf_S2_Pynjoy2016)*1e5
    error_PT_PyNjoy2016 = (float(Kinf_PT) - Kinf_S2_Pynjoy2016)*1e5
    error_SUBG_PyNjoy2016 = (float(Kinf_SUBG) - Kinf_S2_Pynjoy2016)*1e5

    error_AUTO_sss_jeff311 = (float(Kinf_AUTO) - Kinf_S2_sss_jeff311)*1e5
    error_RSE_sss_jeff311 = (float(Kinf_RSE) - Kinf_S2_sss_jeff311)*1e5
    error_PT_sss_jeff311 = (float(Kinf_PT) - Kinf_S2_sss_jeff311)*1e5
    error_SUBG_sss_jeff311 = (float(Kinf_SUBG) - Kinf_S2_sss_jeff311)*1e5

    print(f"Error to S2 for AUTO vs PyNjoy2016 is {error_AUTO_PyNjoy2016:.0f} pcm")
    print(f"Error to S2 for RSE vs PyNjoy2016 is {error_RSE_PyNjoy2016:.0f} pcm")
    print(f"Error to S2 for PT vs PyNjoy2016 is {error_PT_PyNjoy2016:.0f} pcm")
    print(f"Error to S2 for SUBG vs PyNjoy2016 is {error_SUBG_PyNjoy2016:.0f} pcm")

    print(f"Error to S2 for AUTO vs sss_jeff311 is {error_AUTO_sss_jeff311:.0f} pcm")
    print(f"Error to S2 for RSE vs sss_jeff311 is {error_RSE_sss_jeff311:.0f} pcm")
    print(f"Error to S2 for PT vs sss_jeff311 is {error_PT_sss_jeff311:.0f} pcm")
    print(f"Error to S2 for SUBG vs sss_jeff311 is {error_SUBG_sss_jeff311:.0f} pcm")
    data = {
    "AUTO": {"PyNjoy2016": error_AUTO_PyNjoy2016, "sss_jeff311": error_AUTO_sss_jeff311},
    "RSE": {"PyNjoy2016": error_RSE_PyNjoy2016, "sss_jeff311": error_RSE_sss_jeff311},
    "PT": {"PyNjoy2016": error_PT_PyNjoy2016, "sss_jeff311": error_PT_sss_jeff311},
    "SUBG": {"PyNjoy2016": error_SUBG_PyNjoy2016, "sss_jeff311": error_SUBG_sss_jeff311}
    }
    
    # Begin the LaTeX table
    table = r"""
            \begin{table}[h!]
            \centering
            \begin{tabular}{|l|c|c|}
            \hline
            \textbf{Self-Shielding method} & $\Delta k_{eff}$  \textbf{vs S2 PyNjoy2016 (pcm)} & $\Delta k_{eff}$  \textbf{vs S2 oldlib (pcm)} \\
            \hline
            """
                
    # Add the data rows
    for method, errors in data.items():
            table += f"{method} & {errors['PyNjoy2016']:.0f} & {errors['sss_jeff311']:.0f} \\\\ \n"
        
        # End the LaTeX table
    table += r"""\hline
            \end{tabular}
            \caption{D5-S2 bias on eigenvalues depending on the self-shielding method, fill_this}
            \label{table:fill_this}
            \end{table}
            """
    print(table)



# Parse the {input_case}_autop.result file
input_case = "HOM_UOX_Gd157"
correlation = "_CORR" # "_CORR", "_noCORR"
path=os.getcwd()
path_exists = os.path.exists(f"{path}/VDG_errors/{input_case}")
if not path_exists:
    os.makedirs(f"{path}/VDG_errors/{input_case}")
parse_error_abs_rates(f"/home/p117902/working_dir/PolyBWR_project/Version5_wc/Dragon/Linux_aarch64/{input_case}_autop{correlation}.result", input_case, output_case = f"{input_case}_USS_AUTO_inrs1{correlation}")

"""
for case in ["HOM_U5_U8", "HOM_UOX", "HOM_UOX_Gd155", "HOM_UOX_Gd157", "HOM_UOXGd"]:
    for correlation in ["_CORR", "_noCORR"]:
        print(f"$$$- VDG: error analysis for {case} with {correlation} -$$$")
        input_case = case
        path_exists = os.path.exists(f"{path}/VDG_errors/{input_case}")
        if not path_exists:
            os.makedirs(f"{path}/VDG_errors/{input_case}")
        parse_error_abs_rates(f"/home/p117902/working_dir/PolyBWR_project/Version5_wc/Dragon/Linux_aarch64/{input_case}_autop{correlation}.result", input_case, output_case = f"{input_case}_USS_AUTO_inrs1{correlation}")
"""
