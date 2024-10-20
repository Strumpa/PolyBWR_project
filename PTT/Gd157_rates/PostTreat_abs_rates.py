# Python3 script to post treat Serpent output files
# Author : R. Guasch
# Date : 15/10/2024
# Purpose : Post treat Serpent output mdx and det files to compare vs Dragon reaction rates.

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import serpentTools as st

def plot_histogram(values, iso, grmin, grmax):
    print(f"there are {len(values)} values")
    # Check that the number of values is 295
    if len(values) != (grmax - grmin) + 1:
        raise ValueError(f"All input rates lists must contain exactly {grmax - grmin + 1} values.")
    
    # Create bins from 1 to 295
    energy_groups = list(range(grmin, grmax + 1))
    #colors = ['skyblue', 'red']
    # Plot the histogram
    plt.figure(figsize=(10, 6))
    plt.bar(energy_groups, values, edgecolor='black', alpha=0.5, color="skyblue")
    plt.xticks(range(1, 296, 20))  # Set x-ticks with an interval of 20 for readability
    plt.xlabel('Energy group')
    plt.ylabel('Reaction rate')
    plt.title('Histogram of 295 Values')
    #plt.grid(True)
    #plt.tight_layout()
    plt.savefig(f'histogram_{iso}.png')
    plt.show()

def plot_histograms(values_lists, iso, colors, grmin, grmax, legends=None, save_dir=None, save_name=None):
    # Check that the number of values is 295
    for values in values_lists:
        print(len(values))
        if len(values) != grmax - grmin + 1:
            raise ValueError("All input rates lists must contain exactly 295 values.")
    
    # Create bins from 1 to 295
    energy_groups = list(range(grmin, grmax + 1))
    #colors = ['skyblue', 'red']
    # Plot the histogram
    
    plt.figure(figsize=(10, 6))
    for i in range(len(values_lists)):
        print(f"values_lists[{i}] = {values_lists[i]}")
        plt.bar(energy_groups, values_lists[i], edgecolor='black', alpha=0.5, color=colors[i], label=legends[i])
    if "error" in legends[0]:
        plt.ylabel('Relative error (%)')
    elif "ratio" in save_name:
        plt.ylabel('Ratio of absorption rates')
    else:
        plt.ylabel('Absorption rate')
    plt.xticks(range(grmin, grmax+1, 20))  # Set x-ticks with an interval of 20 for readability
    plt.xlabel('Energy group')
    if "error" in legends[0]:
        plt.title(f'Relative errors on absorption rates for {iso}')
    else:
        plt.title(f'Absorption rates for {iso}')
    plt.grid()
    plt.legend()
    #plt.tight_layout()
    if save_dir:
        plt.savefig(f'{save_dir}/{save_name}.png')
    plt.show()
    plt.close()

def parse_Serpent_detector(path_to_serpent_results, library, bu_step):
    # Read detector file
    det = st.read(f"{path_to_serpent_results}/HOM_UOX_Gd157_XS_{library}_mc_det{bu_step}.m")
    # Get detector names
    Gd_det = det.detectors["Gd_det"].tallies
    #print(Gd_det.indexes)
    #print(det_scores.bins.shape)
    n_groups = Gd_det.shape[0]
    n_reactions = Gd_det.shape[1]
    return Gd_det, n_groups, n_reactions

def parse_Serpent_microdepletion(path_to_serpent_results, library, bu_step):
    # Read mdet file
    mdet = st.read(f"{path_to_serpent_results}/HOM_UOX_Gd157_XS_{library}_mc_mdx{bu_step}.m")
    # Get mdet names
    mdet_scores = mdet.detectors.tallies
    print(mdet_scores.shape)
    return

def parse_Dragon_abs_rates(filename, input_case, output_case, grmin=1, grmax=295):
    """
    This function parses the HOM_*_autop.result file and retrieves the relative errors on the absorption rates for U238, Gd155 and Gd157.
    The relative error is computed taking the AUTOSECOL (AUTO:) value as the reference. In HOM_*_autop.x2m, AUTO values are compared to RSE, PT, and SUBG values.
    """
    RATES_ = {}
    ERROR_RATES = {}
    USS_methods = ["RSE", "PT", "SUBG"]
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
                groups, errors, AUTO_rates, USS_rates  = get_groups_rates_and_VDG_errors(parsed_data)
                RATES_[f"AUTO {iso}, {USS_methods[isotope_counter[0]]}"] = AUTO_rates
                RATES_[f"USS {iso}, {USS_methods[isotope_counter[0]]}"] = USS_rates
                ERROR_RATES[f"{iso}, {USS_methods[isotope_counter[0]]}-AUTO"] = errors
                plot_VDG_error_histogram(groups, grmin, grmax, errors, isotopes[0], isotope_counter[isotopes[0]], input_case, output_case)
                parsed_MAX_INT_AVG_error = lines[i+length_VDG_table+5:i+length_VDG_table+9]
                print(f"Parsed MAX_INT_AVG_error for {isotopes[0]}")
                print(parsed_MAX_INT_AVG_error)
                #generate_MAX_INT_AVG_error_table(parsed_MAX_INT_AVG_error, isotopes[0], isotope_counter[isotopes[0]])
                isotope_counter[isotopes[0]] += 1 
    else: 
        for iso in isotopes_to_parse:
            print(f"iso = {iso}")
            for i in range(len(lines)):
                if f" ISOTOPE ='{iso}" in lines[i]:
                    parsed_data = lines[i+3:i+length_VDG_table+3]
                    print(f"iso = {iso}")
                    print(f"ssh_method = {USS_methods[isotope_counter[iso]]}")
                    groups, errors, AUTO_rates, USS_rates  = get_groups_rates_and_VDG_errors(parsed_data)
                    RATES_[f"AUTO {iso}, {USS_methods[isotope_counter[iso]]}"] = AUTO_rates
                    RATES_[f"USS {iso}, {USS_methods[isotope_counter[iso]]}"] = USS_rates
                    ERROR_RATES[f"{iso}, {USS_methods[isotope_counter[iso]]}-AUTO"] = errors
                    print(f"iso = {iso}, counter = {isotope_counter[iso]}")
                    plot_VDG_error_histogram(groups, grmin, grmax, errors, iso, isotope_counter[iso], input_case, output_case)
                    parsed_MAX_INT_AVG_error = lines[i+length_VDG_table+5:i+length_VDG_table+9]
                    print(f"Parsed MAX_INT_AVG_error for {iso}")
                    print(parsed_MAX_INT_AVG_error)
                    #generate_MAX_INT_AVG_error_table(parsed_MAX_INT_AVG_error, iso, isotope_counter[iso])
                    isotope_counter[iso] += 1 
        for i in range(len(lines)):
            if "SUM OF ALL ISOTOPES" in lines[i]:
                print("in ALL ISOTOPES")
                parsed_data = lines[i+3:i+length_VDG_table+3]
                groups, errors, AUTO_rates, USS_rates  = get_groups_rates_and_VDG_errors(parsed_data)
                print(f"iso = ALL, counter = {isotope_counter['ALL']}")
                print(f"Relative errors = {errors}")
                #RATES_[f"AUTO ALL, {USS_methods[isotope_counter["ALL"]]}"] = AUTO_rates
                #RATES_[f"USS ALL, {USS_methods[isotope_counter["ALL"]]}"] = USS_rates
                #ERROR_RATES[f"ALL, {USS_methods[isotope_counter["ALL"]]}-AUTO"] = errors
                plot_VDG_error_histogram(groups, grmin, grmax, errors, "ALL", isotope_counter["ALL"], input_case, output_case)
                parsed_MAX_INT_AVG_error = lines[i+length_VDG_table+5:i+length_VDG_table+9]
                #generate_MAX_INT_AVG_error_table(parsed_MAX_INT_AVG_error, "ALL", isotope_counter["ALL"])
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

    RATES_["AUTO Gd157"] = np.array(RATES_["AUTO Gd157, RSE"])
    RATES_["AUTO U238"] = np.array(RATES_["AUTO U238, RSE"])
    RATES_["AUTO U235"] = np.array(RATES_["AUTO U235, RSE"])
    RATES_["AUTO U234"] = np.array(RATES_["AUTO U234, RSE"])
    RATES_.pop("AUTO U234, RSE")
    RATES_.pop("AUTO U235, RSE")
    RATES_.pop("AUTO U238, RSE")
    RATES_.pop("AUTO Gd157, RSE")
    RATES_.pop("AUTO U234, PT")
    RATES_.pop("AUTO U235, PT")
    RATES_.pop("AUTO U238, PT")
    RATES_.pop("AUTO Gd157, PT")
    RATES_.pop("AUTO U234, SUBG")
    RATES_.pop("AUTO U235, SUBG")
    RATES_.pop("AUTO U238, SUBG")
    RATES_.pop("AUTO Gd157, SUBG")
    print(f"keys = {RATES_.keys()}")
    return RATES_, ERROR_RATES

def get_groups_rates_and_VDG_errors(lines):
    groups = []
    errors = []
    abs_rates_AUTO = []
    abs_rates_USS = []
    for line in lines:
        line = line.strip()
        line = line.split(" ")
        line = [word for word in line if word != ""]
        if line:
            print(line)
            if len(line) == 5:
                group, AUTO_rate, USS_rate, error = line[0], line[1], line[2], line[-2] 
                print(f"error in % is {error}")
            elif len(line) == 4:
                print(f"line is = {line}")
                group,  error = line[0], line[-2]
                print(line[1].split("-"))
                print(f"error in % is {error}")
                AUTO_rate, USS_rate = "-"+line[1].split("-")[1]+line[1].split("-")[2], "-"+line[1].split("-")[3]+line[1].split("-")[4]
            groups.append(int(group))
            errors.append(float(error))
            abs_rates_AUTO.append(float(AUTO_rate))
            abs_rates_USS.append(float(USS_rate))

    return groups, errors, abs_rates_AUTO, abs_rates_USS

def plot_VDG_error_histogram(groups, grmin, grmax, relative_errors, isotope, counter, input_case, output_case):
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

def get_S2_rates_for_iso(S2_rates, iso):
    """
    Get the absorption rates for a specific isotope from the Serpent2 results.
    """
    iso_S2 = np.array(S2_rates[f"{iso} (n,gamma)"][::-1])
    return iso_S2

def get_D5_rates_for_iso(D5_rates, iso, SSH_method):
    """
    Get the absorption rates for a specific isotope from the DRAGON results.
    """
    if SSH_method == "AUTO":
        iso_D5 = np.array(D5_rates[f"AUTO {iso}"])
    else:
        iso_D5 = np.array(D5_rates[f"USS {iso}, {SSH_method}"])
    return iso_D5

def compare_iso_U8_absorption_rates(D5_rates, S2_rates, SSH_method, iso, N_U238, N_iso, grmin=1, grmax=295, input_case="", correlation="", path_save="", local_path=""):
    """
    D5_rates (dict): Dictionary containing the absorption rates from DRAGON/AUTO. key is test's name ("{Module} {isotope} ({method} if Module = USS)")
    S2_rates (dict): Dictionary containing the absorption rates from Serpent. key is test's name
    """
    if SSH_method == "AUTO":
        U8_D5 = np.array(D5_rates["AUTO U238"])
        iso_D5 = np.array(D5_rates[f"AUTO {iso}"])
    else:
        U8_D5 = np.array(D5_rates[f"USS U238, {SSH_method}"])
        iso_D5 = np.array(D5_rates[f"USS {iso}, {SSH_method}"])

    U8_S2 = np.array(S2_rates["U238 (n,gamma)"])*N_U238
    iso_S2 = np.array(S2_rates[f"{iso} (n,gamma)"])*N_iso
    plot_histogram(iso_S2, iso+"_S2", grmin, grmax)
    plot_histogram(iso_D5, iso+"_D5", grmin, grmax)

    print(f"length of U8_D5 = {len(U8_D5)}")

    ratio_D5 = iso_D5/U8_D5
    ratio_S2 = iso_S2/U8_S2

    fig, ax = plt.subplots()
    ax.plot(ratio_D5, label=f"D5 {iso}/U8 absorption rates : {SSH_method}")
    ax.plot(ratio_S2, label=f"Serpent2 {iso}/U8 absorption rates")
    ax.set_xlabel("Energy group")
    ax.set_ylabel("Ratio")
    ax.set_title(f"{iso}/U8 absorption rates")
    ax.legend()
    os.chdir(path_save)
    fig.savefig(f"{iso}_U8_ratio_absorption_rates_{input_case}{correlation}_{SSH_method}.png")
    os.chdir(local_path)

    plot_histograms([ratio_D5, ratio_S2], f"{iso}/U8", ["skyblue", "red"], grmin, grmax, [f"D5 : {iso} abs rates / U8 abs rates", f"S2 : {iso} abs rates / U8 abs rates"], save_dir = path_save, save_name=f"hist_{iso}_U8_ratio_absorption_rates_{input_case}{correlation}_{SSH_method}")

    print(f"ratio_D5 = {ratio_D5}")
    print(f"ratio_S2 = {ratio_S2}")

    error_ratio = (ratio_D5 - ratio_S2)*100/ratio_S2

    rms_error = np.sqrt(np.sum(error_ratio**2)/len(error_ratio))
    print(f"error_ratio = {error_ratio}")
    max_error = np.max(np.abs(error_ratio))
    print(f"rms_error = {rms_error}")
    print(f"For {iso}, ssh_method = {SSH_method} max error = {max_error}")
    print(f"max error = {max_error}")
    print(f"index of max error = {np.argmax(np.abs(error_ratio))}")
    print(f"ratio D5 at max error = {ratio_D5[np.argmax(np.abs(error_ratio))]}")
    print(f"ratio S2 at max error = {ratio_S2[np.argmax(np.abs(error_ratio))]}")

    plot_histograms([error_ratio], f"{iso}/U8", ["skyblue"], grmin, grmax, [f"Relative error on {iso}/U8 abs rate ratios : {SSH_method}"], save_dir = path_save, save_name=f"hist_{iso}_U8_error_on_ratio_absorption_rates_{input_case}{correlation}_{SSH_method}")

    return ratio_D5, ratio_S2, error_ratio


def normalize_absorption_rates(rates, GRMIN, GRMAX):
    """
    Normalize the absorption rates to the total absorption rate.
    """
    if len(rates) > (GRMAX - GRMIN) + 1:
        # extract the rates for the energy groups of interest
        rates = rates[GRMIN-1:GRMAX]
    rates = np.array(rates)/np.sum(rates)
    print(f"Normalized rates = {rates}")
    print(f"Sum of normalized rates = {np.sum(rates)}")
    return rates
    
def compute_group_relative_error(ratesS2, ratesD5, ngroups):
    """
    Compute the relative error on the absorption rates in each group.
    """
    relative_errors = []
    for i in range(ngroups):
        relative_error = (ratesD5[i]-ratesS2[i])*100/ratesS2[i]
        relative_errors.append(relative_error)

    print(f"length of relative errors = {len(relative_errors)}")
    return relative_errors

    

if __name__ == "__main__":
    input_case = "HOM_UOX_Gd157"
    correlation = "_CORR"
    path_S2 = f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{input_case}/XS_study"
    path_D5 = f"/home/p117902/working_dir/PolyBWR_project/Version5_wc/Dragon/Linux_aarch64/{input_case}_autop{correlation}.result"
    local_path = os.getcwd()
    library = "oldlib"
    bu_step = 0
    GRMIN = 52
    GRMAX = 206

    N_U238 = 7.765571e-03
    N_Gd157 = 7.899877e-05
    N_iso = N_Gd157

    #path_save=os.getcwd()
    path_exists_VDG = os.path.exists(f"{local_path}/VDG_errors/{input_case}")
    if not path_exists_VDG:
        os.makedirs(f"{local_path}/VDG_errors/{input_case}")

    path_exists = os.path.exists(f"{local_path}/{input_case}{correlation}")
    if not path_exists:
        os.makedirs(f"{local_path}/{input_case}{correlation}")


    path_save = f"{local_path}/{input_case}{correlation}"
    # Parse detector file for Gd157 absorptions (n,gamma) MT=102
    det_scores, n_ene_groups, n_react = parse_Serpent_detector(path_S2, library, bu_step)

    # recover reaction rates per group and reaction
    Gd_nGamma_rates = []
    U8_nGamma_rates = []
    U8_abs_rates = []
    for group in range(n_ene_groups):
        Gd_nGamma_rates.append(det_scores[group, 0])
        U8_nGamma_rates.append(det_scores[group, 1])
        U8_abs_rates.append(det_scores[group, 2])
    Gd_nGamma_rates_reversed = np.array(Gd_nGamma_rates[::-1])
    U8_nGamma_rates_reversed = np.array(U8_nGamma_rates[::-1])
    U8_abs_rates_reversed = np.array(U8_abs_rates[::-1])
    

    S2_rates = {"Gd157 (n,gamma)": Gd_nGamma_rates_reversed[GRMIN-1:GRMAX], "U238 (n,gamma)": U8_nGamma_rates_reversed[GRMIN-1:GRMAX], "U238 absorption": U8_abs_rates_reversed[GRMIN-1:GRMAX]}
    # Plot reaction rates
    # Histogram of reaction rates
    # plot_histogram(Gd_nGamma_rates.values())

    # Parse DRAGON absoprtion rates, these are microscopic XS rates
    D5_rates, D5_USSvsAUTO = parse_Dragon_abs_rates(path_D5, input_case, output_case = f"{input_case}_USS_AUTO_inrs1{correlation}", grmin=GRMIN, grmax=GRMAX)


    # Normalize the S2 absorption rates
    U8_abs_rates_norm = normalize_absorption_rates(U8_abs_rates_reversed, GRMIN, GRMAX)
    Gd_nGamma_rates_norm = normalize_absorption_rates(Gd_nGamma_rates_reversed, GRMIN, GRMAX)
    U8_nGamma_rates_norm = normalize_absorption_rates(U8_nGamma_rates_reversed, GRMIN, GRMAX)

    ERROR_RATIOS = {}

    for ssh_method in ["AUTO", "RSE", "PT", "SUBG"]:

        # Study U5/U8 absorption rates for both Serpent and Dragon
    
        # Study Gd157/U8 absorption rates for both Serpent and Dragon
        ratio_D5, ratio_S2, error_ratio = compare_iso_U8_absorption_rates(D5_rates, S2_rates, ssh_method, "Gd157", N_U238, N_iso, GRMIN, GRMAX, input_case, correlation, path_save, local_path)
        ERROR_RATIOS[f"Gd157/U8 {ssh_method}"] = error_ratio

        print(f"ratio_D5 = {ratio_D5}")
        print(f"ratio_S2 = {ratio_S2}")

        print(f"ratio of ratios = {ratio_D5/ratio_S2}")
        print(f"ratio of ratios scaled by N_U238/N_Gd157 = {ratio_D5/ratio_S2*(N_U238/N_Gd157)}")
        # Normalized D5/Autosecol absorption rates
        if ssh_method == "AUTO":
            U8_D5_norm = normalize_absorption_rates(D5_rates["AUTO U238"], GRMIN, GRMAX)
            Gd157_D5_norm = normalize_absorption_rates(D5_rates["AUTO Gd157"], GRMIN, GRMAX)
        else:
            U8_D5_norm = normalize_absorption_rates(D5_rates[f"USS U238, {ssh_method}"], GRMIN, GRMAX)
            Gd157_D5_norm = normalize_absorption_rates(D5_rates[f"USS Gd157, {ssh_method}"], GRMIN, GRMAX)

        

        # Study normalized Gd157, U8, U5 absorption rates for both Serpent and Dragon    
        # Plot normalized absorption rates
        plot_histograms([U8_abs_rates_norm, U8_nGamma_rates_norm, U8_D5_norm], "U238", ["skyblue", "red", "green"], GRMIN, GRMAX, legends=["S2 : U8 abs rates", "S2 : U8 nGamma rates", f"D5 : U8 abs rates {ssh_method}"], save_dir = path_save, save_name=f"U8_abs_rates_{ssh_method}")
        plot_histograms([Gd_nGamma_rates_norm, Gd157_D5_norm], "Gd157", ["skyblue", "red"], GRMIN, GRMAX, legends=["S2 : Gd157 nGamma rates", f"D5 : Gd157 abs rates {ssh_method}"], save_dir = path_save, save_name=f"Gd157_abs_rates_{ssh_method}")

        # compute relative error on abs rates in each group
        relative_errors_Gd157 = compute_group_relative_error(Gd_nGamma_rates_norm, Gd157_D5_norm, GRMAX-GRMIN+1)
        relative_errors_U8 = compute_group_relative_error(U8_nGamma_rates_norm, U8_D5_norm, GRMAX-GRMIN+1)
        plot_histograms([relative_errors_Gd157, relative_errors_U8], "Gd157 and U8", ["skyblue", "red"], GRMIN, GRMAX, legends=[f"{ssh_method}-S2 : relative errors Gd157", f"{ssh_method}-S2 : relative errors U8"], save_dir = path_save, save_name=f"relative_errors_{ssh_method}_Gd157_U8")


    # Plot the relative errors on the absorption rates ratios
        
    plot_histograms([ERROR_RATIOS["Gd157/U8 AUTO"], ERROR_RATIOS["Gd157/U8 RSE"], ERROR_RATIOS["Gd157/U8 PT"], ERROR_RATIOS["Gd157/U8 SUBG"]], "Gd157/U8", ["skyblue", "red", "green", "purple"], GRMIN, GRMAX, legends=["error AUTO-S2", "error RSE-S2", "error PT-S2", "error SUBG-S2"], save_dir = path_save, save_name=f"hist_relative_errors_abs_rates_ratios_Gd157_U8_ALL")