# PyGan script to post-treat results for HOM_UOX_Gd15* calculations
# Author: R. Guasch
# Date: 2024-11-12
# Purpose: Post-treat results for HOM_UOX_Gd155/157 calculations :
#         - Extract absorption rates and cross sections from COMPO_MESHES multi-compo 
#         - Extract absorption rates and cross sections from Serpent2 reference
#
# In COMP_MESHES, the results are tabulated according to : 
#           - the energy mesh used, SHEM281/295/315 identified in EDIHOM_xxx directories
#           - the self-shielding method used, identified through the 'SSH' parameters : A, R, P, S for Autosecol, RSE, PT and SUBG respectively

# In Serpent2, the results for absorption rates are obtained from detector results for MT=101 and MT=102 reactions of Gd157/U238
# The cross sections are obtained from microdepletion output files
# Results are obtained for PyNjoy2016 and oldlib librairies, only on the SHEM295 energy mesh for now.
#

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import lifo
import lcm
import cle2000
from MeshHandler import energyMeshHandler as ENEMESH
from PT_Gd157 import postTreatment_rates_XS_D5 as PT_D5
from PT_Gd157 import postTreatment_rates_XS_S2 as PT_S2
from compareD5_S2 import compare_D5_S2_rates_XS as CD5S2

def generate_latex_table(library, mesh, data):
    """
    Generate LaTeX code for a table showing SSH results for a given library and mesh.

    Parameters:
    - library (str): Name of the library (e.g., "PyNjoy2016", "oldlib").
    - mesh (str): Mesh name (e.g., "SHEM295").
    - data (list of dict): List of dictionaries containing the table rows. 
      Each dictionary should have the following keys:
      'SSH Method', 'Error (%)', 'Group', 'u_min', 'u_max', 'E_min', 'E_max'.

    Returns:
    - str: LaTeX code for the table.
    """
    table_header = f"""
\\section*{{Results for {library} Library on {mesh} Mesh}}
\\begin{{table}}[h!]
\\centering
\\begin{{tabular}}{{@{{}}lrrrrrr@{{}}}}
\\toprule
\\textbf{{SSH Method}} & \\textbf{{Error (\\%)}} & \\textbf{{Group}} & \\textbf{{u\\_min}} & \\textbf{{u\\_max}} & \\textbf{{E\\_min (eV)}} & \\textbf{{E\\_max (eV)}} \\\\ \\midrule
"""
    table_rows = ""
    for row in data:
        table_rows += (
            f"{row['SSH Method']} & {row['Error (%)']:.3f} & {row['Group']} & "
            f"{row['u_min']:.3f} & {row['u_max']:.3f} & {row['E_min']:.3f} & {row['E_max']:.3f} \\\\\n"
        )

    table_footer = """\\bottomrule
\\end{tabular}
\\caption{Results for different SSH methods with %s on the %s mesh.}
\\label{tab:%s_%s}
\\end{table}
""" % (library, mesh, library.lower(), mesh.lower())

    return table_header + table_rows + table_footer

# Import COMPO object
path = os.getcwd()
name_compo='COMPO_MESHES_noCORR'
os.chdir("Gd157_rates_XS_proc")
pyCOMPO=lcm.new('LCM_INP',name_compo,impx=0)
os.chdir(path)
# Creation of results directory
path=os.getcwd()
save_path = f"Gd157_Rates_and_XS_results_PyGan"
save_path_comparison = f"{save_path}/comparison_D5_S2"
a=os.path.exists(save_path)
if a==False:
	os.mkdir(save_path)

a=os.path.exists(save_path_comparison)
if a==False:
	os.mkdir(save_path_comparison)
case_name = "HOM_UOX_Gd157"

meshes = ["XMAS172", "SHEM281","SHEM295","SHEM315"]
SSH_methods = {"A":"Autosecol","R":"RSE","P":"PT","S":"SUBG"}
Gd157_ngamma = {"XMAS172":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}},
                "SHEM281":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}},
                "SHEM295":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}},
                "SHEM315":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}}} # build a dictionary of dictionaries to store the rates for each mesh and SSH method
Gd157_abs =  {"XMAS172":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}},
                "SHEM281":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}},
                "SHEM295":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}},
                "SHEM315":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}}}
U238_ngamma = {"XMAS172":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}},
                "SHEM281":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}},
                "SHEM295":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}},
                "SHEM315":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}}}
mesh_objects_D5 = {"XMAS172":None,"SHEM281":None,"SHEM295":None,"SHEM315":None}
mesh_objects_S2 = {"SHEM281":None,"SHEM295":None,"SHEM315":None}

fluxes = {"XMAS172":{"Autosecol":[], "RSE":[], "PT":[], "SUBG":[]}, "SHEM281":{"Autosecol":[], "RSE":[], "PT":[], "SUBG":[]},"SHEM295":{"Autosecol":[], "RSE":[], "PT":[], "SUBG":[]},"SHEM315":{"Autosecol":[], "RSE":[], "PT":[], "SUBG":[]}}

for mesh in meshes:
    if mesh == "SHEM281":
        DIR = "EDIHOM_281"
    elif mesh == "SHEM295":
        DIR = "EDIHOM_295"
    elif mesh == "SHEM315":
        DIR = "EDIHOM_315"
    elif mesh == "XMAS172":
        DIR = "EDIHOM_172"
    energyMESH = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][0]["ENERGY"]
    print(f"The {mesh} mesh has {len(energyMESH)} energy bounds, so {len(energyMESH)-1} energy groups")
    MESH_obj = ENEMESH(mesh,energyMESH,1.0E+07,"eV")
    MESH_obj.printnfgCard()
    mesh_objects_D5[mesh] = MESH_obj
    if mesh in ["SHEM281","SHEM295","SHEM315"]:
        mesh_objects_S2[mesh] = MESH_obj

    SSH_methods_par = pyCOMPO[DIR]["GLOBAL"]['pval00000001'].strip().split("           ")
    for i in range(len(SSH_methods_par)):
        print(f"SSH method {i} = {SSH_methods[SSH_methods_par[i]]}")
        isotope = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][1]['ALIAS'][0:5]
        print(f"isotope = {isotope}")
        print(f"number of calculations = {len(pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'])}")

        N_U238 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESDENS'][0]
        N_Gd157 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESDENS'][1]

        XS_NGAMMA_U238 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][0]['NG'] 
        XS_NGAMMA_Gd157 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][1]['NG'] 

        PHI = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][0]['NWT0']  

        # Reconstruct reaction rates
        NGAMMA_Gd157_rates = np.array(XS_NGAMMA_Gd157)*np.array(PHI)*N_Gd157
        NGAMMA_U238_rates = np.array(XS_NGAMMA_U238)*np.array(PHI)*N_U238

        # Store the results in the dictionaries

        Gd157_ngamma[mesh][SSH_methods[SSH_methods_par[i]]]["rates"] = NGAMMA_Gd157_rates
        U238_ngamma[mesh][SSH_methods[SSH_methods_par[i]]]["rates"] = NGAMMA_U238_rates

        Gd157_ngamma[mesh][SSH_methods[SSH_methods_par[i]]]["XS"] = XS_NGAMMA_Gd157
        U238_ngamma[mesh][SSH_methods[SSH_methods_par[i]]]["XS"] = XS_NGAMMA_U238

        fluxes[mesh][SSH_methods[SSH_methods_par[i]]] = PHI


# Call post treatment class

DRAGON5_case = PT_D5(case_name, mesh_objects_D5, SSH_methods.values(), save_path)
DRAGON5_case.set_reaction_data("Gd157_ngamma",Gd157_ngamma)
DRAGON5_case.set_reaction_data("U238_ngamma",U238_ngamma)
DRAGON5_case.set_fluxes(fluxes)

# Plot the results

for mesh in meshes:
    DRAGON5_case.plot_fluxes(mesh)
    DRAGON5_case.plot_XS(mesh,"Gd157_ngamma")
    DRAGON5_case.plot_XS(mesh,"U238_ngamma")
    DRAGON5_case.plot_reaction_rates(mesh,"Gd157_ngamma")
    DRAGON5_case.plot_reaction_rates(mesh,"U238_ngamma")
    if mesh == "SHEM281":
        grmin = 1
        grmax = 281
    elif mesh == "SHEM295":
        grmin = 1
        grmax = 295
    elif mesh == "SHEM315":
        grmin = 1
        grmax = 315
    elif mesh == "XMAS172":
        grmin = 1
        grmax = 172
    D5_SSH_methods = ["RSE","PT","SUBG"]
    if mesh == "XMAS172":
        D5_SSH_methods = ["RSE","SUBG"] # not enough groups for PT : 250 groups required
    DRAGON5_case.compute_relative_differences_XS("Gd157_ngamma",mesh, "Autosecol", D5_SSH_methods)
    DRAGON5_case.compute_relative_differences_XS("U238_ngamma",mesh, "Autosecol", D5_SSH_methods)
    DRAGON5_case.compute_relative_differences_Rates("Gd157_ngamma",mesh, "Autosecol", D5_SSH_methods)
    DRAGON5_case.compute_relative_differences_Rates("U238_ngamma",mesh, "Autosecol", D5_SSH_methods)

    DRAGON5_case.plot_histogram_relative_differences_XS("Gd157_ngamma",mesh, D5_SSH_methods, grmin,grmax)
    DRAGON5_case.plot_relative_differences_XS("Gd157_ngamma", mesh, D5_SSH_methods)
    DRAGON5_case.plot_histogram_relative_differences_XS("U238_ngamma", mesh, D5_SSH_methods,grmin, grmax)
    DRAGON5_case.plot_relative_differences_XS("U238_ngamma", mesh, D5_SSH_methods)

    DRAGON5_case.plot_histogram_relative_differences_Rates("Gd157_ngamma", mesh, D5_SSH_methods, grmin, grmax)
    DRAGON5_case.plot_relative_differences_Rates("Gd157_ngamma", mesh, D5_SSH_methods)
    DRAGON5_case.plot_histogram_relative_differences_Rates("U238_ngamma", mesh, D5_SSH_methods, grmin, grmax)
    DRAGON5_case.plot_relative_differences_Rates("U238_ngamma", mesh, D5_SSH_methods)
    
    for method in D5_SSH_methods:
        DRAGON5_case.plot_histogram_relative_differences_XS("Gd157_ngamma",mesh, [method],grmin,grmax)
        DRAGON5_case.plot_histogram_relative_differences_XS("U238_ngamma",mesh, [method],grmin,grmax)
        DRAGON5_case.plot_histogram_relative_differences_Rates("Gd157_ngamma",mesh, [method],grmin,grmax)
        DRAGON5_case.plot_histogram_relative_differences_Rates("U238_ngamma",mesh, [method],grmin,grmax)

# Post-treatment for SERPENT2 results
# Make a case disjunction for XS and rates : due to amount of time/no necessity to obtain XS at all BU steps : 
# only BU=0 is considered for XS
# All BU are considered for rates 
SERPENT2_case = PT_S2(case_name, mesh_objects_S2, ["PyNjoy2016", "oldlib"], range(0,1), save_path)
SERPENT2_case.parse_S2_outputs()

for mesh_name in ["SHEM281","SHEM295","SHEM315"]: 
    # Plot XS and reaction rates at BU=0 for PyNjoy2016 and oldlib on SHEM281/295/315
    # not all bu steps are available for all meshes with oldlib, it took too long to run the simulations for oldlib
    SERPENT2_case.plot_XS(mesh_name,"Gd157_ngamma",bu_step=0)
    SERPENT2_case.plot_XS(mesh_name,"U238_ngamma",bu_step=0)
    SERPENT2_case.plot_reaction_rates(mesh_name,"Gd157_ngamma",bu_step=0)
    SERPENT2_case.plot_reaction_rates(mesh_name,"U238_ngamma",bu_step=0)

    # Compute relative differences
    SERPENT2_case.compute_relative_differences_XS("Gd157_ngamma", mesh_name, "PyNjoy2016", "oldlib")
    SERPENT2_case.compute_relative_differences_XS("U238_ngamma", mesh_name, "PyNjoy2016", "oldlib")
    SERPENT2_case.compute_relative_differences_Rates("Gd157_ngamma", mesh_name, "PyNjoy2016", "oldlib")
    SERPENT2_case.compute_relative_differences_Rates("U238_ngamma", mesh_name, "PyNjoy2016", "oldlib")

    # Plot relative differences

    SERPENT2_case.plot_relative_differences_XS("Gd157_ngamma", mesh_name, bu_step = 0)
    SERPENT2_case.plot_relative_differences_XS("U238_ngamma", mesh_name, bu_step = 0)
    SERPENT2_case.plot_relative_differences_Rates("Gd157_ngamma", mesh_name, bu_step = 0)
    SERPENT2_case.plot_relative_differences_Rates("U238_ngamma", mesh_name, bu_step = 0)

comparison_D5_S2 = CD5S2("HOM_UOX_Gd157_rates_XS_study", DRAGON5_case, SERPENT2_case, save_path_comparison)
all_ssh_methods = ["Autosecol","RSE","PT","SUBG"]
diff_data_rates = {"PyNjoy2016":{}, "oldlib":{}}
diff_data_XS = {"PyNjoy2016":{}, "oldlib":{}}
for mesh_name in ["SHEM281","SHEM295","SHEM315"]:
    comparison_D5_S2.plot_rates_D5_S2("Gd157_ngamma", mesh_name, all_ssh_methods)
    comparison_D5_S2.plot_XS_D5_S2("Gd157_ngamma", mesh_name, all_ssh_methods)
    # renormalize reaction rates to the same value for both D5 and S2
    comparison_D5_S2.renorm_rates("Gd157_ngamma", mesh_name)
    comparison_D5_S2.compare_reaction_rates_and_XS("Gd157_ngamma", mesh_name, all_ssh_methods)

    # Plot XS and rates for all ssh methods on the same plot
    comparison_D5_S2.plot_diff_rates_D5_S2("Gd157_ngamma", mesh_name, "PyNjoy2016", all_ssh_methods)
    comparison_D5_S2.plot_diff_rates_D5_S2("Gd157_ngamma", mesh_name, "oldlib", all_ssh_methods)
    comparison_D5_S2.plot_diff_XS_D5_S2("Gd157_ngamma", mesh_name, "PyNjoy2016", all_ssh_methods)
    comparison_D5_S2.plot_diff_XS_D5_S2("Gd157_ngamma", mesh_name, "oldlib", all_ssh_methods)

    # plot on each individual plot : 1 per ssh method
    for ssh_method in all_ssh_methods:
        comparison_D5_S2.plot_diff_rates_D5_S2("Gd157_ngamma", mesh_name, "PyNjoy2016", [ssh_method])
        comparison_D5_S2.plot_diff_rates_D5_S2("Gd157_ngamma", mesh_name, "oldlib", [ssh_method])
        comparison_D5_S2.plot_diff_XS_D5_S2("Gd157_ngamma", mesh_name, "PyNjoy2016", [ssh_method])
        comparison_D5_S2.plot_diff_XS_D5_S2("Gd157_ngamma", mesh_name, "oldlib", [ssh_method])


    # find top differences
    diff_data_rates["PyNjoy2016"][mesh_name] = comparison_D5_S2.find_top_differences_rates("PyNjoy2016", "Gd157_ngamma", mesh_name, 5)
    diff_data_rates["oldlib"][mesh_name] = comparison_D5_S2.find_top_differences_rates("oldlib", "Gd157_ngamma", mesh_name, 5)
    diff_data_XS["PyNjoy2016"][mesh_name] = comparison_D5_S2.find_top_differences_XS("PyNjoy2016", "Gd157_ngamma", mesh_name, 5)
    diff_data_XS["oldlib"][mesh_name] = comparison_D5_S2.find_top_differences_XS("oldlib", "Gd157_ngamma", mesh_name, 5)


# Analyse the results
# What needs to be checked and analyzed : 
    # - What are the most important differences in the reaction rates and cross sections
    # - Are they consistent across the different energy meshes : ie where are they located on the energy spectrum
    # - Are differences on rates and XS fully correlated : if so can focus on XS only 
    # - Are differences on rates and XS consistent across the different self-shielding methods
    # - Are differences on rates and XS consistent across the different S2 libaries

# Generate LaTeX tables for the top differences
for mesh_name in ["SHEM281","SHEM295","SHEM315"]:
    for library in ["PyNjoy2016", "oldlib"]:
        diff_data_rates_table = generate_latex_table(library, mesh_name, diff_data_rates[library][mesh_name])
        diff_data_XS_table = generate_latex_table(library, mesh_name, diff_data_XS[library][mesh_name])

        with open(f"{save_path_comparison}/top_diff_rates_{library}_{mesh_name}.tex", "w") as f:
            f.write(diff_data_rates_table)
        
        with open(f"{save_path_comparison}/top_diff_XS_{library}_{mesh_name}.tex", "w") as f:
            f.write(diff_data_XS_table)

print("Post-treatment done")

