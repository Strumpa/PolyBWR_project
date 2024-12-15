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
import re
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import lifo
import lcm
import cle2000
from MeshHandler import energyMeshHandler as ENEMESH
from PT_D5S2 import postTreatment_rates_XS_D5 as PT_D5
from PT_D5S2 import postTreatment_rates_XS_S2 as PT_S2
from D5_vs_S2 import compare_D5_S2_rates_XS as CD5S2

def generate_latex_table_EMESH(library, mesh, data):
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

def generate_latex_table_autolib(library, autolib, data):
    """
    Generate LaTeX code for a table showing SSH results for a given library and mesh.

    Parameters:
    - library (str): Name of the Serpent2 library (e.g., "PyNjoy2016", "oldlib").
    - autolib (str): Autolib name (e.g., "J311_SHEM295", "J311_SHEM295F" ...).
    - data (list of dict): List of dictionaries containing the table rows. 
      Each dictionary should have the following keys:
      'SSH Method', 'Error (%)', 'Group', 'u_min', 'u_max', 'E_min', 'E_max'.

    Returns:
    - str: LaTeX code for the table.
    """
    table_header = f"""
\\section*{{Results for S2 {library} Library vs D5 on {autolib}}}
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
""" % (library, autolib, library.lower(), autolib.lower())

    return table_header + table_rows + table_footer

def generate_latex_table_IRSET(library, IRSET, data):
    """
    Generate LaTeX code for a table showing SSH results for a given library and mesh.

    Parameters:
    - library (str): Name of the Serpent2 library (e.g., "PyNjoy2016", "oldlib").
    - IRSET (str): IRSET option name (e.g., "IRSET 101", ).
    - data (list of dict): List of dictionaries containing the table rows. 
      Each dictionary should have the following keys:
      'SSH Method', 'Error (%)', 'Group', 'u_min', 'u_max', 'E_min', 'E_max'.

    Returns:
    - str: LaTeX code for the table.
    """
    table_header = f"""
\\section*{{Results for S2 {library} Library vs D5 with {IRSET}}}
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
""" % (library, IRSET, library.lower(), IRSET.lower())

    return table_header + table_rows + table_footer



def post_treat_HOM_UOX_Gd157_MESHES(correlation_model = False):
    # Path to Serpent2 results
    path_to_S2_results = "/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/HOM_UOX_Gd157/XS_rates_study"
    # Import COMPO object
    path = os.getcwd()
    if correlation_model:
        name_compo='COMPO_MESHES_CORR'
    else:
        name_compo='COMPO_MESHES_noCORR'
    os.chdir("Gd157_rates_XS_proc")
    pyCOMPO=lcm.new('LCM_INP',name_compo,impx=0)
    os.chdir(path)
    # Creation of results directory
    path=os.getcwd()
    save_path = f"Gd157_Rates_and_XS_results_PyGan/EnergyMeshes_study"
    save_path_comparison = f"{save_path}/comparison_D5_S2"
    a=os.path.exists(save_path)
    if a==False:
        os.makedirs(save_path)

    a=os.path.exists(save_path_comparison)
    if a==False:
        os.makedirs(save_path_comparison)
    if correlation_model:
        case_name = "HOM_UOX_Gd157_CORR"
    else:
        case_name = "HOM_UOX_Gd157_noCORR"

    meshes = ["XMAS172", "SHEM281","SHEM295","SHEM315"]
    SSH_methods = {"A":"Autosecol","R":"RSE","P":"PT","S":"SUBG"}
    Gd157_ngamma = {}
    U238_ngamma = {}

    mesh_objects = {"XMAS172":None,"SHEM281":None,"SHEM295":None,"SHEM315":None}

    fluxes = {}

    for mesh in meshes:
        Gd157_ngamma[mesh] = {}
        U238_ngamma[mesh] = {}
        fluxes[mesh] = {}
        if mesh == "XMAS172":
            DIR = "EDIHOM_172"
        elif mesh == "SHEM281":
            DIR = "EDIHOM_281"
        elif mesh == "SHEM295":
            DIR = "EDIHOM_295"
        elif mesh == "SHEM315":
            DIR = "EDIHOM_315"
        energyMESH = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][0]["ENERGY"]
        print(f"The {mesh} mesh has {len(energyMESH)} energy bounds, so {len(energyMESH)-1} energy groups")
        MESH_obj = ENEMESH(mesh,energyMESH,1.0E+07,"eV")
        #MESH_obj.printnfgCard()
        mesh_objects[mesh] = MESH_obj

        SSH_methods_par = re.findall(r"[A-Za-z]+", pyCOMPO[DIR]["GLOBAL"]['pval00000001'].strip())
        for i in range(len(SSH_methods_par)):
            Gd157_ngamma[mesh][SSH_methods[SSH_methods_par[i]]] = {"rates":[], "XS":[]}
            U238_ngamma[mesh][SSH_methods[SSH_methods_par[i]]] = {"rates":[], "XS":[]}
            fluxes[mesh][SSH_methods[SSH_methods_par[i]]] = []
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
            if mesh == "XMAS172":
                if len(Gd157_ngamma[mesh][SSH_methods[SSH_methods_par[i]]]['XS']) != 172 or len(U238_ngamma[mesh][SSH_methods[SSH_methods_par[i]]]['XS']) != 172:
                    print(f"Error : XS data for {mesh} and {SSH_methods[SSH_methods_par[i]]} is not complete")
            if mesh == "SHEM281":
                if len(Gd157_ngamma[mesh][SSH_methods[SSH_methods_par[i]]]['XS']) != 281 or len(U238_ngamma[mesh][SSH_methods[SSH_methods_par[i]]]['XS']) != 281:
                    print(f"Error : XS data for {mesh} and {SSH_methods[SSH_methods_par[i]]} is not complete")
            if mesh == "SHEM295":
                if len(Gd157_ngamma[mesh][SSH_methods[SSH_methods_par[i]]]['XS']) != 295 or len(U238_ngamma[mesh][SSH_methods[SSH_methods_par[i]]]['XS']) != 295:
                    print(f"Error : XS data for {mesh} and {SSH_methods[SSH_methods_par[i]]} is not complete")
            if mesh == "SHEM315":
                if len(Gd157_ngamma[mesh][SSH_methods[SSH_methods_par[i]]]['XS']) != 315 or len(U238_ngamma[mesh][SSH_methods[SSH_methods_par[i]]]['XS']) != 315:
                    print(f"Error : XS data for {mesh} and {SSH_methods[SSH_methods_par[i]]} is not complete")

    # Call post treatment class
    
    DRAGON5_EMESH = PT_D5(case_name, "EMESH", mesh_objects, ["Gd157_ngamma", "U238_ngamma"], meshes, SSH_methods.values(), save_path)
    DRAGON5_EMESH.set_reaction_data("Gd157_ngamma",Gd157_ngamma)
    DRAGON5_EMESH.set_reaction_data("U238_ngamma",U238_ngamma)
    DRAGON5_EMESH.set_fluxes(fluxes)

    # Plot the results

    for mesh_name in meshes:
        print(f"Plotting results for mesh_name {mesh_name}")
        DRAGON5_EMESH.plot_fluxes(mesh_name, mesh_name)
        DRAGON5_EMESH.plot_XS(mesh_name,mesh_name,"Gd157_ngamma")
        DRAGON5_EMESH.plot_XS(mesh_name,mesh_name,"U238_ngamma")
        DRAGON5_EMESH.plot_reaction_rates(mesh_name,mesh_name,"Gd157_ngamma")
        DRAGON5_EMESH.plot_reaction_rates(mesh_name,mesh_name,"U238_ngamma")
        if mesh_name == "SHEM281":
            grmin = 1
            grmax = 281
        elif mesh_name == "SHEM295":
            grmin = 1
            grmax = 295
        elif mesh_name == "SHEM315":
            grmin = 1
            grmax = 315
        elif mesh_name == "XMAS172":
            grmin = 1
            grmax = 172
        D5_SSH_methods = ["RSE","PT","SUBG"]
        if mesh_name == "XMAS172":
            D5_SSH_methods = ["RSE","SUBG"] # not enough groups for PT : 250 groups required
        DRAGON5_EMESH.compute_relative_differences_XS("Gd157_ngamma", mesh_name, "Autosecol", D5_SSH_methods)
        DRAGON5_EMESH.compute_relative_differences_XS("U238_ngamma", mesh_name, "Autosecol", D5_SSH_methods)
        DRAGON5_EMESH.compute_relative_differences_Rates("Gd157_ngamma", mesh_name, "Autosecol", D5_SSH_methods)
        DRAGON5_EMESH.compute_relative_differences_Rates("U238_ngamma", mesh_name, "Autosecol", D5_SSH_methods)

        DRAGON5_EMESH.plot_histogram_relative_differences_XS("Gd157_ngamma", mesh_name, D5_SSH_methods, grmin, grmax)
        DRAGON5_EMESH.plot_relative_differences_XS("Gd157_ngamma", mesh_name, mesh_name, D5_SSH_methods)
        DRAGON5_EMESH.plot_histogram_relative_differences_XS("U238_ngamma", mesh_name, D5_SSH_methods, grmin, grmax)
        DRAGON5_EMESH.plot_relative_differences_XS("U238_ngamma", mesh_name, mesh_name, D5_SSH_methods)

        DRAGON5_EMESH.plot_histogram_relative_differences_Rates("Gd157_ngamma", mesh_name, D5_SSH_methods, grmin, grmax)
        DRAGON5_EMESH.plot_relative_differences_Rates("Gd157_ngamma", mesh_name, mesh_name, D5_SSH_methods)
        DRAGON5_EMESH.plot_histogram_relative_differences_Rates("U238_ngamma", mesh_name, D5_SSH_methods, grmin, grmax)
        DRAGON5_EMESH.plot_relative_differences_Rates("U238_ngamma", mesh_name, mesh_name, D5_SSH_methods)
        
        for method in D5_SSH_methods:
            DRAGON5_EMESH.plot_relative_differences_XS("Gd157_ngamma", mesh_name, mesh_name, [method])
            DRAGON5_EMESH.plot_relative_differences_XS("U238_ngamma", mesh_name, mesh_name, [method])
            DRAGON5_EMESH.plot_relative_differences_Rates("Gd157_ngamma", mesh_name, mesh_name, [method])
            DRAGON5_EMESH.plot_relative_differences_Rates("U238_ngamma", mesh_name, mesh_name, [method])

    # Post-treatment for SERPENT2 results
    # Make a case disjunction for XS and rates : due to amount of time/no necessity to obtain XS at all BU steps : 
    # only BU=0 is considered for XS
    # All BU are considered for rates 
    S2_case_name = "HOM_UOX_Gd157"
    SERPENT2_case = PT_S2(S2_case_name, mesh_objects, ["PyNjoy2016", "oldlib"], range(0,1), save_path)
    SERPENT2_case.parse_S2_outputs(path_to_S2_results)

    for mesh_name in meshes: 
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
    # (case_name, study_type, D5_case, S2_case, S2_libs, compo_keywords, isotopes, self_shielding_methods, save_path)
    comparison_D5_S2 = CD5S2("HOM_UOX_Gd157_rates_XS_study", "EMESH", DRAGON5_EMESH, SERPENT2_case, 
                            S2_libs = ["oldlib", "PyNjoy2016"], compo_keywords = meshes, isotopes = ["Gd157", "U238"], self_shielding_methods = SSH_methods.values(),
                            save_path =   save_path_comparison)
    all_ssh_methods = ["Autosecol","RSE","PT","SUBG"]
    diff_data_rates = {"PyNjoy2016":{}, "oldlib":{}}
    diff_data_XS = {"PyNjoy2016":{}, "oldlib":{}}
    for mesh_name in meshes:
        if mesh_name == "XMAS172":
            all_ssh_methods = ["Autosecol","RSE","SUBG"] # not enough groups for PT : 250 groups required
        else:
            all_ssh_methods = ["Autosecol","RSE","PT","SUBG"]
        comparison_D5_S2.plot_rates_D5_S2("Gd157_ngamma", mesh_name, mesh_name, all_ssh_methods)
        comparison_D5_S2.plot_XS_D5_S2("Gd157_ngamma", mesh_name, mesh_name, all_ssh_methods)
        # renormalize reaction rates to the same value for both D5 and S2
        comparison_D5_S2.renorm_rates("Gd157_ngamma", mesh_name, mesh_name)
        comparison_D5_S2.compare_reaction_rates_and_XS("Gd157_ngamma", mesh_name, all_ssh_methods)

        # Plot XS and rates for all ssh methods on the same plot
        comparison_D5_S2.plot_diff_rates_D5_S2("Gd157_ngamma", mesh_name, mesh_name, "PyNjoy2016", all_ssh_methods)
        comparison_D5_S2.plot_diff_rates_D5_S2("Gd157_ngamma", mesh_name, mesh_name, "oldlib", all_ssh_methods)
        comparison_D5_S2.plot_diff_XS_D5_S2("Gd157_ngamma", mesh_name, mesh_name, "PyNjoy2016", all_ssh_methods)
        comparison_D5_S2.plot_diff_XS_D5_S2("Gd157_ngamma", mesh_name, mesh_name, "oldlib", all_ssh_methods)

        # plot on each individual plot : 1 per ssh method
        for ssh_method in all_ssh_methods:
            comparison_D5_S2.plot_diff_rates_D5_S2("Gd157_ngamma", mesh_name, mesh_name, "PyNjoy2016", [ssh_method])
            comparison_D5_S2.plot_diff_rates_D5_S2("Gd157_ngamma", mesh_name, mesh_name, "oldlib", [ssh_method])
            comparison_D5_S2.plot_diff_XS_D5_S2("Gd157_ngamma", mesh_name, mesh_name, "PyNjoy2016", [ssh_method])
            comparison_D5_S2.plot_diff_XS_D5_S2("Gd157_ngamma", mesh_name, mesh_name, "oldlib", [ssh_method])


        # find top differences
        diff_data_rates["PyNjoy2016"][mesh_name] = comparison_D5_S2.find_top_differences_rates("PyNjoy2016", "Gd157_ngamma", mesh_name, mesh_name, 5)
        diff_data_rates["oldlib"][mesh_name] = comparison_D5_S2.find_top_differences_rates("oldlib", "Gd157_ngamma", mesh_name, mesh_name, 5)
        diff_data_XS["PyNjoy2016"][mesh_name] = comparison_D5_S2.find_top_differences_XS("PyNjoy2016", "Gd157_ngamma", mesh_name, mesh_name, 5)
        diff_data_XS["oldlib"][mesh_name] = comparison_D5_S2.find_top_differences_XS("oldlib", "Gd157_ngamma", mesh_name, mesh_name, 5)


    # Analyse the results
    # What needs to be checked and analyzed : 
        # - What are the most important differences in the reaction rates and cross sections
        # - Are they consistent across the different energy meshes : ie where are they located on the energy spectrum
        # - Are differences on rates and XS fully correlated : if so can focus on XS only 
        # - Are differences on rates and XS consistent across the different self-shielding methods
        # - Are differences on rates and XS consistent across the different S2 libaries

    # Generate LaTeX tables for the top differences
    for mesh_name in meshes:
        for library in ["PyNjoy2016", "oldlib"]:
            diff_data_rates_table = generate_latex_table_EMESH(library, mesh_name, diff_data_rates[library][mesh_name])
            diff_data_XS_table = generate_latex_table_EMESH(library, mesh_name, diff_data_XS[library][mesh_name])

            with open(f"{save_path_comparison}/top_diff_rates_{library}_{mesh_name}.tex", "w") as f:
                f.write(diff_data_rates_table)
            
            with open(f"{save_path_comparison}/top_diff_XS_{library}_{mesh_name}.tex", "w") as f:
                f.write(diff_data_XS_table)

    print("Post-treatment of HOM_UOX_Gd157 Energy meshes study done")
    return


def post_treat_HOM_UOX_Gd157_AUTOlib():

    # Path to Serpent2 results
    path_to_S2_results = "/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/HOM_UOX_Gd157/XS_rates_study"
    # Import COMPO object
    path = os.getcwd()
    name_compo='COMPO_295_AUTOLIB'
    os.chdir("Gd157_rates_XS_proc")
    pyCOMPO=lcm.new('LCM_INP',name_compo,impx=0)
    os.chdir(path)
    # Creation of results directory
    path=os.getcwd()
    save_path = f"Gd157_Rates_and_XS_results_PyGan/AUTOlib_study"
    save_path_comparison = f"{save_path}/comparison_D5_S2"
    a=os.path.exists(save_path)
    if a==False:
        os.makedirs(save_path)

    a=os.path.exists(save_path_comparison)
    if a==False:
        os.makedirs(save_path_comparison)
    case_name = "HOM_UOX_Gd157"
    meshes = ["SHEM295"]
    AUTOlibs = ["J311_295", "J311_295F", "J311_295FF"]
    SSH_methods = {"A":"Autosecol","R":"RSE","P":"PT","S":"SUBG"}
    D5_SSH_methods = ["RSE","PT","SUBG"]
    Gd157_ngamma = {}
    U238_ngamma = {}
    fluxes = {}
    mesh_objects = {"SHEM295":None}
    for autolib in AUTOlibs:
        Gd157_ngamma[autolib] = {}
        U238_ngamma[autolib] = {}
        fluxes[autolib] = {}
        if autolib == "J311_295":
            DIR = "EDIHOM_295"
        elif autolib == "J311_295F":
            DIR = "EDIHOM_295F"
        elif autolib == "J311_295FF":
            DIR = "EDIHOM_295FF"
        energyMESH = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][0]["ENERGY"]
        print(f"The {autolib} mesh has {len(energyMESH)} energy bounds, so {len(energyMESH)-1} energy groups")
        MESH_obj = ENEMESH("SHEM295",energyMESH,1.0E+07,"eV")
        mesh_objects["SHEM295"] = MESH_obj
        SSH_methods_par = re.findall(r"[A-Za-z]+", pyCOMPO[DIR]["GLOBAL"]['pval00000001'].strip())
        for i in range(len(SSH_methods_par)):
            Gd157_ngamma[autolib][SSH_methods[SSH_methods_par[i]]] = {"rates":[], "XS":[]}
            U238_ngamma[autolib][SSH_methods[SSH_methods_par[i]]] = {"rates":[], "XS":[]}
            fluxes[autolib][SSH_methods[SSH_methods_par[i]]] = []

            N_U238 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESDENS'][0]
            N_Gd157 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESDENS'][1]

            XS_NGAMMA_U238 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][0]['NG'] 
            XS_NGAMMA_Gd157 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][1]['NG'] 

            PHI = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][0]['NWT0']  
            # Reconstruct reaction rates
            NGAMMA_Gd157_rates = np.array(XS_NGAMMA_Gd157)*np.array(PHI)*N_Gd157
            NGAMMA_U238_rates = np.array(XS_NGAMMA_U238)*np.array(PHI)*N_U238

            # Store the results in the dictionaries

            Gd157_ngamma[autolib][SSH_methods[SSH_methods_par[i]]]["rates"] = NGAMMA_Gd157_rates
            U238_ngamma[autolib][SSH_methods[SSH_methods_par[i]]]["rates"] = NGAMMA_U238_rates

            Gd157_ngamma[autolib][SSH_methods[SSH_methods_par[i]]]["XS"] = XS_NGAMMA_Gd157
            U238_ngamma[autolib][SSH_methods[SSH_methods_par[i]]]["XS"] = XS_NGAMMA_U238

            fluxes[autolib][SSH_methods[SSH_methods_par[i]]] = PHI

    # Call post treatment class
    DRAGON5_AUTOLIB = PT_D5(case_name, "AUTOlib", mesh_objects, ["Gd157_ngamma", "U238_ngamma"], AUTOlibs, SSH_methods.values(), save_path)
    DRAGON5_AUTOLIB.set_reaction_data("Gd157_ngamma",Gd157_ngamma)
    DRAGON5_AUTOLIB.set_reaction_data("U238_ngamma",U238_ngamma)
    DRAGON5_AUTOLIB.set_fluxes(fluxes)

    # Plot the results
    for autolib in AUTOlibs:
        DRAGON5_AUTOLIB.plot_fluxes(autolib, "SHEM295")
        DRAGON5_AUTOLIB.plot_XS(autolib, "SHEM295", "Gd157_ngamma")
        DRAGON5_AUTOLIB.plot_XS(autolib, "SHEM295", "U238_ngamma")
        DRAGON5_AUTOLIB.plot_reaction_rates(autolib, "SHEM295", "Gd157_ngamma")
        DRAGON5_AUTOLIB.plot_reaction_rates(autolib, "SHEM295", "U238_ngamma")
        DRAGON5_AUTOLIB.compute_relative_differences_XS("Gd157_ngamma", autolib, "Autosecol", D5_SSH_methods)
        DRAGON5_AUTOLIB.compute_relative_differences_Rates("Gd157_ngamma", autolib, "Autosecol", D5_SSH_methods)
        for method in D5_SSH_methods:
            DRAGON5_AUTOLIB.plot_relative_differences_XS("Gd157_ngamma", autolib, "SHEM295", [method])
            DRAGON5_AUTOLIB.plot_relative_differences_Rates("Gd157_ngamma", autolib, "SHEM295", [method])
    # Post-treatment for SERPENT2 results
    # Make a case disjunction for XS and rates : due to amount of time/no necessity to obtain XS at all BU steps :
    # only BU=0 is considered for XS
    # All BU are considered for rates
    S2_case_name = "HOM_UOX_Gd157"
    SERPENT2_case = PT_S2(S2_case_name, mesh_objects, ["PyNjoy2016", "oldlib"], range(0,1), save_path)
    SERPENT2_case.parse_S2_outputs(path_to_S2_results)

    # Compare D5 results for different AUTOlibs wwith S2 results
    comparison_D5_S2 = CD5S2("HOM_UOX_Gd157_rates_XS_study", "AUTOLIB", DRAGON5_AUTOLIB, SERPENT2_case, 
                            S2_libs = ["oldlib", "PyNjoy2016"], compo_keywords = AUTOlibs, isotopes = ["Gd157", "U238"], self_shielding_methods = SSH_methods.values(),
                            save_path =   save_path_comparison)
    all_ssh_methods = ["Autosecol","RSE","PT","SUBG"]
    diff_data_rates = {"PyNjoy2016":{}, "oldlib":{}}
    diff_data_XS = {"PyNjoy2016":{}, "oldlib":{}}
    for autolib in AUTOlibs:
        comparison_D5_S2.plot_rates_D5_S2("Gd157_ngamma", autolib, "SHEM295", all_ssh_methods)
        comparison_D5_S2.plot_XS_D5_S2("Gd157_ngamma", autolib, "SHEM295", all_ssh_methods)
        # renormalize reaction rates to the same value for both D5 and S2
        comparison_D5_S2.renorm_rates("Gd157_ngamma", autolib, "SHEM295")
        comparison_D5_S2.compare_reaction_rates_and_XS("Gd157_ngamma", autolib, all_ssh_methods)

        # Plot XS and rates for all ssh methods on the same plot
        comparison_D5_S2.plot_diff_rates_D5_S2("Gd157_ngamma", autolib, "SHEM295", "PyNjoy2016", all_ssh_methods)
        comparison_D5_S2.plot_diff_rates_D5_S2("Gd157_ngamma", autolib, "SHEM295", "oldlib", all_ssh_methods)
        comparison_D5_S2.plot_diff_XS_D5_S2("Gd157_ngamma", autolib, "SHEM295", "PyNjoy2016", all_ssh_methods)
        comparison_D5_S2.plot_diff_XS_D5_S2("Gd157_ngamma", autolib, "SHEM295", "oldlib", all_ssh_methods)

        # plot on each individual plot : 1 per ssh method
        for ssh_method in all_ssh_methods:
            comparison_D5_S2.plot_diff_rates_D5_S2("Gd157_ngamma", autolib, "SHEM295", "PyNjoy2016", [ssh_method])
            comparison_D5_S2.plot_diff_rates_D5_S2("Gd157_ngamma", autolib, "SHEM295", "oldlib", [ssh_method])
            comparison_D5_S2.plot_diff_XS_D5_S2("Gd157_ngamma", autolib, "SHEM295", "PyNjoy2016", [ssh_method])
            comparison_D5_S2.plot_diff_XS_D5_S2("Gd157_ngamma", autolib, "SHEM295", "oldlib", [ssh_method])


        # find top differences
        diff_data_rates["PyNjoy2016"][autolib] = comparison_D5_S2.find_top_differences_rates("PyNjoy2016", "Gd157_ngamma", autolib, "SHEM295", 5)
        diff_data_rates["oldlib"][autolib] = comparison_D5_S2.find_top_differences_rates("oldlib", "Gd157_ngamma", autolib, "SHEM295", 5)
        diff_data_XS["PyNjoy2016"][autolib] = comparison_D5_S2.find_top_differences_XS("PyNjoy2016", "Gd157_ngamma", autolib, "SHEM295", 5)
        diff_data_XS["oldlib"][autolib] = comparison_D5_S2.find_top_differences_XS("oldlib", "Gd157_ngamma", autolib, "SHEM295", 5)
    # Analyse the results
    # What needs to be checked and analyzed : 
        # - What are the most important differences in the reaction rates and cross sections
        # - Are they consistent across the different energy meshes : ie where are they located on the energy spectrum
        # - Are differences on rates and XS fully correlated : if so can focus on XS only 
        # - Are differences on rates and XS consistent across the different self-shielding methods
        # - Are differences on rates and XS consistent across the different S2 libaries

    # Generate LaTeX tables for the top differences
    for autolib in AUTOlibs:
        for library in ["PyNjoy2016", "oldlib"]:
            diff_data_rates_table = generate_latex_table_autolib(library, autolib, diff_data_rates[library][autolib])
            diff_data_XS_table = generate_latex_table_autolib(library, autolib, diff_data_XS[library][autolib])

            with open(f"{save_path_comparison}/top_diff_rates_{library}_{autolib}.tex", "w") as f:
                f.write(diff_data_rates_table)
            
            with open(f"{save_path_comparison}/top_diff_XS_{library}_{autolib}.tex", "w") as f:
                f.write(diff_data_XS_table)

    print("Post-treatment of HOM_UOX_Gd157 Autolib study done")
    return


def post_treat_HOM_UOX_Gd157_IRSET():
    # Path to Serpent2 results
    path_to_S2_results = "/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/HOM_UOX_Gd157/XS_rates_study"
    # Import COMPO object
    path = os.getcwd()
    name_compo='newCOMPO_Gd157_IRSET'
    os.chdir("Gd157_rates_XS_proc")
    pyCOMPO=lcm.new('LCM_INP',name_compo,impx=0)
    os.chdir(path)
    # Creation of results directory
    path=os.getcwd()
    save_path = f"Gd157_Rates_and_XS_results_PyGan/IRSET_study"
    save_path_comparison = f"{save_path}/comparison_D5_S2"
    a=os.path.exists(save_path)
    if a==False:
        os.makedirs(save_path)

    a=os.path.exists(save_path_comparison)
    if a==False:
        os.makedirs(save_path_comparison)
    case_name = "HOM_UOX_Gd157"
    meshes = ["SHEM295"]
    IRSETs = ["noIRSET", "noIRSETF", "noIRSETFF", "IRSET_101", "IRSET_101F", "IRSET_101FF", "IRSETNone", "IRSETNoneFF"]
    SSH_methods = {"A":"Autosecol","R":"RSE","P":"PT","S":"SUBG","MC":"PTMC","SL":"PTSL"}
    D5_SSH_methods = ["RSE","PT","SUBG", "PTMC", "PTSL"]
    
    Gd157_ngamma = {}
    U238_ngamma =  {}
    fluxes = {}
    mesh_objects = {"SHEM295":None}

    for IRSET in IRSETs:
        Gd157_ngamma[IRSET] = {}
        U238_ngamma[IRSET] = {}
        fluxes[IRSET] = {}
        DIR = IRSET
        energyMESH = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][0]["ENERGY"]
        print(f"The {IRSET} mesh has {len(energyMESH)} energy bounds, so {len(energyMESH)-1} energy groups")
        MESH_obj = ENEMESH("SHEM295",energyMESH,1.0E+07,"eV")
        mesh_objects["SHEM295"] = MESH_obj
        SSH_methods_par = re.findall(r"[A-Za-z]+", pyCOMPO[DIR]["GLOBAL"]['pval00000001'].strip())
        for i in range(len(SSH_methods_par)):
            Gd157_ngamma[IRSET][SSH_methods[SSH_methods_par[i]]] = {"rates":[], "XS":[]}
            U238_ngamma[IRSET][SSH_methods[SSH_methods_par[i]]] = {"rates":[], "XS":[]}
            fluxes[IRSET][SSH_methods[SSH_methods_par[i]]] = []
            
            isotope = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][1]['ALIAS'][0:5]

            N_U238 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESDENS'][0]
            N_Gd157 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESDENS'][1]

            XS_NGAMMA_U238 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][0]['NG'] 
            XS_NGAMMA_Gd157 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][1]['NG'] 

            PHI = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][0]['NWT0']  

            # Reconstruct reaction rates
            NGAMMA_Gd157_rates = np.array(XS_NGAMMA_Gd157)*np.array(PHI)*N_Gd157
            NGAMMA_U238_rates = np.array(XS_NGAMMA_U238)*np.array(PHI)*N_U238

            # Store the results in the dictionaries

            Gd157_ngamma[IRSET][SSH_methods[SSH_methods_par[i]]]["rates"] = NGAMMA_Gd157_rates
            U238_ngamma[IRSET][SSH_methods[SSH_methods_par[i]]]["rates"] = NGAMMA_U238_rates

            Gd157_ngamma[IRSET][SSH_methods[SSH_methods_par[i]]]["XS"] = XS_NGAMMA_Gd157
            U238_ngamma[IRSET][SSH_methods[SSH_methods_par[i]]]["XS"] = XS_NGAMMA_U238

            fluxes[IRSET][SSH_methods[SSH_methods_par[i]]] = PHI
    
    # Call post treatment class
    DRAGON5_IRSET = PT_D5(case_name, "IRSET", mesh_objects, ["Gd157_ngamma", "U238_ngamma"], IRSETs, SSH_methods.values(), save_path)
    DRAGON5_IRSET.set_reaction_data("Gd157_ngamma",Gd157_ngamma)
    DRAGON5_IRSET.set_reaction_data("U238_ngamma",U238_ngamma)
    DRAGON5_IRSET.set_fluxes(fluxes)

    # Plot the results
    for irset in IRSETs:
        DRAGON5_IRSET.plot_fluxes(irset, "SHEM295")
        DRAGON5_IRSET.plot_XS(irset, "SHEM295", "Gd157_ngamma")
        DRAGON5_IRSET.plot_XS(irset, "SHEM295", "U238_ngamma")
        DRAGON5_IRSET.plot_reaction_rates(irset, "SHEM295", "Gd157_ngamma")
        DRAGON5_IRSET.plot_reaction_rates(irset, "SHEM295", "U238_ngamma")
        if irset == "noIRSET" or irset == "noIRSETF" or irset == "noIRSETFF":
            D5_SSH_methods = ["Autosecol","RSE","PT","SUBG", "PTMC","PTSL"] # Similarly as previously, IRSET not defined for TONE, SHI and SUBG
        else:
            D5_SSH_methods = ["Autosecol","RSE","PT", "PTMC", "PTSL"]
        DRAGON5_IRSET.compute_relative_differences_XS("Gd157_ngamma", irset, "Autosecol", D5_SSH_methods)
        DRAGON5_IRSET.compute_relative_differences_Rates("Gd157_ngamma", irset, "Autosecol", D5_SSH_methods)
        for method in D5_SSH_methods:
            DRAGON5_IRSET.plot_relative_differences_XS("Gd157_ngamma", irset, "SHEM295", [method])
            DRAGON5_IRSET.plot_relative_differences_Rates("Gd157_ngamma", irset, "SHEM295", [method])

    # Post-treatment for SERPENT2 results
    # Make a case disjunction for XS and rates : due to amount of time/no necessity to obtain XS at all BU steps :
    # only BU=0 is considered for XS
    # All BU are considered for rates
    S2_case_name = "HOM_UOX_Gd157"
    SERPENT2_case = PT_S2(S2_case_name, mesh_objects, ["PyNjoy2016", "oldlib"], range(0,1), save_path)
    SERPENT2_case.parse_S2_outputs(path_to_S2_results)

    # Compare D5 results for different IRSET settings with S2 results
    comparison_D5_S2 = CD5S2("HOM_UOX_Gd157_rates_XS_study", "IRSET", DRAGON5_IRSET, SERPENT2_case, 
                            S2_libs = ["oldlib", "PyNjoy2016"], compo_keywords = IRSETs, isotopes = ["Gd157", "U238"], self_shielding_methods = SSH_methods.values(),
                            save_path =   save_path_comparison)
    all_ssh_methods = ["Autosecol","RSE","PT","SUBG"]
    diff_data_rates = {"PyNjoy2016":{}, "oldlib":{}}
    diff_data_XS = {"PyNjoy2016":{}, "oldlib":{}}
    for irset in IRSETs:
        if irset == "noIRSET" or irset == "noIRSETF" or irset == "noIRSETFF":
            all_ssh_methods = ["Autosecol","RSE","PT","SUBG", "PTMC","PTSL"] # Similarly as previously, IRSET not defined for TONE, SHI and SUBG
        else:
            all_ssh_methods = ["Autosecol","RSE","PT", "PTMC", "PTSL"]
        comparison_D5_S2.plot_rates_D5_S2("Gd157_ngamma", irset, "SHEM295", all_ssh_methods)
        comparison_D5_S2.plot_XS_D5_S2("Gd157_ngamma", irset, "SHEM295", all_ssh_methods)
        # renormalize reaction rates to the same value for both D5 and S2
        comparison_D5_S2.renorm_rates("Gd157_ngamma", irset, "SHEM295")
        comparison_D5_S2.compare_reaction_rates_and_XS("Gd157_ngamma", irset, all_ssh_methods)

        # Plot XS and rates for all ssh methods on the same plot
        comparison_D5_S2.plot_diff_rates_D5_S2("Gd157_ngamma", irset, "SHEM295", "PyNjoy2016", all_ssh_methods)
        comparison_D5_S2.plot_diff_rates_D5_S2("Gd157_ngamma", irset, "SHEM295", "oldlib", all_ssh_methods)
        comparison_D5_S2.plot_diff_XS_D5_S2("Gd157_ngamma", irset, "SHEM295", "PyNjoy2016", all_ssh_methods)
        comparison_D5_S2.plot_diff_XS_D5_S2("Gd157_ngamma", irset, "SHEM295", "oldlib", all_ssh_methods)

        # plot on each individual plot : 1 per ssh method
        for ssh_method in all_ssh_methods:
            comparison_D5_S2.plot_diff_rates_D5_S2("Gd157_ngamma", irset, "SHEM295", "PyNjoy2016", [ssh_method])
            comparison_D5_S2.plot_diff_rates_D5_S2("Gd157_ngamma", irset, "SHEM295", "oldlib", [ssh_method])
            comparison_D5_S2.plot_diff_XS_D5_S2("Gd157_ngamma", irset, "SHEM295", "PyNjoy2016", [ssh_method])
            comparison_D5_S2.plot_diff_XS_D5_S2("Gd157_ngamma", irset, "SHEM295", "oldlib", [ssh_method])


        # find top differences
        diff_data_rates["PyNjoy2016"][irset] = comparison_D5_S2.find_top_differences_rates("PyNjoy2016", "Gd157_ngamma", irset, "SHEM295", 5)
        diff_data_rates["oldlib"][irset] = comparison_D5_S2.find_top_differences_rates("oldlib", "Gd157_ngamma", irset, "SHEM295", 5)
        diff_data_XS["PyNjoy2016"][irset] = comparison_D5_S2.find_top_differences_XS("PyNjoy2016", "Gd157_ngamma", irset, "SHEM295", 5)
        diff_data_XS["oldlib"][irset] = comparison_D5_S2.find_top_differences_XS("oldlib", "Gd157_ngamma", irset, "SHEM295", 5)
    # Analyse the results
    # What needs to be checked and analyzed : 
        # - What are the most important differences in the reaction rates and cross sections
        # - Are they consistent across the different energy meshes : ie where are they located on the energy spectrum
        # - Are differences on rates and XS fully correlated : if so can focus on XS only 
        # - Are differences on rates and XS consistent across the different self-shielding methods
        # - Are differences on rates and XS consistent across the different S2 libaries

    # Generate LaTeX tables for the top differences
    for irset in IRSETs:
        for library in ["PyNjoy2016", "oldlib"]:
            diff_data_rates_table = generate_latex_table_IRSET(library, irset, diff_data_rates[library][irset])
            diff_data_XS_table = generate_latex_table_IRSET(library, irset, diff_data_XS[library][irset])

            with open(f"{save_path_comparison}/top_diff_rates_{library}_{irset}.tex", "w") as f:
                f.write(diff_data_rates_table)
            
            with open(f"{save_path_comparison}/top_diff_XS_{library}_{irset}.tex", "w") as f:
                f.write(diff_data_XS_table)

    print("Post-treatment of HOM_UOX_Gd157 IRSET study done")
    return


def post_treat_HOM_UOX_Gd157_TONE_SHI_IRSET():
    # Path to Serpent2 results
    path_to_S2_results = "/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/HOM_UOX_Gd157/XS_rates_study"
    # Import COMPO object
    path = os.getcwd()
    name_compo='COMPO_Gd157_IR_TONE_SHI'
    os.chdir("Gd157_rates_XS_proc")
    pyCOMPO=lcm.new('LCM_INP',name_compo,impx=0)
    os.chdir(path)
    # Creation of results directory
    path=os.getcwd()
    save_path = f"Gd157_Rates_and_XS_results_PyGan/TONE_SHI_IRSET_study"
    save_path_comparison = f"{save_path}/comparison_D5_S2"
    a=os.path.exists(save_path)
    if a==False:
        os.makedirs(save_path)

    a=os.path.exists(save_path_comparison)
    if a==False:
        os.makedirs(save_path_comparison)
    case_name = "HOM_UOX_Gd157"
    meshes = ["SHEM295"]
    IRSETs = ["noIRSET", "IRSET_97", "IRSET_101", "IRSET_105", "IRSETNone"]
    SSH_methods = {"A":"Autosecol","R":"RSE","P":"PT","S":"SUBG","SHI": "SHI","T":"TONE"}
    D5_SSH_methods = ["RSE","PT","SUBG", "SHI","TONE"]
    Gd157_nGamma = {}
    U238_nGamma = {}
    fluxes = {}
    mesh_objects = {"SHEM295":None}
    for IRSET in IRSETs:
        Gd157_nGamma[IRSET] = {}
        U238_nGamma[IRSET] = {}
        fluxes[IRSET] = {}
        DIR = IRSET
        energyMESH = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][0]["ENERGY"]
        print(f"The {IRSET} mesh has {len(energyMESH)} energy bounds, so {len(energyMESH)-1} energy groups")
        MESH_obj = ENEMESH("SHEM295",energyMESH,1.0E+07,"eV")
        mesh_objects["SHEM295"] = MESH_obj
        SSH_methods_par = re.findall(r"[A-Za-z]+", pyCOMPO[DIR]["GLOBAL"]['pval00000001'].strip())
        for i in range(len(SSH_methods_par)):
            Gd157_nGamma[IRSET][SSH_methods[SSH_methods_par[i]]] = {"rates":[], "XS":[]}
            U238_nGamma[IRSET][SSH_methods[SSH_methods_par[i]]] = {"rates":[], "XS":[]}
            fluxes[IRSET][SSH_methods[SSH_methods_par[i]]] = []

            isotope = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][1]['ALIAS'][0:5]

            N_U238 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESDENS'][0]
            N_Gd157 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESDENS'][1]

            XS_NGAMMA_U238 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][0]['NG'] 
            XS_NGAMMA_Gd157 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][1]['NG'] 

            PHI = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][0]['NWT0']  

            # Reconstruct reaction
            NGAMMA_Gd157_rates = np.array(XS_NGAMMA_Gd157)*np.array(PHI)*N_Gd157
            NGAMMA_U238_rates = np.array(XS_NGAMMA_U238)*np.array(PHI)*N_U238

            # Store the results in the dictionaries

            Gd157_nGamma[IRSET][SSH_methods[SSH_methods_par[i]]]["rates"] = NGAMMA_Gd157_rates
            U238_nGamma[IRSET][SSH_methods[SSH_methods_par[i]]]["rates"] = NGAMMA_U238_rates

            Gd157_nGamma[IRSET][SSH_methods[SSH_methods_par[i]]]["XS"] = XS_NGAMMA_Gd157
            U238_nGamma[IRSET][SSH_methods[SSH_methods_par[i]]]["XS"] = XS_NGAMMA_U238

            fluxes[IRSET][SSH_methods[SSH_methods_par[i]]] = PHI

    # Call post treatment class
            
    DRAGON5_IRSET_extended = PT_D5(case_name, "IRSET_TONE_SHI", mesh_objects, ["Gd157_ngamma", "U238_ngamma"], IRSETs, SSH_methods.values(), save_path)
    DRAGON5_IRSET_extended.set_reaction_data("Gd157_ngamma",Gd157_nGamma)
    DRAGON5_IRSET_extended.set_reaction_data("U238_ngamma",U238_nGamma)
    DRAGON5_IRSET_extended.set_fluxes(fluxes)

    # Plot the results
    
    for irset in IRSETs:
        DRAGON5_IRSET_extended.plot_fluxes(irset, "SHEM295")
        DRAGON5_IRSET_extended.plot_XS(irset, "SHEM295", "Gd157_ngamma")
        DRAGON5_IRSET_extended.plot_XS(irset, "SHEM295", "U238_ngamma")
        DRAGON5_IRSET_extended.plot_reaction_rates(irset, "SHEM295", "Gd157_ngamma")
        DRAGON5_IRSET_extended.plot_reaction_rates(irset, "SHEM295", "U238_ngamma")
        if irset != "noIRSET": 
            D5_SSH_methods = ["RSE","PT"]
        else:
            D5_SSH_methods = ["RSE","PT","SUBG", "SHI","TONE"]   # IRSET not available for TONE, SHI and SUBG     
        DRAGON5_IRSET_extended.compute_relative_differences_XS("Gd157_ngamma", irset, "Autosecol", D5_SSH_methods)
        DRAGON5_IRSET_extended.compute_relative_differences_Rates("Gd157_ngamma", irset, "Autosecol", D5_SSH_methods)
        for method in D5_SSH_methods:
            DRAGON5_IRSET_extended.plot_relative_differences_XS("Gd157_ngamma", irset, "SHEM295", [method])
            DRAGON5_IRSET_extended.plot_relative_differences_Rates("Gd157_ngamma", irset, "SHEM295", [method])

    # Post-treatment for SERPENT2 results
    # Make a case disjunction for XS and rates : due to amount of time/no necessity to obtain XS at all BU steps :
    # only BU=0 is considered for XS
    # All BU are considered for rates
    S2_case_name = "HOM_UOX_Gd157"
    SERPENT2_case = PT_S2(S2_case_name, mesh_objects, ["PyNjoy2016", "oldlib"], range(0,1), save_path)
    SERPENT2_case.parse_S2_outputs(path_to_S2_results)

    # Compare D5 results for different IRSET settings with S2 results
    comparison_D5_S2 = CD5S2("HOM_UOX_Gd157_rates_XS_study", "IRSET", DRAGON5_IRSET_extended, SERPENT2_case, 
                            S2_libs = ["oldlib", "PyNjoy2016"], compo_keywords = IRSETs, isotopes = ["Gd157", "U238"], self_shielding_methods = SSH_methods.values(),
                            save_path =   save_path_comparison)
    all_ssh_methods = ["Autosecol","RSE","PT","SUBG", "SHI","TONE"]
    diff_data_rates = {"PyNjoy2016":{}, "oldlib":{}}
    diff_data_XS = {"PyNjoy2016":{}, "oldlib":{}}
    for irset in IRSETs:
        if irset != "noIRSET":
            all_ssh_methods = ["Autosecol","RSE","PT"]
        else:
            all_ssh_methods = ["Autosecol","RSE","PT","SUBG", "SHI","TONE"] # Similarly as previously, IRSET not defined for TONE, SHI and SUBG
        comparison_D5_S2.plot_rates_D5_S2("Gd157_ngamma", irset, "SHEM295", all_ssh_methods)
        comparison_D5_S2.plot_XS_D5_S2("Gd157_ngamma", irset, "SHEM295", all_ssh_methods)
        # renormalize reaction rates to the same value for both D5 and S2
        comparison_D5_S2.renorm_rates("Gd157_ngamma", irset, "SHEM295")
        comparison_D5_S2.compare_reaction_rates_and_XS("Gd157_ngamma", irset, all_ssh_methods)

        # Plot XS and rates for all ssh methods on the same plot
        comparison_D5_S2.plot_diff_rates_D5_S2("Gd157_ngamma", irset, "SHEM295", "PyNjoy2016", all_ssh_methods)
        comparison_D5_S2.plot_diff_rates_D5_S2("Gd157_ngamma", irset, "SHEM295", "oldlib", all_ssh_methods)
        comparison_D5_S2.plot_diff_XS_D5_S2("Gd157_ngamma", irset, "SHEM295", "PyNjoy2016", all_ssh_methods)
        comparison_D5_S2.plot_diff_XS_D5_S2("Gd157_ngamma", irset, "SHEM295", "oldlib", all_ssh_methods)

        # plot on each individual plot : 1 per ssh method
        for ssh_method in all_ssh_methods:
            comparison_D5_S2.plot_diff_rates_D5_S2("Gd157_ngamma", irset, "SHEM295", "PyNjoy2016", [ssh_method])
            comparison_D5_S2.plot_diff_rates_D5_S2("Gd157_ngamma", irset, "SHEM295", "oldlib", [ssh_method])
            comparison_D5_S2.plot_diff_XS_D5_S2("Gd157_ngamma", irset, "SHEM295", "PyNjoy2016", [ssh_method])
            comparison_D5_S2.plot_diff_XS_D5_S2("Gd157_ngamma", irset, "SHEM295", "oldlib", [ssh_method])


        # find top differences
        diff_data_rates["PyNjoy2016"][irset] = comparison_D5_S2.find_top_differences_rates("PyNjoy2016", "Gd157_ngamma", irset, "SHEM295", 5)
        diff_data_rates["oldlib"][irset] = comparison_D5_S2.find_top_differences_rates("oldlib", "Gd157_ngamma", irset, "SHEM295", 5)
        diff_data_XS["PyNjoy2016"][irset] = comparison_D5_S2.find_top_differences_XS("PyNjoy2016", "Gd157_ngamma", irset, "SHEM295", 5)
        diff_data_XS["oldlib"][irset] = comparison_D5_S2.find_top_differences_XS("oldlib", "Gd157_ngamma", irset, "SHEM295", 5)

        
    # Analyse the results
    # What needs to be checked and analyzed : 
        # - What are the most important differences in the reaction rates and cross sections
        # - Are they consistent across the different energy meshes : ie where are they located on the energy spectrum
        # - Are differences on rates and XS fully correlated : if so can focus on XS only 
        # - Are differences on rates and XS consistent across the different self-shielding methods
        # - Are differences on rates and XS consistent across the different S2 libaries

    # Generate LaTeX tables for the top differences
    for irset in IRSETs:
        for library in ["PyNjoy2016", "oldlib"]:
            diff_data_rates_table = generate_latex_table_IRSET(library, irset, diff_data_rates[library][irset])
            diff_data_XS_table = generate_latex_table_IRSET(library, irset, diff_data_XS[library][irset])

            with open(f"{save_path_comparison}/top_diff_rates_{library}_{irset}.tex", "w") as f:
                f.write(diff_data_rates_table)
            
            with open(f"{save_path_comparison}/top_diff_XS_{library}_{irset}.tex", "w") as f:
                f.write(diff_data_XS_table)

    # Plot in self-shielding region : for SHEM295, by default between groups 31 and 206
    # for now focus on SHEM295, Gd157_ngamma, IRSET = noIRSET. Compare Autosecol, SHI and Tone in a first plot
    # Plot the reaction rates
    comparison_D5_S2.plot_zoom_rates_and_errors("PyNjoy2016", "Gd157_ngamma", "noIRSET", "SHEM295", ["Autosecol", "SHI", "TONE"], grmin=52, grmax=206)
    comparison_D5_S2.plot_zoom_XS_and_errors("oldlib", "Gd157_ngamma", "noIRSET", "SHEM295", ["Autosecol", "SHI", "TONE"], grmin=52, grmax=206)
    # Plot the cross sections 
    comparison_D5_S2.plot_zoom_XS_and_errors("PyNjoy2016", "Gd157_ngamma", "noIRSET", "SHEM295", ["Autosecol", "SHI", "TONE"], grmin=52, grmax=206)
    comparison_D5_S2.plot_zoom_XS_and_errors("oldlib", "Gd157_ngamma", "noIRSET", "SHEM295", ["Autosecol", "SHI", "TONE"], grmin=52, grmax=206)
    
    # For AUTO: only, zoom into groups 70 to 206
    # Plot the reaction rates
    comparison_D5_S2.plot_zoom_rates_and_errors("PyNjoy2016", "Gd157_ngamma", "noIRSET", "SHEM295", ["Autosecol"], grmin=70, grmax=206)
    comparison_D5_S2.plot_zoom_rates_and_errors("oldlib", "Gd157_ngamma", "noIRSET", "SHEM295", ["Autosecol"], grmin=70, grmax=206)
    # Plot the cross sections
    comparison_D5_S2.plot_zoom_XS_and_errors("PyNjoy2016", "Gd157_ngamma", "noIRSET", "SHEM295", ["Autosecol"], grmin=70, grmax=206)
    comparison_D5_S2.plot_zoom_XS_and_errors("oldlib", "Gd157_ngamma", "noIRSET", "SHEM295", ["Autosecol"], grmin=70, grmax=206)

    # For TONE: only, zoom into groups 70 to 206
    # Plot the reaction rates
    comparison_D5_S2.plot_zoom_rates_and_errors("PyNjoy2016", "Gd157_ngamma", "noIRSET", "SHEM295", ["TONE"], grmin=70, grmax=206)
    comparison_D5_S2.plot_zoom_rates_and_errors("oldlib", "Gd157_ngamma", "noIRSET", "SHEM295", ["TONE"], grmin=70, grmax=206)
    # Plot the cross sections
    comparison_D5_S2.plot_zoom_XS_and_errors("PyNjoy2016", "Gd157_ngamma", "noIRSET", "SHEM295", ["TONE"], grmin=70, grmax=206)
    comparison_D5_S2.plot_zoom_XS_and_errors("oldlib", "Gd157_ngamma", "noIRSET", "SHEM295", ["TONE"], grmin=70, grmax=206)

    # For SHI: only, zoom into groups 70 to 206
    # Plot the reaction rates
    comparison_D5_S2.plot_zoom_rates_and_errors("PyNjoy2016", "Gd157_ngamma", "noIRSET", "SHEM295", ["SHI"], grmin=70, grmax=206)
    comparison_D5_S2.plot_zoom_rates_and_errors("oldlib", "Gd157_ngamma", "noIRSET", "SHEM295", ["SHI"], grmin=70, grmax=206)
    # Plot the cross sections
    comparison_D5_S2.plot_zoom_XS_and_errors("PyNjoy2016", "Gd157_ngamma", "noIRSET", "SHEM295", ["SHI"], grmin=70, grmax=206)
    comparison_D5_S2.plot_zoom_XS_and_errors("oldlib", "Gd157_ngamma", "noIRSET", "SHEM295", ["SHI"], grmin=70, grmax=206)



    # For AUTO: / RSE / PT, zoom into groups 70 to 206, noIRSET
    # Plot the reaction rates
    comparison_D5_S2.plot_zoom_rates_and_errors("PyNjoy2016", "Gd157_ngamma", "noIRSET", "SHEM295", ["Autosecol","RSE","PT"], grmin=70, grmax=206)
    comparison_D5_S2.plot_zoom_rates_and_errors("oldlib", "Gd157_ngamma", "noIRSET", "SHEM295", ["Autosecol","RSE","PT"], grmin=70, grmax=206)
    # Plot the cross sections
    comparison_D5_S2.plot_zoom_XS_and_errors("PyNjoy2016", "Gd157_ngamma", "noIRSET", "SHEM295", ["Autosecol","RSE","PT"], grmin=70, grmax=206)
    comparison_D5_S2.plot_zoom_XS_and_errors("oldlib", "Gd157_ngamma", "noIRSET", "SHEM295", ["Autosecol","RSE","PT"], grmin=70, grmax=206)

    # For AUTO: / RSE / PT, zoom into groups 70 to 206, IRSET_101
    # Plot the reaction rates
    comparison_D5_S2.plot_zoom_rates_and_errors("PyNjoy2016", "Gd157_ngamma", "IRSET_101", "SHEM295", ["Autosecol","RSE","PT"], grmin=70, grmax=206)
    comparison_D5_S2.plot_zoom_rates_and_errors("oldlib", "Gd157_ngamma", "IRSET_101", "SHEM295", ["Autosecol","RSE","PT"], grmin=70, grmax=206)
    # Plot the cross sections
    comparison_D5_S2.plot_zoom_XS_and_errors("PyNjoy2016", "Gd157_ngamma", "IRSET_101", "SHEM295", ["Autosecol","RSE","PT"], grmin=70, grmax=206)
    comparison_D5_S2.plot_zoom_XS_and_errors("oldlib", "Gd157_ngamma", "IRSET_101", "SHEM295", ["Autosecol","RSE","PT"], grmin=70, grmax=206)

    
    print("Post-treatment of HOM_UOX_Gd157 IRSET extended study done")
    return



    

if __name__ == "__main__":
    print("Post-treating HOM_UOX_Gd157 results")
    post_treating_MESHES_study = True
    post_treating_AUTOlib_study = True
    post_treating_IRSET_study = True
    post_treating_TONE_SHI_IRSET_study = True
    if post_treating_MESHES_study:
        print("Post-treating HOM_UOX_Gd157_MESHES")
        post_treat_HOM_UOX_Gd157_MESHES()
    if post_treating_AUTOlib_study:
        print("Post-treating HOM_UOX_Gd157_AUTOlib")
        post_treat_HOM_UOX_Gd157_AUTOlib()
    if post_treating_IRSET_study:
        print("Post-treating HOM_UOX_Gd157_IRSET")
        post_treat_HOM_UOX_Gd157_IRSET() # Check Autosecol for noIRSETF, should work. Re-running the calculation for noIRSETF/Autosecol
    if post_treating_TONE_SHI_IRSET_study:
        print("Post-treating HOM_UOX_Gd157_IRtest")
        post_treat_HOM_UOX_Gd157_TONE_SHI_IRSET()
