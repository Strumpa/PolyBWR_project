##############################################################################
# Post treatment of Dragon5 calculations for BWR assembly benchmarks.
# Cases generated with starterDD + glow
# Assemblies : 
#   - ATRIUM-10
#   - GE14 (DOM case)
##############################################################################
#
# ----------------------------------------------------------------------------
#                                  IMPORTS                                   |
# ----------------------------------------------------------------------------

import os, shutil, sys
import lifo
import lcm
import cle2000
from collections import defaultdict
import random as rd
import numpy as np
import matplotlib
import matplotlib.cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from typing import List, Union
import serpentTools as st
from Serpent2_parsing import  parse_S2_pin_mat_det
from Dragon5_cpo_parsing import  parse_COMPO
from plotters import plot_pinwise_errors_BWR_assembly, plot_spectrum_comparison, plot_U238_abs_rates, plot_fluxes_U238_rates_perturbation
import pandas as pd

from starterDD.MaterialProperties.material_mixture import MaterialMixture, Composition
from starterDD.MaterialProperties.material_mixture import parse_all_compositions_from_yaml
from starterDD.GeometryAnalysis.tdt_parser import read_material_mixture_indices_from_tdt_file
from starterDD.DDModel.DragonModel import CartesianAssemblyModel, FuelPinModel
from starterDD.DDModel.helpers import associate_material_to_rod_ID
matplotlib.use('Agg')

        
def analyze_computational_schemes(assembly_id, computational_scheme, evaluation="endfb8r1", compo_name=None):
    """
    Analyze the impact of different computational schemes on flux and rates accuracy for a BWR assembly.
    Compares Dragon5 results with Serpent2 reference.
    """
    ssh_sol = "IC"
    ssh_method = "RSE"
    fission_isotopes = ["U235", "U238"]
    deltas = {}
    U238_rates_295g = {}
    fluxes_295g = {}
    times = [0,0,0]
    # Create a starterDD assemmbly model to recover the geometry, material and compositions for the cases of interest.
    path_to_yaml_geometry = f"{os.environ['HOME']}/working_dir/PolyBWR_project/Version5_wc/PyGan/data/starterDD_postproc_proc/configs/{assembly_id}/GEOM.yaml"
    #"/home/apollon/p117902/working_dir/PolyBWR_project/Version5_wc/PyGan/data/starterDD_postproc_proc/configs/GE14_DOM/material_compositions.yaml"
    path_to_yaml_compositions = f"{os.environ['HOME']}/working_dir/PolyBWR_project/Version5_wc/PyGan/data/starterDD_postproc_proc/configs/{assembly_id}/material_compositions.yaml"
    tdt_file_name = f"{assembly_id}_SSH_IC_TISO.dat"
    path_to_tdt = os.environ['tdt_data']
    
    compositions = parse_all_compositions_from_yaml(path_to_yaml_compositions)
    ROD_to_material = associate_material_to_rod_ID(path_to_yaml_compositions,
                                               path_to_yaml_geometry)
    
    assembly_model = CartesianAssemblyModel(name=assembly_id,
                                    tdt_file=path_to_tdt + "/" + tdt_file_name,
                                    geometry_description_yaml=path_to_yaml_geometry)
    
    assembly_model.set_rod_ID_to_material_mapping(ROD_to_material)
    assembly_model.set_uniform_temperatures(fuel_temperature=900.0, gap_temperature=600.0, coolant_temperature=600.0,moderator_temperature=600.0, structural_temperature=600.0, )
    assembly_model.analyze_lattice_description(build_pins=True)
    assembly_model.set_material_compositions(compositions)
    assembly_model.number_fuel_material_mixtures_by_pin()
    
    print(f"Assembly model for {assembly_id} created. Number of fuel pins = {assembly_model.get_postprocessing_lattice_info()['n_fuel_pins']}")
    if compo_name is not None:
        print(f"Using provided compo name : {compo_name}")
        COMPO_name = compo_name
    else:
        if assembly_id == "GE14_DOM" or assembly_id == "GE14_DOM-C":
            COMPO_name = f"CPO_{assembly_id}_{computational_scheme}_glow_Zr_Fe_Gd152"    
        elif assembly_id == "AT10" or assembly_id == "AT10-C":
            COMPO_name = f"CPO_{assembly_id}_{computational_scheme}_glow_Zr_Fe"
   
    keff_D5, fiss_rates_D5, fiss_over_abs_rates_D5, D5_abs_U238, FLUX_groups_D5, energy_mesh = parse_COMPO(COMPO_name, assembly_model, fission_isotopes)
    keff_S2, fission_rates_S2, ngamma_rates_S2, fiss_over_abs_S2, fluxes_S2, S2_abs_U238 = parse_S2_pin_mat_det(assembly_model, assembly_id, evaluation, fission_isotopes)
    #keff_S2, fission_rates_S2, ngamma_rates_S2 = parse_S2_pin_mat_det(name_case=S2_name_case, XS_lib_S2="endfb8r1_pynjoy2012_kerma", fission_isotopes=fission_isotopes, ngamma_isotopes=["U238", "Gd155", "Gd157"], bu=0)

    print(fission_rates_S2)
    # Normalise to nCells 
    nCells = assembly_model.get_postprocessing_lattice_info()["n_fuel_pins"]
    fiss_rates_D5[0] =  fiss_rates_D5[0] * nCells / np.sum(fiss_rates_D5[0])
    fission_rates_S2[0] = fission_rates_S2[0] * nCells / np.sum(fission_rates_S2[0])
    fiss_rates_D5[1] =  fiss_rates_D5[1] * nCells / np.sum(fiss_rates_D5[1])
    fission_rates_S2[1] = fission_rates_S2[1] * nCells / np.sum(fission_rates_S2[1])

    delta_therm = fiss_rates_D5[0] - fission_rates_S2[0]
    delta_fast = fiss_rates_D5[1] - fission_rates_S2[1]
    
    delta_rel_therm = [(fiss_rates_D5[0][i] - fission_rates_S2[0][i])*100 / fission_rates_S2[0][i] if fission_rates_S2[0][i]!=0 else 0.0 for i in range(len(fission_rates_S2[0]))]
    delta_rel_fast = [(fiss_rates_D5[1][i] - fission_rates_S2[1][i])*100 / fission_rates_S2[1][i] if fission_rates_S2[1][i]!=0 else 0.0 for i in range(len(fission_rates_S2[1]))]

    #print(f"For geometry name : {geometry_type}")
    print(f"For compo name : {COMPO_name}")
    print(f"For computational scheme : {computational_scheme}, ssh solution door : {ssh_sol}")
    print(f"keff_D5 = {keff_D5}, keff_S2 = {keff_S2}")
    delta_keff = (keff_D5 - keff_S2)*1e5
    print(f"error (D5-S2) on keff = {delta_keff:.1f}")
    
    deltas[f"{ssh_method}_{ssh_sol}"] = {}
    deltas[f"{ssh_method}_{ssh_sol}"]["keffD5"] = keff_D5
    deltas[f"{ssh_method}_{ssh_sol}"]["keffS2"] = keff_S2
    deltas[f"{ssh_method}_{ssh_sol}"]["delta_keff"] = delta_keff

    # Compute RMS, MAX and AVG errors on fission rates
    rms_therm = np.sqrt(np.mean(np.array(delta_rel_therm)**2))
    rms_fast = np.sqrt(np.mean(np.array(delta_rel_fast)**2))
    avg_therm = np.mean(np.abs(np.array(delta_rel_therm)))
    avg_fast = np.mean(np.abs(np.array(delta_rel_fast)))
    max_therm = np.max(np.abs(np.array(delta_rel_therm)))
    max_fast = np.max(np.abs(np.array(delta_rel_fast)))
    print(f"RMS errors on fission rates (thermal, fast) = {rms_therm:.2f} %, {rms_fast:.2f} %")
    print(f"AVG errors on fission rates (thermal, fast) = {avg_therm:.2f} %, {avg_fast:.2f} %")
    print(f"MAX errors on fission rates (thermal, fast) = {max_therm:.2f} %, {max_fast:.2f} %")
    deltas[f"{ssh_method}_{ssh_sol}"]["rms_therm"] = rms_therm
    deltas[f"{ssh_method}_{ssh_sol}"]["avg_therm"] = avg_therm
    deltas[f"{ssh_method}_{ssh_sol}"]["max_therm"] = max_therm
    deltas[f"{ssh_method}_{ssh_sol}"]["rms_fast"] = rms_fast
    deltas[f"{ssh_method}_{ssh_sol}"]["avg_fast"] = avg_fast
    deltas[f"{ssh_method}_{ssh_sol}"]["max_fast"] = max_fast
    total_time = times[0]
    uss_time = times[1]
    asm_flu_time = times[2]
    print(total_time)
    num_hours = int(float(total_time)/3600)
    mins = (float(total_time)/3600 - float(num_hours)) * 60 
    print(f"{delta_keff:.0f} & {rms_fast:.2f}\\% & {rms_therm:.2f}\\% & {max_fast:.2f}\\% &  {max_therm:.2f}\\% & $\\approx {num_hours}h{mins:.00f}$")
    print(f"USS time = {uss_time} s, ASM+FLU time = {asm_flu_time} s")
    # errors_rates, assembly_model, assembly_id, name_compo, calculation_opt, fig_name, evaluation,
    plot_pinwise_errors_BWR_assembly(np.array([delta_rel_therm,delta_rel_fast]), assembly_model, assembly_id, COMPO_name, calculation_opt = f"{ssh_method}_{ssh_sol}_{computational_scheme}", fig_name="fission_rates_diff", evaluation=evaluation)
    
    plot_spectrum_comparison(energy_mesh, FLUX_groups_D5, fluxes_S2, assembly_id, COMPO_name, calculation_opt = f"{ssh_method}_{ssh_sol}_{computational_scheme}", evaluation=evaluation)
    # energy_mesh, D5_U238_abs_rates, S2_U238_abs_rates, name_case, CPO_name, calculation_opt, evaluation
    plot_U238_abs_rates(energy_mesh, D5_abs_U238, S2_abs_U238, assembly_id, COMPO_name, calculation_opt = f"{ssh_method}_{ssh_sol}_{computational_scheme}", evaluation=evaluation)

    delta_fiss_over_abs_therm = fiss_over_abs_rates_D5[0] - fiss_over_abs_S2[0]
    print("thermal fiss/abs rates comparison :")
    print(f"fiss_over_abs_rates_D5[0] = {fiss_over_abs_rates_D5[0]}")
    print(f"fiss_over_abs_S2[0] = {fiss_over_abs_S2[0]}")
    print(f"delta_fiss_over_abs_therm = {delta_fiss_over_abs_therm}")
    delta_fiss_over_abs_fast = fiss_over_abs_rates_D5[1] - fiss_over_abs_S2[1]
    print("fast fiss/abs rates comparison :")
    print(f"fiss_over_abs_rates_D5[1] = {fiss_over_abs_rates_D5[1]}")
    print(f"fiss_over_abs_S2[1] = {fiss_over_abs_S2[1]}")
    print(f"delta_fiss_over_abs_fast = {delta_fiss_over_abs_fast}")
    delta_rel_fiss_over_abs_therm = [(fiss_over_abs_rates_D5[0][i] - fiss_over_abs_S2[0][i])*100 / fiss_over_abs_S2[0][i] if fiss_over_abs_S2[0][i]!=0 else 0.0 for i in range(len(fiss_over_abs_S2[0]))]
    delta_rel_fiss_over_abs_fast = [(fiss_over_abs_rates_D5[1][i] - fiss_over_abs_S2[1][i])*100 / fiss_over_abs_S2[1][i] if fiss_over_abs_S2[1][i]!=0 else 0.0 for i in range(len(fiss_over_abs_S2[1]))]
    plot_pinwise_errors_BWR_assembly(np.array([delta_fiss_over_abs_therm,delta_fiss_over_abs_fast]), assembly_model, assembly_id, COMPO_name, calculation_opt = f"{ssh_method}_{ssh_sol}_{computational_scheme}", fig_name="fiss_over_abs_diff", evaluation=evaluation)
    
    U238_rates_295g= {"D5": D5_abs_U238, "S2": S2_abs_U238, "energy_mesh": energy_mesh} 
    fluxes_295g = {"D5": FLUX_groups_D5, "S2": fluxes_S2, "energy_mesh": energy_mesh}
    
    return deltas, U238_rates_295g, fluxes_295g
    
def format_delta_in_latex_table(deltas, assembly_ids):
    for assembly_id in assembly_ids:
        print(f"\\multirow{{2}}{{*}}{{{assembly_id}}} & 1L_MOC & {deltas[assembly_id]['1L_MOC']['RSE_IC']['delta_keff']:.0f} & {deltas[assembly_id]['1L_MOC']['RSE_IC']['rms_fast']:.2f}\\% & {deltas[assembly_id]['1L_MOC']['RSE_IC']['rms_therm']:.2f}\\% & {deltas[assembly_id]['1L_MOC']['RSE_IC']['max_fast']:.2f}\\% &  {deltas[assembly_id]['1L_MOC']['RSE_IC']['max_therm']:.2f}\\% \\\\")
        print(f"& 2L_IC_MOC & {deltas[assembly_id]['2L_IC_MOC']['RSE_IC']['delta_keff']:.0f} & {deltas[assembly_id]['2L_IC_MOC']['RSE_IC']['rms_fast']:.2f}\\% & {deltas[assembly_id]['2L_IC_MOC']['RSE_IC']['rms_therm']:.2f}\\% & {deltas[assembly_id]['2L_IC_MOC']['RSE_IC']['max_fast']:.2f}\\% &  {deltas[assembly_id]['2L_IC_MOC']['RSE_IC']['max_therm']:.2f}\\% \\\\")

def compute_cross_worth(deltas, assembly_id, assembly_id_controlled):
    delta_controlled = deltas[assembly_id_controlled]
    delta_uncontrolled = deltas[assembly_id]
    ssh_opt = "RSE_IC"
    print(f"Computing cross worth for assembly {assembly_id} with reference to controlled case {assembly_id_controlled} for ssh solution {ssh_opt}")
    for scheme in delta_controlled.keys():
        rho_cross_D5 = (1/delta_uncontrolled[scheme][ssh_opt]["keffD5"] - 1/delta_controlled[scheme][ssh_opt]["keffD5"])*1e5
        rho_cross_S2 = (1/delta_uncontrolled[scheme][ssh_opt]["keffS2"] - 1/delta_controlled[scheme][ssh_opt]["keffS2"])*1e5
        print(f"Cross worth for scheme {scheme} : rho_cross_D5 = {rho_cross_D5:.1f} pcm, rho_cross_S2 = {rho_cross_S2:.1f} pcm")
        print(f"Cross worth difference (D5-S2) for scheme {scheme} : delta_rho_cross = {(rho_cross_D5 - rho_cross_S2):.1f} pcm")
        
            
    

if __name__ == "__main__":

    GE14_assembly_ids = []#["GE14_DOM", "GE14_DOM-C"]
    AT10_assembly_ids = []#["AT10", "AT10-C"]
    calculation_schemes = []#["1L_MOC", "2L_IC_MOC"]
    deltas = {}
    rates_and_spectra = {}
    #analyze_main_flux_geometry_refinement()
    #analyze_computational_schemes("GE14_DOM", "1L_MOC")
    #analyze_computational_schemes("AT10", "1L_MOC")
    #analyze_computational_schemes("AT10", "00_RSE_IC_1L_MOC_ANIS6")
    for assembly_id in GE14_assembly_ids:
        for scheme in calculation_schemes:
            delta_scheme, U238_rates, spectra = analyze_computational_schemes(assembly_id, scheme)
            if assembly_id not in deltas.keys():
                deltas[assembly_id] = {}
            if assembly_id not in rates_and_spectra.keys():
                rates_and_spectra[assembly_id] = {}
            deltas[assembly_id][scheme] = delta_scheme
            rates_and_spectra[assembly_id][scheme] = {"U238_rates": U238_rates, "spectra": spectra}
    for assembly_id in AT10_assembly_ids:
        for scheme in calculation_schemes:
            if assembly_id not in deltas.keys():
                deltas[assembly_id] = {}
            if assembly_id not in rates_and_spectra.keys():
                rates_and_spectra[assembly_id] = {}
            delta_scheme, U238_rates, spectra = analyze_computational_schemes(assembly_id, scheme)
            deltas[assembly_id][scheme] = delta_scheme
            rates_and_spectra[assembly_id][scheme] = {"U238_rates": U238_rates, "spectra": spectra}
            
            
            
    #df = pd.DataFrame.from_dict({(i,j): deltas[i][j] for i in deltas.keys() for j in deltas[i].keys()},orient='index')
    #df.to_csv("postproc_computational_schemes_comparison.csv")
    
    #print(deltas)
    #format_delta_in_latex_table(deltas, GE14_assembly_ids)
    #format_delta_in_latex_table(deltas, AT10_assembly_ids)
    #compute_cross_worth(deltas, "GE14_DOM", "GE14_DOM-C")
    #compute_cross_worth(deltas, "AT10", "AT10-C")
    
    #plot_fluxes_U238_rates_perturbation(rates_and_spectra, GE14_assembly_ids, scheme="1L_MOC", evaluation="endfb8r1")
    #plot_fluxes_U238_rates_perturbation(rates_and_spectra, AT10_assembly_ids, scheme="1L_MOC", evaluation="endfb8r1")
    
    ## Test taget compo with specific self-shielding options

    #print(deltas)
    # GE14 DOM cases
    # 1L MOC
    #deltas, U238_rates_295g, fluxes_295g = analyze_computational_schemes("GE14_DOM", "1L_MOC", evaluation="endfb8r1", compo_name="CPO_GE14_DOM_IC_1L_MOC_ANIS4_finer_water_rods_8_moder_3016_corners4")
    deltas, U238_rates_295g, fluxes_295g = analyze_computational_schemes("GE14_DOM", "1L_MOC", evaluation="endfb8r1", compo_name="CPO_GE14_DOM_IC_1L_MOC_ANIS4_finer_water_rods_8_moder_308_corners4")
    
    #deltas, U238_rates_295g, fluxes_295g = analyze_computational_schemes("GE14_DOM", "1L_MOC", evaluation="endfb8r1", compo_name="CPO_GE14_DOM_IC_1L_MOC_ANIS6_finer_water_rods")
    #deltas, U238_rates_295g, fluxes_295g = analyze_computational_schemes("GE14_DOM", "1L_MOC", evaluation="endfb8r1", compo_name="CPO_GE14_DOM_IC_1L_MOC_finer_water_rods_8_moder_3010")
    #deltas, U238_rates_295g, fluxes_295g = analyze_computational_schemes("GE14_DOM", "1L_MOC", evaluation="endfb8r1", compo_name="CPO_GE14_DOM_IC_1L_MOC_glow_Zr_Fe_Gd152_Cr_Hf_finer_water_rods")
    # 2L IC+MOC
    #deltas, U238_rates_295g, fluxes_295g = analyze_computational_schemes("GE14_DOM", "2L_IC_MOC", evaluation="endfb8r1", compo_name="CPO_GE14_DOM_IC_2L_IC_MOC_finer_water_rods_8_moder_3016_corners4")
    deltas, U238_rates_295g, fluxes_295g = analyze_computational_schemes("GE14_DOM", "2L_IC_MOC", evaluation="endfb8r1", compo_name="CPO_GE14_DOM_IC_2L_IC_MOC_finer_water_rods_8_moder_308_corners4")
    #deltas, U238_rates_295g, fluxes_295g = analyze_computational_schemes("GE14_DOM", "2L_IC_MOC", evaluation="endfb8r1", compo_name="CPO_GE14_DOM_IC_2L_IC_MOC_finer_water_rods_8_moder_3016")
    #deltas, U238_rates_295g, fluxes_295g = analyze_computational_schemes("GE14_DOM", "2L_IC_MOC", evaluation="endfb8r1", compo_name="CPO_GE14_DOM_2L_IC_MOC_glow_Zr_Fe_Gd152_Cr_Hf_finer_water_rods")
    #deltas, U238_rates_295g, fluxes_295g = analyze_computational_schemes("GE14_DOM", "2L_IC_MOC", evaluation="endfb8r1", compo_name="CPO_GE14_DOM_2L_IC_MOC_glow_Zr_Fe_Gd152_Cr_Hf_finer_water_rods")
    #deltas, U238_rates_295g, fluxes_295g = analyze_computational_schemes("GE14_DOM", "2L_IC_MOC", evaluation="endfb8r1", compo_name="CPO_GE14_DOM_IC_2L_MOC_finer_water_rods_8_moder_3010")
    #deltas, U238_rates_295g, fluxes_295g = analyze_computational_schemes("GE14_DOM", "2L_IC_MOC", evaluation="endfb8r1", compo_name="CPO_GE14_DOM_2L_IC_MOC_glow_Zr_Fe_Gd152_Cr_Hf")
    
    # GE14 DOM-C cases
    # 1L MOC
    #deltas, U238_rates_295g, fluxes_295g = analyze_computational_schemes("GE14_DOM-C", "1L_MOC", evaluation="endfb8r1", compo_name="CPO_GE14_DOM-C_1L_MOC_water_rods_cross_submesh_extend_and_split_tubes")
    deltas, U238_rates_295g, fluxes_295g = analyze_computational_schemes("GE14_DOM-C", "1L_MOC", evaluation="endfb8r1", compo_name="CPO_GE14_DOM-C_1L_MOC_water_rods_cross_submesh_extend")
    
    # 2L IC+MOC
    deltas, U238_rates_295g, fluxes_295g = analyze_computational_schemes("GE14_DOM-C", "2L_IC_MOC", evaluation="endfb8r1", compo_name="CPO_GE14_DOM-C_2L_IC_MOC_water_rods_cross_submesh_extend")
    #deltas, U238_rates_295g, fluxes_295g = analyze_computational_schemes("GE14_DOM-C", "2L_IC_MOC", evaluation="endfb8r1", compo_name="CPO_GE14_DOM-C_2L_IC_MOC_glow_Zr_Fe_Gd152_Cr_Hf")
    #deltas, U238_rates_295g, fluxes_295g = analyze_computational_schemes("GE14_DOM-C", "2L_IC_MOC", evaluation="endfb8r1", compo_name="CPO_GE14_DOM-C_2L_IC_MOC_water_rods_cross_submesh_extend_and_split_tube")