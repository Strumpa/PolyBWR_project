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
from plotters import plot_pinwise_errors_BWR_assembly, plot_spectrum_comparison, plot_U238_abs_rates
import pandas as pd

from starterDD.MaterialProperties.material_mixture import MaterialMixture, Composition
from starterDD.MaterialProperties.material_mixture import parse_all_compositions_from_yaml
from starterDD.GeometryAnalysis.tdt_parser import read_material_mixture_indices_from_tdt_file
from starterDD.DDModel.DragonModel import CartesianAssemblyModel, FuelPinModel
from starterDD.DDModel.helpers import associate_material_to_rod_ID
matplotlib.use('Agg')

        
def analyze_computational_schemes(assembly_id, computational_scheme, evaluation="endfb8r1"):
    """
    Analyze the impact of different computational schemes on flux and rates accuracy for a BWR assembly.
    Compares Dragon5 results with Serpent2 reference.
    """
    ssh_sol = "IC"
    ssh_method = "RSE"
    fission_isotopes = ["U235", "U238"]
    deltas = {}
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
    
    if assembly_id == "GE14_DOM":
        COMPO_name = f"CPO_{assembly_id}_{computational_scheme}_glow_Zr_Fe_Gd152"    
    elif assembly_id == "AT10":
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
    print(f"For computational scheme : {computational_scheme}, ssh solution door : {ssh_sol}")
    print(f"keff_D5 = {keff_D5}, keff_S2 = {keff_S2}")
    delta_keff = (keff_D5 - keff_S2)*1e5
    print(f"error (D5-S2) on keff = {delta_keff:.1f}")
    deltas[computational_scheme] = {}
    deltas[computational_scheme][f"{ssh_method}_{ssh_sol}"] = {}
    deltas[computational_scheme][f"{ssh_method}_{ssh_sol}"]["keff"] = delta_keff

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
    deltas[computational_scheme][f"{ssh_method}_{ssh_sol}"]["rms_therm"] = rms_therm
    deltas[computational_scheme][f"{ssh_method}_{ssh_sol}"]["avg_therm"] = avg_therm
    deltas[computational_scheme][f"{ssh_method}_{ssh_sol}"]["max_therm"] = max_therm
    deltas[computational_scheme][f"{ssh_method}_{ssh_sol}"]["rms_fast"] = rms_fast
    deltas[computational_scheme][f"{ssh_method}_{ssh_sol}"]["avg_fast"] = avg_fast
    deltas[computational_scheme][f"{ssh_method}_{ssh_sol}"]["max_fast"] = max_fast
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

    print(deltas)
    
    
    

if __name__ == "__main__":


    #analyze_main_flux_geometry_refinement()
    #analyze_computational_schemes("GE14_DOM", "1L_MOC")
    #analyze_computational_schemes("AT10", "1L_MOC")
    analyze_computational_schemes("AT10", "00_RSE_IC_1L_MOC_ANIS6")
    #analyze_computational_schemes("GE14_DOM", "2L_IC_MOC")
    analyze_computational_schemes("AT10", "2L_IC_MOC")
    