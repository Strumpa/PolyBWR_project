### PyGan script optimize flux calculations for the ATRIUM-10 BWR fuel assembly
# Author : R. Guasch
# Date : 2025-06-02, updated 2025-06-04
## 2025-06-04 : implmenting mix numbering per region

# Import modules from python / python-CLE2000 API / python-LCM API
import os
import shutil
import numpy as np
import pandas as pd
import lcm 
import cle2000
import lifo
import time

# Import custom modules for the ATRIUM-10 BWR fuel assembly
from create_geo import createGeo
from tracking_operator import trackFluxGeomSALT
from tracking_operator import trackSSHGeomSALT
from mix_handling import createLib
from self_shielding import selfShieldingUSS
from flux_calculation import fluxCalculationMOC
from edi_compo import ediCompo

save_dir = "ATRIUM10_results"
if not os.path.exists(save_dir):
    os.makedirs(save_dir)


########################################################### PARAMETER SELECTION ##########################################################################
exec = True

# Options from DRAGON calculation setup.
# Geometry parameters : ATRIUM-10 BWR fuel assembly
refinement_opt_name = "finest_geom" # "default", "fine1", "fine2", "fine3", "coolant_ring", "finest_on_Gd", "finest_on_Gd_coolant_ring"

refinement_options = {} 
# posible keys to "moderator" entry
# "default" -> coarse refinement of water in and out of assembly, windmills around cells
#"fine1" -> finer refinement in inner box and outer water. Ibn thermal group : still ~-4% error on fission rates in assembly corners, 1.4% RMS and 1% average,
#"fine2" -> splity in side water increased, same as fine1 apart from that
#"fine3" -> refining moderator in central box : to test
refinement_options["moderator"] = "fine2"
# Possible keys to "UOX_cells" entry :
# SECT_4_{n}, where n is the number of non discretized rings in cell
# coolant_ring_SECT_4_{n}, where n is the number of non discretized rings in cell : added coolant ring 
refinement_options["UOX_cells"] = "SECT_4_6" #"coolant_ring_SECT_4_6"
# Possible keys to "Gd_cells" entry :
# SECT_4_{n}, where n is the number of non discretized rings in cell
# coolant_ring_SECT_4_{n}, where n is the number of non discretized rings in cell : added coolant ring 
refinement_options["Gd_cells"] = "SECT_4_0" #"coolant_ring_SECT_4_0"


if refinement_opt_name == "finest_geom": 
    refinement_options["moderator"] = "fine2" 
    refinement_options["UOX_cells"] = "coolant_ring_SECT_4_0"
    refinement_options["Gd_cells"] = "coolant_ring_SECT_4_0"
    num_angles = 24
    line_density = 150.0
    batch = 3000
    
elif refinement_opt_name == "fine2" or refinement_opt_name == "fine1" or refinement_opt_name == "default":
    refinement_options["moderator"] = refinement_opt_name
    refinement_options["UOX_cells"] = "SECT_4_6"
    refinement_options["Gd_cells"] = "SECT_4_8"
    num_angles = 24
    line_density = 75.0
    batch = 750
    
elif refinement_opt_name == "coolant_ring":
    refinement_options["moderator"] = "fine2" 
    refinement_options["UOX_cells"] = "coolant_ring_SECT_4_6"
    refinement_options["Gd_cells"] = "coolant_ring_SECT_4_8"
    num_angles = 24
    line_density = 75.0
    batch = 750
    
elif refinement_opt_name == "finest_on_Gd_coolant_ring": 
    refinement_options["moderator"] = "fine2" 
    refinement_options["UOX_cells"] = "coolant_ring_SECT_4_6"
    refinement_options["Gd_cells"] = "coolant_ring_SECT_4_0"
    num_angles = 24
    line_density = 140.0
    batch = 3000
elif refinement_opt_name == "finest_on_Gd": 
    refinement_options["moderator"] = "fine2" 
    refinement_options["UOX_cells"] = "SECT_4_6"
    refinement_options["Gd_cells"] = "SECT_4_0"
    num_angles = 24
    line_density = 140.0
    batch = 1500
else:
    print(f"This combination of options {refinement_options} is not supported yet, please review options are add specific save name to support results post-treatment")

if refinement_options["moderator"] == "default":
    split_water_in_moderator_box = 10
    split_moderator_box = 2
    split_water_around_moderator_box = 2
    split_intra_assembly_coolant = 4
    split_assembly_box = 2
    split_out_assembly_moderator = [5,5]
elif refinement_options["moderator"] == "fine1":
    split_water_in_moderator_box = 20
    split_moderator_box = 2
    split_water_around_moderator_box = 2
    split_intra_assembly_coolant = 4
    split_assembly_box = 2
    split_out_assembly_moderator = [10,10]
elif refinement_options["moderator"] == "fine2":
    split_water_in_moderator_box = 20
    split_moderator_box = 2
    split_water_around_moderator_box = 2
    split_intra_assembly_coolant = 4
    split_assembly_box = 2
    split_out_assembly_moderator = [10,30]
elif refinement_options["moderator"] == "fine3":
    split_water_in_moderator_box = 30
    split_moderator_box = 2
    split_water_around_moderator_box = 2
    split_intra_assembly_coolant = 4
    split_assembly_box = 2
    split_out_assembly_moderator = [10,30]

mix_numbering_option =  "number_mix_families_per_region" # "number_mix_families_per_region" , "number_mix_families_per_enrichment"
name_geom = "AT10_ASSBLY"

# Tracking parameters : main flux geometry
# (24, 75.0) --> SALTLC: Global RMS, maximum and average errors (%) on region volumes :  0.15415    0.89377    0.00094, below 1% on max : satisfactory.
# For SECT 4 0 with coolant ring : (24, 150)
# SALTLC: Global RMS, maximum and average errors (%) on region volumes :     0.13502     0.95490    -0.00102, below 1% on max : satisfactory.
# Could be optimized further, but not necessary for the moment.
#num_angles = 24
#line_density = 150.0
#batch = 1500 # 750 was found to be "optimal" with 20 omp procs and (24, 75.0) tracking parameters for MOC. 
reflection_type = "TSPC"
anisotropy_level = 2 # Level of anisotropy for the tracking, can be 1 (isotropic), 2 (linearly anisotropic), 3 (anisotropy order P_2), or 4 (anisotropy order P_3). 
solution_door = "MOC"  # Flag to indicate whether the tracking should be modified for a MOC solution, or not., else : set it to CP (Collision Probability) tracking.
moc_angular_quandrature = "GAUS"
nmu = 4  # Number of polar angles for MOC tracking : conservation ensured up to the order of P_{nmu-1} scattering : # nmu = 4 -> P3 scattering
postscript_file = f"AT10_FIG_MAIN_{refinement_opt_name}.ps"

# Tracking parameters : self-shielding geometry
num_angles_ssh = 8
line_density_ssh = 25.0
solution_door_ssh = "PIJ" # "PIJ"
if solution_door_ssh == "PIJ":
    reflection_type_ssh = "TSPC"
elif solution_door_ssh == "IC":
    reflection_type_ssh = "TISO"
else: 
    print(f"Selected self-shielding solution door is not supported")
batch_ssh = 200
postscript_file_ssh = "AT10_FIG_SSH.ps"

# Parameters for the LIBRARY creation
draglib_name = "J311_295" # "endfb8r1_295" # "J311_295"
self_shielding_method = "RSE"  # Method to be used for self-shielding calculations, "PT" for Mathematical Probability Tables, "SUBG" for Physical Probaility tables, "RSE" for Resonant Spectrum Expansion.
resonance_correlation = "NOCORR"  # Specify if the resonance correlation model should be applied. Only available for "RSE" and "PT". This will use a correlation model to treat reonances of U238, Pu240 and Gd157.
transport_correction = "APOL"
composition_option = "AT10_void_0"  # Specify which composition of mixes should be used for the LIBRARY creation. For now "AT10_void_0" and "AT10_void_40" are available.

# USS: call parameters :
ssh_option = "default"  # Option to specify specific groupings of self-shielding regions, to be tested. "default" is default from USS: based on LIB: data, "RSE" is for Resonant Spectrum Expansion method, and groups all U8, U5 and Zr isotopes in a single self-shielding region --> to be tested.

# EDI: and COMPO: calls to save the results to a MULTICOMPO file
if transport_correction == "NONE":
    name_compo = f"_CPO_n{num_angles}_ld{int(line_density)}_n{num_angles_ssh}_ld{int(line_density_ssh)}_{reflection_type}_{anisotropy_level}_{solution_door_ssh}_{solution_door}_{moc_angular_quandrature}_{nmu}_{batch}_{batch_ssh}"
else:
    name_compo = f"_CPO_n{num_angles}_ld{int(line_density)}_n{num_angles_ssh}_ld{int(line_density_ssh)}_{reflection_type}_{anisotropy_level}_{solution_door_ssh}_{solution_door}_{moc_angular_quandrature}_{nmu}_{batch}_{batch_ssh}_{transport_correction}"
########################################################### END OF PARAMETER SELECTION #####################################################################

########################################################### GEO: : geometry creation ########################################################################
start_time = time.time()
geo_flx, geo_ssh, connectivity_dict = createGeo(name_geom, split_water_in_moderator_box, split_moderator_box, split_water_around_moderator_box, 
                                                split_intra_assembly_coolant, split_assembly_box, split_out_assembly_moderator,
                                                mix_numbering_option, refinement_options,
                                                reflection_type_ssh, reflection_type)
time_creating_geo = time.time() - start_time
current_time = time.time()


######################################################## SALT: : main flux geometry tracking ################################################################

# Track the geometry using the SALT: module
track_lcm, track_binary, figure = trackFluxGeomSALT(geo_flx, num_angles, line_density, reflection_type, anisotropy_level, solution_door, moc_angular_quandrature, nmu, batch, postscript_file)
time_tracking_flux = time.time() - current_time
current_time = time.time()

######################################################## SALT: : self-shielding geometry tracking ###########################################################

# Track the self-shielding geometry
track_lcm_ssh, track_binary_ssh, figure_ssh = trackSSHGeomSALT(geo_ssh, num_angles_ssh, line_density_ssh, reflection_type_ssh, batch_ssh, postscript_file_ssh, solution_door_ssh)
time_tracking_ssh = time.time() - current_time
current_time = time.time()

if exec:
    ###################################################### LIB: : LIBRARY creation ############################################################
    # Create the LIBRARY according to the options selected.
    lib_lcm = createLib(mix_numbering_option, draglib_name, anisotropy_level, self_shielding_method, resonance_correlation, transport_correction, composition_option, connectivity_dict)
    time_creating_lib = time.time() - current_time
    current_time = time.time()


    ###################################################### USS: : self-shielding calculations ##########################################
    # Perform the self-shielding calculations
    lib_ssh = selfShieldingUSS(mix_numbering_option, lib_lcm, track_lcm_ssh, track_binary_ssh, name_geom, ssh_option, connectivity_dict)
    # Time taken for self-shielding calculations
    time_self_shielding = time.time() - current_time
    current_time = time.time()


    ###################################################### ASM: and FLU: calls ######################################################

    # Main flux calculation 
    if solution_door == "PIJ":
        keff, flux_lcm = fluxCalculationPIJ(track_lcm, track_binary, lib_ssh)
    elif solution_door == "MOC":
        keff, flux_lcm = fluxCalculationMOC(track_lcm, track_binary, lib_ssh)
    elif solution_door == "IC":
        keff, flux_lcm = fluxCalculationIC(track_lcm, track_binary, lib_ssh) # to be updated for IC solution door
    else:
        print(f"Selected solution door {solution_door} is not supported for flux calculation")
    print(f"Effective multiplication factor from flux calculation : {keff}")
    # Time taken for flux calculation
    time_flux_calculation = time.time() - current_time
    current_time = time.time()

    #################################################### EDI: and COMPO: ############################################################
    # Call EDI: and COMPO: to save the results to a MULTICOMPO file
    compo_lcm = ediCompo(mix_numbering_option, flux_lcm, lib_ssh, track_lcm, name_compo, save_option="SAVE", mix_connectivity_dict=connectivity_dict)
    # Time taken for EDI: and COMPO: calls
    time_edi_compo = time.time() - current_time
    current_time = time.time()


    ######################################################### EXPORTS ############################################################
    if mix_numbering_option == "number_mix_families_per_region":
        numbering_save_opt = "region_num"
    elif mix_numbering_option == "number_mix_families_per_enrichment":
        numbering_save_opt = "enrich_num"
    # Save the MULTICOMPO to a specific directory
    save_dir_case = f"{save_dir}/{refinement_opt_name}_{composition_option}_{draglib_name}_{self_shielding_method}_{resonance_correlation}_{numbering_save_opt}"
    if not os.path.exists(save_dir_case):
        os.makedirs(save_dir_case)
    # Save the LCM objects to the specified directory
    shutil.copyfile(name_compo, f"{save_dir_case}/{name_compo}")


    # Save times per module in the same save_dir_case directory, in a pd.DataFrame
    times_dict = {
        "GEO: call (s)": time_creating_geo,
        "SALT: call [flux] (s)": time_tracking_flux,
        "SALT: call [ssh] (s)": time_tracking_ssh,
        f"LIB: call ({self_shielding_method}) (s)": time_creating_lib,
        "USS: call (s)": time_self_shielding,
        "ASM: + FLU: calls (s)": time_flux_calculation,
        "EDI: + COMPO: calls (s)": time_edi_compo
    }

    times_df = pd.DataFrame.from_dict(times_dict, orient='index', columns=['Time (s)'])
    # Save the DataFrame to a CSV file in the specified directory
    times_df.to_csv(f"{save_dir_case}/df_module_times_n{num_angles}_ld{int(line_density)}_n2{num_angles_ssh}_ld2{int(line_density_ssh)}_{reflection_type}_{anisotropy_level}_{solution_door}_{moc_angular_quandrature}_{nmu}_{batch}_{batch_ssh}.csv")