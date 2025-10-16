### PyGan script to perform burnup calculations for an ATRIUM-10 BWR fuel assembly
# Author : R. Guasch
# Date : modified from ATRIUM_10_opti on 18/09/2025

# Import modules from python / python-CLE2000 API / python-LCM API
import os
import shutil
import numpy as np
import pandas as pd
import lcm 
import cle2000
import lifo
import time

# Import burnup lists from getLists.py
from getLists import getLists 

# Import custom modules for the ATRIUM-10 BWR fuel assembly
from create_geo_levels import createGeoLevels 
from tracking_operator import trackSSHGeomSALT, trackLVL1GeomSALT, trackFluxGeomSALT
from mix_handling import createLib
from burnupSchemes import burnup2LScheme, burnup1LScheme

from GEO_C_SALT import GEO_C_SALT
from MIX_C import MIX_C

save_dir = "ATRIUM10_evol_results"
if not os.path.exists(save_dir):
    os.makedirs(save_dir)


########################################################### PARAMETER SELECTION ##########################################################################
exec = True

# Options from DRAGON calculation setup.
name_geom = "AT10_24UOX" # "AT10_ASSBLY", "AT10_24UOX" for tests
burnup_points = "ATRIUM_10_S2_BU" # "ATRIUM_10_S2_BU", "CASMO-5"
### Get burnup lists from getLists.py
BU_lists = getLists(burnup_points)
BU_pts = BU_lists[0]
SelfShielding_pts = BU_lists[1]
COMPO_pts = BU_lists[2]
RATES_pts = BU_lists[3]
StepList = lcm.new('LCM','StepList')
StepList['BUList']    = np.array(BU_pts, dtype='f')
StepList['SSHList'] = np.array(SelfShielding_pts, dtype='f')
StepList['CompoList'] = np.array(COMPO_pts, dtype='f')
StepList['RatesList'] = np.array(RATES_pts, dtype='f')
StepList.close() # close without erasing
print(f"Burnup points (MWd/kgU) : {BU_pts}")
print(f"Self-shielding calculation points (MWd/kgU) : {SelfShielding_pts}")
print(f"COMPO points (MWd/kgU) : {COMPO_pts}")
print(f"RATES points (MWd/kgU) : {RATES_pts}")

### Select computational scheme used for flux calculation in depletion loop.

flux_computational_scheme = "2L_PIJ_MOC" #"2L_PIJ_MOC" # "1L_MOC" : 1 Level direct flux calculation with MOC on 295 groups
# "2L_IC_MOC" : 1st level 295gr IC flux calculation on coarse geometry, 2nd level 26gr MOC flux calculation on fine geometry
# "2L_PIJ_MOC" : 1st level 295gr PIJ flux calculation on coarse geometry, 2nd level 26gr MOC flux calculation on fine geometry



# Geometry parameters : ATRIUM-10 BWR fuel assembly
refinement_opt_name = "finest_geom" # "default", "fine1", "fine2", "fine3", "coolant_ring", "finest_on_Gd", "finest_on_Gd_coolant_ring"

refinement_options = {} 
# posible keys to "moderator" entry
# "default" -> coarse refinement of water in and out of assembly, windmills around cells
#"fine1" -> finer refinement in inner box and outer water. Ibn thermal group : still ~-4% error on fission rates in assembly corners, 1.4% RMS and 1% average,
#"fine2" -> splity in side water increased, same as fine1 apart from that
#"fine3" -> refining moderator in central box : to test
# Possible keys to "UOX_cells" entry :
# SECT_4_{n}, where n is the number of non discretized rings in cell
# coolant_ring_SECT_4_{n}, where n is the number of non discretized rings in cell : added coolant ring 
# Possible keys to "Gd_cells" entry :
# SECT_4_{n}, where n is the number of non discretized rings in cell
# coolant_ring_SECT_4_{n}, where n is the number of non discretized rings in cell : added coolant ring 


if refinement_opt_name == "finest_geom": 
    refinement_options["moderator"] = "fine2" 
    refinement_options["UOX_cells"] = "coolant_ring_SECT_4_0"
    refinement_options["Gd_cells"] = "coolant_ring_SECT_4_0"
    num_angles = 24
    line_density = 150.0
    batch = 1000
    split_water_in_moderator_box = 20
    split_moderator_box = 2
    split_water_around_moderator_box = 2
    split_intra_assembly_coolant = 4
    split_assembly_box = 2
    split_out_assembly_moderator = [10,30]


mix_numbering_option =  "number_mix_families_per_region" # "number_mix_families_per_region" , "number_mix_families_per_enrichment"


# Tracking parameters : main flux geometry
# (24, 75.0) --> SALTLC: Global RMS, maximum and average errors (%) on region volumes :  0.15415    0.89377    0.00094, below 1% on max : satisfactory.
# For SECT 4 0 with coolant ring : (24, 150)
# SALTLC: Global RMS, maximum and average errors (%) on region volumes :     0.13502     0.95490    -0.00102, below 1% on max : satisfactory.
# Could be optimized further, but not necessary for the moment.
#num_angles = 24
#line_density = 150.0
#batch = 1500 # 750 was found to be "optimal" with 20 omp procs and (24, 75.0) tracking parameters for MOC. 
reflection_type_lvl2 = "TSPC"
anisotropy_level = 4 # Level of anisotropy for the tracking, can be 1 (isotropic), 2 (linearly anisotropic), 3 (anisotropy order P_2), or 4 (anisotropy order P_3). 
solution_door = "MOC"  # Flag to indicate whether the tracking should be modified for a MOC solution, or not., else : set it to CP (Collision Probability) tracking.
moc_angular_quandrature = "GAUS"
nmu = 4  # Number of polar angles for MOC tracking : conservation ensured up to the order of P_{nmu-1} scattering : # nmu = 4 -> P3 scattering
postscript_file = f"AT10_FIG_MAIN_{refinement_opt_name}.ps"

# Tracking parameters : self-shielding geometry
num_angles_ssh = 8
line_density_ssh = 25.0
solution_door_ssh = "PIJ" # "PIJ", "IC"
if solution_door_ssh == "PIJ":
    reflection_type_ssh = "TSPC"
elif solution_door_ssh == "IC":
    reflection_type_ssh = "TISO"
else: 
    print(f"Selected self-shielding solution door is not supported")
batch_ssh = 200
postscript_file_ssh = "AT10_FIG_SSH.ps"

# Tracking/flux calculation parameters --> Level 1 geometry :
if "2L" in flux_computational_scheme:
    if "PIJ" in flux_computational_scheme:
        solution_door_lvl1 = "PIJ"
    elif "IC" in flux_computational_scheme:
        solution_door_lvl1 = "IC"
    else:
        print(f"Selected computational scheme {flux_computational_scheme} is not supported for first level calculation")
    num_angles_lvl1 = 12
    line_density_lvl1 = 45.0
    batch_lvl1 = 750
    if solution_door_lvl1 == "PIJ":
        reflection_type_lvl1 = "TSPC"
    elif solution_door_lvl1 == "IC":
        reflection_type_lvl1 = "TISO"
    else:
        print(f"Selected first level solution door is not supported")
    postscript_file_lvl1 = "AT10_FIG_LVL1.ps"
else:
    # Dummy values if 1L scheme is selected
    solution_door_lvl1 = None
    num_angles_lvl1 = None
    line_density_lvl1 = None
    batch_lvl1 = None
    reflection_type_lvl1 = "TSPC"
    postscript_file_lvl1 = None

# Parameters for the LIBRARY creation
draglib_name = "J311_295" # "endfb8r1_295" # "J311_295"
self_shielding_method = "PT"  # Method to be used for self-shielding calculations, "PT" for Mathematical Probability Tables, "SUBG" for Physical Probaility tables, "RSE" for Resonant Spectrum Expansion.
resonance_correlation = "NOCORR"  # Specify if the resonance correlation model should be applied. Only available for "RSE" and "PT". This will use a correlation model to treat reonances of U238, Pu240 and Gd157.
transport_correction = "NONE"  # Specify if a transport correction should be applied in the self-shielding calculation. "NONE", "APOL"
composition_option = "AT10_void_0"  # Specify which composition of mixes should be used for the LIBRARY creation. For now "AT10_void_0" and "AT10_void_40" are available.

# USS: call parameters :
ssh_option = "default"  # Option to specify specific groupings of self-shielding regions, to be tested. "default" is default from USS: based on LIB: data, "RSE" is for Resonant Spectrum Expansion method, and groups all U8, U5 and Zr isotopes in a single self-shielding region --> to be tested.
# EVO: call parameters :
rates_normalisation_option = "NONE" # "EDEP0"

# EDI: and COMPO: calls to save the results to a MULTICOMPO file
if transport_correction == "NONE":
    name_compo = f"_CPO_n{num_angles}_ld{int(line_density)}_n{num_angles_ssh}_ld{int(line_density_ssh)}_{reflection_type_lvl2}_{anisotropy_level}_{solution_door_ssh}_{solution_door}_{moc_angular_quandrature}_{nmu}"
else:
    name_compo = f"_CPO_n{num_angles}_ld{int(line_density)}_n{num_angles_ssh}_ld{int(line_density_ssh)}_{reflection_type_lvl2}_{anisotropy_level}_{solution_door_ssh}_{solution_door}_{moc_angular_quandrature}_{nmu}_{transport_correction}"
########################################################### END OF PARAMETER SELECTION #####################################################################

########################################################### GEO: : geometry creation ########################################################################
# Create the geometry levels according to the options selected.
if "ASSBLY" in name_geom:

    start_time = time.time()

    geo_flx, geo_lvl1, geo_ssh, connectivity_dict = createGeoLevels(name_geom, split_water_in_moderator_box, split_moderator_box, split_water_around_moderator_box, 
                                                    split_intra_assembly_coolant, split_assembly_box, split_out_assembly_moderator,
                                                    mix_numbering_option, refinement_options,
                                                    reflection_type_ssh, reflection_type_lvl1, reflection_type_lvl2)
    time_creating_geo = time.time() - start_time
    current_time = time.time()
elif name_geom == "AT10_24UOX":
    geo_flx, geo_lvl1, geo_ssh = GEO_C_SALT()


######################################################## SALT: : self-shielding geometry tracking ###########################################################

# Track the self-shielding geometry
track_lcm_ssh, track_binary_ssh, figure_ssh = trackSSHGeomSALT(geo_ssh, num_angles_ssh, line_density_ssh, reflection_type_ssh, batch_ssh, postscript_file_ssh, solution_door_ssh)

######################################################## SAL: : first level geometry tracking ###########################################################################

if "2L" in flux_computational_scheme:
    ## tracking for 1st level PIJ / IC calculation
    track_lcm_lvl1, track_binary_lvl1, figure_lvl1 = trackLVL1GeomSALT(geo_lvl1, num_angles_lvl1, line_density_lvl1, reflection_type_lvl1, batch_lvl1, postscript_file_lvl1, solution_door_lvl1)

######################################################## SALT: : main flux geometry tracking ################################################################

# Track the geometry using the SALT: module
track_lcm, track_binary, figure = trackFluxGeomSALT(geo_flx, num_angles, line_density, reflection_type_lvl2, anisotropy_level, solution_door, moc_angular_quandrature, nmu, batch, postscript_file)

###################################################### LIB: : LIBRARY creation ############################################################
# Create the LIBRARY according to the options selected.
if "ASSBLY" in name_geom:
    lib_lcm = createLib(mix_numbering_option, draglib_name, anisotropy_level, self_shielding_method, resonance_correlation, transport_correction, composition_option, rates_normalisation_option, connectivity_dict)
else:
    lib_lcm = MIX_C(draglib_name, "PT_N")
if exec:
    if "2L" in flux_computational_scheme:
        pyCOMPO = burnup2LScheme(lib_lcm, track_lcm, track_binary, track_lcm_lvl1, track_binary_lvl1, track_lcm_ssh, track_binary_ssh, StepList, name_compo, solution_door_ssh, solution_door_lvl1, name_geom)
    elif "1L" in flux_computational_scheme:
        pyCOMPO = burnup1LScheme(lib_lcm, track_lcm, track_binary, track_lcm_ssh, track_binary_ssh, StepList, name_compo, solution_door_ssh, name_geom)
    ######################################################### EXPORTS ############################################################
    if mix_numbering_option == "number_mix_families_per_region":
        numbering_save_opt = "region_num"
    elif mix_numbering_option == "number_mix_families_per_enrichment":
        numbering_save_opt = "enrich_num"
    # Save the MULTICOMPO to a specific directory
    save_dir_case = f"{save_dir}/{flux_computational_scheme}_norm_{rates_normalisation_option}/{refinement_opt_name}_{composition_option}_{draglib_name}_{self_shielding_method}_{resonance_correlation}_{numbering_save_opt}"
    if not os.path.exists(save_dir_case):
        os.makedirs(save_dir_case)
    # Save the LCM objects to the specified directory
    shutil.copyfile(name_compo, f"{save_dir_case}/{name_compo}")