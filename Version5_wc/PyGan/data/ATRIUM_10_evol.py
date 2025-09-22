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
from create_geo import createGeo
from tracking_operator import trackFluxGeomSALT
from tracking_operator import trackSSHGeomSALT
from mix_handling import createLib
from self_shielding import selfShieldingUSS
from flux_calculation import fluxCalculationMOC, fluxCalculationPIJ, fluxCalculationIC
from edi_compo import ediCompoRates, ediCompoBU, saveCpoToFile
from burnup_step import BU_A, normalizeRates
save_dir = "ATRIUM10_evol_results"
if not os.path.exists(save_dir):
    os.makedirs(save_dir)


########################################################### PARAMETER SELECTION ##########################################################################
exec = True

# Options from DRAGON calculation setup.
name_geom = "AT10_ASSBLY"
burnup_points = "ATRIUM_10_S2_BU" # "ATRIUM_10_S2_BU", "CASMO-5"
### Get burnup lists from getLists.py
BU_lists = getLists(burnup_points)
BU_pts = BU_lists[0]
SelfShielding_pts = BU_lists[1]
COMPO_pts = BU_lists[2]
RATES_pts = BU_lists[3]
print(f"Burnup points (MWd/kgU) : {BU_pts}")
print(f"Self-shielding calculation points (MWd/kgU) : {SelfShielding_pts}")
print(f"COMPO points (MWd/kgU) : {COMPO_pts}")
print(f"RATES points (MWd/kgU) : {RATES_pts}")

### Select computational scheme used for flux calculation in depletion loop.

flux_computational_scheme = "1L_MOC" # "1L_MOC" : 1 Level direct flux calculation with MOC on 295 groups
# To be implemented : 
# "2L_IC_MOC : 1st level 295gr IC flux calculation on coarse geometry, 2nd level 26gr MOC flux calculation on fine geometry
# "2L_PIJ_MOC : 1st level 295gr PIJ flux calculation on coarse geometry, 2nd level 26gr MOC flux calculation on fine geometry



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
    batch = 3000
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
reflection_type = "TSPC"
anisotropy_level = 2 # Level of anisotropy for the tracking, can be 1 (isotropic), 2 (linearly anisotropic), 3 (anisotropy order P_2), or 4 (anisotropy order P_3). 
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

# Parameters for the LIBRARY creation
draglib_name = "J311_295" # "endfb8r1_295" # "J311_295"
self_shielding_method = "RSE"  # Method to be used for self-shielding calculations, "PT" for Mathematical Probability Tables, "SUBG" for Physical Probaility tables, "RSE" for Resonant Spectrum Expansion.
resonance_correlation = "NOCORR"  # Specify if the resonance correlation model should be applied. Only available for "RSE" and "PT". This will use a correlation model to treat reonances of U238, Pu240 and Gd157.
transport_correction = "APOL"  # Specify if a transport correction should be applied in the self-shielding calculation. "NONE", "APOL"
composition_option = "AT10_void_0"  # Specify which composition of mixes should be used for the LIBRARY creation. For now "AT10_void_0" and "AT10_void_40" are available.

# USS: call parameters :
ssh_option = "default"  # Option to specify specific groupings of self-shielding regions, to be tested. "default" is default from USS: based on LIB: data, "RSE" is for Resonant Spectrum Expansion method, and groups all U8, U5 and Zr isotopes in a single self-shielding region --> to be tested.

# EDI: and COMPO: calls to save the results to a MULTICOMPO file
if transport_correction == "NONE":
    name_compo_rates = f"_CPO_R_n{num_angles}_ld{int(line_density)}_n{num_angles_ssh}_ld{int(line_density_ssh)}_{reflection_type}_{anisotropy_level}_{solution_door_ssh}_{solution_door}_{moc_angular_quandrature}_{nmu}"
    name_compo_BU = f"_CPO_BU_n{num_angles}_ld{int(line_density)}_n{num_angles_ssh}_ld{int(line_density_ssh)}_{reflection_type}_{anisotropy_level}_{solution_door_ssh}_{solution_door}_{moc_angular_quandrature}_{nmu}"
else:
    name_compo_rates = f"_CPO_R_n{num_angles}_ld{int(line_density)}_n{num_angles_ssh}_ld{int(line_density_ssh)}_{reflection_type}_{anisotropy_level}_{solution_door_ssh}_{solution_door}_{moc_angular_quandrature}_{nmu}_{transport_correction}"
    name_compo_BU = f"_CPO_BU_n{num_angles}_ld{int(line_density)}_n{num_angles_ssh}_ld{int(line_density_ssh)}_{reflection_type}_{anisotropy_level}_{solution_door_ssh}_{solution_door}_{moc_angular_quandrature}_{nmu}_{transport_correction}"
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

######################################################## SALT: : self-shielding geometry tracking ###########################################################

# Track the self-shielding geometry
track_lcm_ssh, track_binary_ssh, figure_ssh = trackSSHGeomSALT(geo_ssh, num_angles_ssh, line_density_ssh, reflection_type_ssh, batch_ssh, postscript_file_ssh, solution_door_ssh)

if exec:
    ###################################################### LIB: : LIBRARY creation ############################################################
    # Create the LIBRARY according to the options selected.
    lib_lcm = createLib(mix_numbering_option, draglib_name, anisotropy_level, self_shielding_method, resonance_correlation, transport_correction, composition_option, connectivity_dict)

    for i in range(len(BU_pts)):
        burnup_point = BU_pts[i]
        print(f"Starting burnup point at {burnup_point} MWd/kgU")
        next_burnup_point = BU_pts[i+1] if i < len(BU_pts)-1 else None
        
        if burnup_point == 0.0:
            print("Starting zero burnup calculation")
            ###################################################### USS: : self-shielding calculations ##########################################
            # Perform the self-shielding calculations
            lib_ssh = None
            lib_ssh = selfShieldingUSS(mix_numbering_option, lib_lcm, lib_ssh, track_lcm_ssh, track_binary_ssh, name_geom, ssh_option, solution_door_ssh, burnup_point)


            ###################################################### ASM: and FLU: calls ######################################################
            if flux_computational_scheme == "1L_MOC":
                keff, flux_lcm, sys_lcm = fluxCalculationMOC(track_lcm, track_binary, lib_ssh, None, None, burnup_point)

            burnup_lcm, lib_ssh = normalizeRates(flux_lcm, lib_ssh, None, track_lcm, rates_extr="EXTR", power_density=30.0, burnup_step=burnup_point)
            #################################################### EDI: and COMPO: ############################################################
            if burnup_point in RATES_pts:
                # Call EDI: and COMPO: to save the results to a MULTICOMPO file
                #compo_lcm, flux_lcm, self_shielded_microlib, track_lcm, BU
                compo_lcm_rates = ediCompoRates(None, flux_lcm, lib_ssh, track_lcm, burnup_lcm, burnup_point)
            if burnup_point in COMPO_pts:
                compo_lcm_BU = ediCompoBU(None, flux_lcm, lib_ssh, track_lcm, burnup_lcm, burnup_point)
            
            #################################################### BU_A: : burnup step ############################################################
            # Perform the burnup step using the BU_A procedure
            burnup_lcm, lib_ssh = BU_A(lib_ssh, track_lcm, burnup_lcm, flux_lcm, burnup_point, next_burnup_point, depl_sol_option="RUNG", glob_option="NOGL", sat_option="SAT", rates_extr="EXTR", edep_mode="")
        elif burnup_point > 0.0 and next_burnup_point is not None:
            print(f"Starting burnup calculation at {burnup_point} MWd/kgU")
            if burnup_point in SelfShielding_pts:
                print(f"Starting burnup calculation at {burnup_point} MWd/kgU with self-shielding")
                ###################################################### USS: : self-shielding calculations ##########################################
                # Perform the self-shielding calculations
                lib_ssh = selfShieldingUSS(mix_numbering_option, lib_lcm, lib_ssh, track_lcm_ssh, track_binary_ssh, name_geom, ssh_option, solution_door_ssh, burnup_point)
                
            ###################################################### ASM: and FLU: calls ######################################################
            if flux_computational_scheme == "1L_MOC":
                keff, flux_lcm, sys_lcm = fluxCalculationMOC(track_lcm, track_binary, lib_ssh, flux_lcm, sys_lcm, burnup_point)
                
                
            ### Reaction rates normalization at BoS of the current step
            burnup_lcm, lib_ssh = normalizeRates(flux_lcm, lib_ssh, burnup_lcm, track_lcm, rates_extr="EXTR", power_density=30.0,  burnup_step=burnup_point)
            
            
            #################################################### EDI: and COMPO: for RATES ############################################################
            if burnup_point in RATES_pts:
                # Call EDI: and COMPO: to save the results to a MULTICOMPO file
                # compo_lcm, flux_lcm, self_shielded_microlib, track_lcm, burnup_lcm, BU
                compo_lcm_rates = ediCompoRates(compo_lcm_rates, flux_lcm, lib_ssh, track_lcm, burnup_lcm, burnup_point)

            #################################################### EDI: and COMPO: for ISODENS ############################################################
            if burnup_point in COMPO_pts:
                # Call EDI: and COMPO: to save the results to a MULTICOMPO file
                compo_lcm_BU = ediCompoBU(compo_lcm_BU, flux_lcm, lib_ssh, track_lcm, burnup_lcm, burnup_point)
                
            #################################################### BU_A: : burnup step ############################################################
            
            if burnup_point in BU_pts:
                # Perform the burnup step using the BU_A procedure
                burnup_lcm, lib_ssh = BU_A(lib_ssh, track_lcm, burnup_lcm, flux_lcm, burnup_point, next_burnup_point, depl_sol_option="RUNG", glob_option="NOGL", sat_option="SAT", rates_extr="EXTR", edep_mode="")

        elif next_burnup_point is None:
            print(f"Reached last burnup point at {burnup_point} MWd/kgU, no further calculation")
            # Normlize rates and save the last COMPO: file
            burnup_lcm, lib_ssh = normalizeRates(flux_lcm, lib_ssh, burnup_lcm, track_lcm, rates_extr="EXTR", power_density=30.0,  burnup_step=burnup_point)
            if burnup_point in RATES_pts:
                compo_lcm_rates = ediCompoRates(compo_lcm_rates, flux_lcm, lib_ssh, track_lcm, burnup_lcm, burnup_point)
            if burnup_point in COMPO_pts:
                compo_lcm_BU = ediCompoBU(compo_lcm_BU, flux_lcm, lib_ssh, track_lcm, burnup_lcm, burnup_point)
            saveCpoToFile(compo_lcm_BU, compo_lcm_rates, name_compo_BU, name_compo_rates)

    ######################################################### EXPORTS ############################################################
    if mix_numbering_option == "number_mix_families_per_region":
        numbering_save_opt = "region_num"
    elif mix_numbering_option == "number_mix_families_per_enrichment":
        numbering_save_opt = "enrich_num"
    # Save the MULTICOMPO to a specific directory
    save_dir_case = f"{save_dir}/{flux_computational_scheme}/{refinement_opt_name}_{composition_option}_{draglib_name}_{self_shielding_method}_{resonance_correlation}_{numbering_save_opt}"
    if not os.path.exists(save_dir_case):
        os.makedirs(save_dir_case)
    # Save the LCM objects to the specified directory
    shutil.copyfile(name_compo_BU, f"{save_dir_case}/{name_compo_BU}")
    shutil.copyfile(name_compo_rates, f"{save_dir_case}/{name_compo_rates}")