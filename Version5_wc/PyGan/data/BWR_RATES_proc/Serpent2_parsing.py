## Python3 file to group Serpent2 parsing functions
import os
import numpy as np
import serpentTools as st
import re

def parse_Serpent2_lattice_det(path_to_S2, name_case, XS_lib_S2, bcond, edepmode, pcc, bu):
    """
    Lattice detector post-treatment for Serpent2 simulations
    path_to_S2 : (str) path to the Serpent2 output files
    name_case : (str) name of the Serpent2 output file
    XS_lib_S2 : (str) name of the Serpent2 XS library
    bcond : (int) boundary condition : (1) : vacuum (2) : reflective, (3) : periodic
    edepmode : (int) energy deposition mode
    pcc : (int) predictor corrector option
    bu : (int) burnup step
    """

    # energy deposition mode 0 (default) : Constant energy deposition per fission.
    # at bu = 0
    if edepmode == 0:
        detector = st.read(f"{path_to_S2}/{name_case}_{XS_lib_S2}_mc_det{bu}.m")
        res = st.read(f"{path_to_S2}/{name_case}_{XS_lib_S2}_mc_res.m")
    else:
        detector = st.read(f"{path_to_S2}/{name_case}_edep{edepmode}_mc_det{bu}.m")
        res = st.read(f"{path_to_S2}/{name_case}_edep{edepmode}_mc_res.m")
    # retrieve keff
    keff = res.resdata["absKeff"].T[0]
        
    ### _pins_2G is a lattice detector
    # For 2x2 lattice :
    #print(detector.detectors["_pins_2G"].tallies.shape) # shape = (2, 4, 10) = (energy group, cell, repsonse to mt ?)
    # Requested 2G energy groups, 4 cells/pins in the 2x2 lattice and 10 tallies / reactions (MT) in the detector
    # Reactions : MT=102 U235, U238, Pu239, Pu241, Xe135, Sm149 (microscopic rates), MT = -6 : fission "macroscopic" rates for U235, U238, Pu239, Pu241
    # For 3x3 lattice :
    # Requested 2G energy groups, 9 cells/pins in the 3x3 lattice and 10 tallies / reactions (MT) in the detector
    # Reactions : MT=102 U235, U238, Pu239, Pu241, Xe135, Sm149 (microscopic rates), MT = -6 : fission "macroscopic" rates for U235, U238, Pu239, Pu241
    # Expected detector shape : (2, 9, 10) 
    ngroups = detector.detectors["_pins_2G"].tallies.shape[0]
    ncells = detector.detectors["_pins_2G"].tallies.shape[1]
    ntallies = detector.detectors["_pins_2G"].tallies.shape[2]
    tally_index_to_react = {0: "U235_ngamma", 1: "U238_ngamma", 2: "Pu239_ngamma", 3: "Pu241_ngamma", 4: "Xe135_ngamma", 5: "Sm149_ngamma", 6: "U235_fiss", 7: "U238_fiss", 8: "Pu239_fiss", 9: "Pu241_fiss"}
    if name_case == "AT10_2x2_UOX":
                #      C1           C2           C2             C4
        N_U235 = [5.67035E-04, 7.560370E-04, 7.560370E-04, 1.051340E-03] # U235 densities in each cell, in 1. positive x, 2. positive y order
    elif name_case == "bench_3x3_UOX":
                       
        N_U235 = [5.67035E-04, 7.560370E-04, 9.686590E-04, #C1, C2, C3  
                  7.560370E-04, 1.051340E-03, 5.67035E-04, #C2, C4, C1 
                  9.686590E-04, 5.67035E-04, 1.169460E-03] #C3, C1, C6
    reaction_rates = {}
    U235_fiss_rates_S2_lat_det = {}
    # Extracting the detector response
    # Strange that responses are not 0 for Pu isotopes fission rates at bu = 0 ?
    for i in range(ngroups):
        for j in range(ncells):
            if f"cell{j+1}" not in U235_fiss_rates_S2_lat_det.keys():
                U235_fiss_rates_S2_lat_det[f"cell{j+1}"] = {}
            for k in range(ntallies):
                print(f"Energy group {i+1}, Cell {j+1}, Reaction {tally_index_to_react[k]} : {detector.detectors['_pins_2G'].tallies[i,j,k]}")
                #print(detector.detectors['_pins_2G'].tallies[i,j,k])
                if tally_index_to_react[k] not in reaction_rates.keys():
                    reaction_rates[tally_index_to_react[k]] = {}
                reaction_rates[tally_index_to_react[k]][f"cell{j+1}_G{i+1}"] = detector.detectors['_pins_2G'].tallies[i,j,k]    
                # For now extracting only the fission rates of U235
                if k == 6:
                    #reaction_rates[f"U235_fiss_cell{j+1}_G{i+1}"] = detector.detectors['_pins_2G'].tallies[i,j,k]*N_U235[f"C{j+1}"]
                    U235_fiss_rates_S2_lat_det[f"cell{j+1}"][f"G{i+1}"] = detector.detectors['_pins_2G'].tallies[i,j,k]*N_U235[j]

    # Symmetrizing the fission rates for the 2x2 lattice
    if name_case == "AT10_2x2_UOX":
        U235_fiss_rates_S2_lat_det["cell2"]["G1"] = (U235_fiss_rates_S2_lat_det["cell2"]["G1"] + U235_fiss_rates_S2_lat_det["cell3"]["G1"]) / 2
        U235_fiss_rates_S2_lat_det["cell3"]["G1"] = U235_fiss_rates_S2_lat_det["cell2"]["G1"]
        U235_fiss_rates_S2_lat_det["cell2"]["G2"] = (U235_fiss_rates_S2_lat_det["cell2"]["G2"] + U235_fiss_rates_S2_lat_det["cell3"]["G2"]) / 2
        U235_fiss_rates_S2_lat_det["cell3"]["G2"] = U235_fiss_rates_S2_lat_det["cell2"]["G2"]
    
    return keff, U235_fiss_rates_S2_lat_det


def parse_Serpent2_material_det(name_case, XS_lib_S2, bcond, edepmode, pcc, bu):
    """
    Material detector post-treatment for Serpent2 simulations
    path_to_S2 : (str) path to the Serpent2 output files
    name_case : (str) name of the Serpent2 output file
    XS_lib_S2 : (str) name of the Serpent2 XS library
    bcond : (int) boundary condition : (1) : vacuum (2) : reflective, (3) : periodic
    edepmode : (int) energy deposition mode
    pcc : (int) predictor corrector option
    bu : (int) burnup step
    """
    # energy deposition mode 0 (default) : Constant energy deposition per fission.
    # at bu = 0
    if edepmode == 0:
        detector = st.read(f"{os.environ['SERPENT_RESULTS']}/{name_case}_{XS_lib_S2}_mc_det{bu}.m")
        res = st.read(f"{os.environ['SERPENT_RESULTS']}/{name_case}_{XS_lib_S2}_mc_res.m")
    else:
        detector = st.read(f"{os.environ['SERPENT_RESULTS']}/{name_case}_edep{edepmode}_mc_det{bu}.m")
        res = st.read(f"{os.environ['SERPENT_RESULTS']}/{name_case}_edep{edepmode}_mc_res.m")
    # retrieve keff
    keff = res.resdata["absKeff"].T[0]

    # Extracting the detector response

    n_groups = detector.detectors['_C1_2G'].tallies.shape[0]
    n_reactions = detector.detectors['_C1_2G'].tallies.shape[1]
    tally_index_to_react = {0: "U235_ngamma", 1: "U238_ngamma", 2: "Pu239_ngamma", 3: "Pu241_ngamma", 4: "Xe135_ngamma", 5: "Sm149_ngamma", 6: "U235_fiss", 7: "U238_fiss", 8: "Pu239_fiss", 9: "Pu241_fiss"}
    reaction_rates = {}
    vol = np.pi * 0.4435 ** 2
    U235_fiss_rates_S2_mat_det = {}
    
    if name_case == "AT10_2x2_UOX":
        """
        1 2 
        2 4 
        """
        N_U235 = {"C1" : 5.67035E-04, "C2" : 7.560370E-04, "C4": 1.051340E-03}
        vol_factor = {"C1": 1, "C2": 2, "C4": 1}
        cells_in_problem = ["C1", "C2", "C4"]
        reg_idx_to_cell = {0 : "C1", 1 : "C2", 
                           2 : "C2", 3 : "C4"}
        n_reg = 4


    elif name_case == "bench_3x3_UOX":       
        """
        1 2 3 
        2 4 1
        3 1 6
        """
        N_U235 = {"C1" : 5.67035E-04, "C2" : 7.560370E-04, "C3" : 9.686590E-04, "C4": 1.051340E-03, "C6" : 1.169460E-03} 
        vol_factor = {"C1": 3, "C2": 2, "C3": 2, "C4": 1, "C6": 1}
        cells_in_problem = ["C1", "C2", "C3", "C4", "C6"]
        reg_idx_to_cell = {0 : "C1", 1 : "C2", 2 : "C3", 
                           3 : "C2", 4 : "C4", 5 : "C1", 
                           6 : "C3", 7 : "C1", 8 : "C6"}
        n_reg = 9
    
    elif name_case == "AT10_3x3_UOX_Gd":
        """
        1 2 3 
        2 4 7 
        3 7 6
        """
        N_U235 = {"C1" : 5.67035E-04, "C2" : 7.560370E-04, "C3" : 9.686590E-04, "C4": 1.051340E-03, "C6" : 1.169460E-03, "C7": 9.945290E-04} 
        vol_factor = {"C1": 1, "C2": 2, "C3": 2, "C4": 1, "C6": 1, "C7": 2}
        cells_in_problem = ["C1", "C2", "C3", "C4", "C6","C7"]
        reg_idx_to_cell = {0 : "C1", 1 : "C2", 2 : "C3", 
                           3 : "C2", 4 : "C4", 5 : "C7", 
                           6 : "C3", 7 : "C7", 8 : "C6"}
        n_reg = 9
    
    for cell in cells_in_problem:
        cell_scores = detector.detectors[f'_{cell}_2G'].tallies
        for gr in range(n_groups):
            if f"{cell}" not in U235_fiss_rates_S2_mat_det.keys():
                U235_fiss_rates_S2_mat_det[f"{cell}"] = {}
            for k in range(n_reactions):
                print(f"Cell {cell}, Energy group {gr+1}, Reaction {tally_index_to_react[k]} : {cell_scores[gr,k]}")
                if tally_index_to_react[k] not in reaction_rates.keys():
                    reaction_rates[tally_index_to_react[k]] = {}
                reaction_rates[tally_index_to_react[k]][f"{cell}_G{gr+1}"] = cell_scores[gr,k]
                # For now extracting only the fission rates of U235
                if k == 6:
                    #reaction_rates[f"U235_fiss_cell{j+1}_G{i+1}"] = detector.detectors['_pins_2G'].tallies[i,j,k]*N_U235[f"C{j+1}"]
                    U235_fiss_rates_S2_mat_det[f"{cell}"][f"G{gr+1}"] = cell_scores[gr,k] * N_U235[cell] * vol / vol_factor[cell]
    
    # Symmetrizing the fission rates for the 3x3 lattice
        
    symmetrized_U235_fiss_rates_S2_mat_det = {}
    for i in range(n_reg):
        cell = reg_idx_to_cell[i]
        if cell not in symmetrized_U235_fiss_rates_S2_mat_det.keys():
            symmetrized_U235_fiss_rates_S2_mat_det[i] = {}
        for gr in range(n_groups):
            symmetrized_U235_fiss_rates_S2_mat_det[i][f"G{gr+1}"] = U235_fiss_rates_S2_mat_det[reg_idx_to_cell[i]][f"G{gr+1}"]
    print(f"U235_fiss_rates_S2_mat_det = {U235_fiss_rates_S2_mat_det}")
    print(f"symmetrized_U235_fiss_rates_S2_mat_det = {symmetrized_U235_fiss_rates_S2_mat_det}")
    
    
    return keff, symmetrized_U235_fiss_rates_S2_mat_det


def parse_S2_ASSBLY_rates(name_case, XS_lib_S2, fission_isotopes, n_gamma_isotopes, bu, unfold_symmetry):
    """
    Serpent2 assembly detector post-treatment for Assembly reaction rates :
    name_case : (str) name of the Serpent2 case
    XS_lib_S2 : (str) name of the Serpent2 XS library
    fission_isotopes : (list) list of fission isotopes
    n_gamma_isotopes : (list) list of (n,gamma) isotopes
    bu : (int) burnup step
    unfold_symmetry : (bool) if True, unfold the symmetry of the assembly : create a 10x10 diagonally symmetric grid

    returns : reaction rates for fission and (n,gamma) reactions for the specified isotopes
    """
    detector = st.read(f"{os.environ['SERPENT_RESULTS']}/{name_case}_{XS_lib_S2}_mc_det{bu}.m")
    res = st.read(f"{os.environ['SERPENT_RESULTS']}/{name_case}_{XS_lib_S2}_mc_res.m")
    keff = res.resdata["absKeff"].T[0]

    print(f"keff = {keff}")
    print(detector.detectors.keys())
    if unfold_symmetry:
        sym_factor = 1
    else:
        sym_factor = 0.5

    # Extracting the detector response
    n_groups = detector.detectors['_pin1_2G'].tallies.shape[0]
    n_reactions = detector.detectors['_pin1_2G'].tallies.shape[1]
    tally_index_to_react = { 0: "U235_ngamma", 1 : "U238_ngamma", 2 : "Pu239_ngamma", 3 : "Pu241_ngamma",
                             4 : "Gd155_ngamma", 5: "Gd157_ngamma", 6: "Xe135_ngamma", 7 : "Sm149_ngamma", 
                             8 : "U235_fiss", 9 : "U238_fiss", 10 : "Pu239_fiss", 11 : "Pu241_fiss"}
    number_of_each_mix = {"pin1":4, "pin2":8, "pin3":10, "pin4":20, "pin5":6, "pin6":27, "pin7":14, "pin8":2}
    vol = np.pi * 0.4435 ** 2
    N_iso = {"pin1" : {"U234":  5.15910E-06,
                    "U235":  5.67035E-04,
                    "U238":  2.27631E-02,
                    "O16" :  4.66705E-02},
            "pin2" : {"O16": 4.667480e-02,
                    "U238": 2.257430e-02,
                    "U234": 7.039170e-06,
                    "U235": 7.560370e-04},
            "pin3": {"U235": 9.686590e-04,
                    "U234": 9.163680e-06,
                    "U238": 2.236200e-02,
                    "O16": 4.667960e-02}, 
            "pin4": {"O16": 4.668150e-02,
                    "U238": 2.227940e-02,
                    "U234": 9.991530e-06,
                    "U235": 1.051340e-03},
            "pin5": {"U234": 1.058330e-05,
                    "U235": 1.110400e-03,
                    "U238": 2.222040e-02,
                    "O16": 4.668280e-02}, 
            "pin6":{"U234": 1.117530e-05,
                    "U235": 1.169460e-03,
                    "O16": 4.668410e-02,
                    "U238": 2.216140e-02}, 
            "pin7":{"Gd160": 2.994740e-04,
                    "Gd157": 2.143990e-04,
                    "Gd158": 3.403000e-04,
                    "Gd156": 2.804310e-04,
                    "U238": 2.107540e-02,
                    "O16": 4.621410e-02,
                    "Gd155": 2.027540e-04,
                    "U234": 9.451580e-06,
                    "Gd154": 2.986510e-05,
                    "U235": 9.945290e-04}, 
            "pin8": {"O16": 4.621230e-02,
                    "U238": 2.115350e-02,
                    "Gd156": 2.804310e-04,
                    "Gd158": 3.403000e-04,
                    "U235": 9.163120e-04,
                    "Gd154": 2.986510e-05,
                    "U234": 8.668470e-06,
                    "Gd155": 2.027540e-04,
                    "Gd157": 2.143990e-04,
                    "Gd160": 2.994740e-04}
        }
    fission_rates = {}
    ngamma_rates={}


    for pin in N_iso.keys():
        pin_scores = detector.detectors[f'_{pin}_2G'].tallies
        
        if f"C{pin[-1]}" not in fission_rates.keys():
            fission_rates[f"C{pin[-1]}"] = {}
        if f"C{pin[-1]}" not in ngamma_rates.keys():
            ngamma_rates[f"C{pin[-1]}"] = {}
        ngamma_rates[f"C{pin[-1]}"]["U235"] = pin_scores[:,0] * N_iso[pin]["U235"] * vol * sym_factor / number_of_each_mix[pin]
        ngamma_rates[f"C{pin[-1]}"]["U238"] = pin_scores[:,1] * N_iso[pin]["U238"] * vol * sym_factor / number_of_each_mix[pin] 
        if pin in ["pin7", "pin8"]:
            ngamma_rates[f"C{pin[-1]}"]["Gd155"] = pin_scores[:,4] * N_iso[pin]["Gd155"] * vol * sym_factor / number_of_each_mix[pin]
            ngamma_rates[f"C{pin[-1]}"]["Gd157"] = pin_scores[:,5] * N_iso[pin]["Gd157"] * vol * sym_factor / number_of_each_mix[pin]
        fission_rates[f"C{pin[-1]}"]["U235"] = pin_scores[:,8] * N_iso[pin]["U235"] * vol * sym_factor / number_of_each_mix[pin]
        fission_rates[f"C{pin[-1]}"]["U238"] = pin_scores[:,9] * N_iso[pin]["U238"] * vol * sym_factor / number_of_each_mix[pin]

    summed_fission_rates_over_isos = sum_S2rates_over_iso(fission_rates)
    summed_ngamma_rates_over_isos = sum_S2rates_over_iso(ngamma_rates)
    fission_rates["TOT"] = summed_fission_rates_over_isos
    ngamma_rates["TOT"] = summed_ngamma_rates_over_isos


    return keff, fission_rates, ngamma_rates


def parse_S2_ASSBLY_rates_lat_det(name_case, XS_lib_S2, fission_isotopes, n_gamma_isotopes, bu):
    """
    Serpent2 assembly detector post-treatment for Assembly reaction rates :
    use lattice detector definition
    name_case : (str) name of the Serpent2 case
    XS_lib_S2 : (str) name of the Serpent2 XS library
    fission_isotopes : (list) list of fission isotopes
    n_gamma_isotopes : (list) list of (n,gamma) isotopes
    bu : (int) burnup step
    unfold_symmetry : (bool) if True, unfold the symmetry of the assembly : create a 10x10 diagonally symmetric grid

    returns : reaction rates for fission and (n,gamma) reactions for the specified isotopes
    """
    detector = st.read(f"{os.environ['SERPENT_RESULTS']}/{name_case}_{XS_lib_S2}_mc_det{bu}.m")
    res = st.read(f"{os.environ['SERPENT_RESULTS']}/{name_case}_{XS_lib_S2}_mc_res.m")
    keff = res.resdata["absKeff"].T[0]

    print(f"keff = {keff}")
    print(detector.detectors.keys())

    # Extracting the detector response
    ngroups = detector.detectors["_lattice_2G"].tallies.shape[0] # 2 = n_groups
    ncells = detector.detectors["_lattice_2G"].tallies.shape[1] # 144 = 12*12 = ncells in lattice but pins are only 10 * 10 : water "pins" around fuel lattice.
    ntallies = detector.detectors["_lattice_2G"].tallies.shape[2] # 7 reactions
    print(f"lattice detector shape : {detector.detectors['_lattice_2G'].tallies.shape}")
    tally_index_to_react = { 0: "U235_ngamma", 1 : "U238_ngamma",
                             2 : "Gd155_ngamma", 3: "Gd157_ngamma", 
                             4 : "U234_fiss", 5 : "U235_fiss", 6 : "U238_fiss"}
    vol = np.pi * 0.4435 ** 2
    N_iso = {"pin1" : {"U234":  5.15910E-06, "U235":  5.67035E-04, "U238":  2.27631E-02, "O16" :  4.66705E-02, "Pu239":0.0, "Pu241":0.0, "Gd155":0.0, "Gd157":0.0, "Xe135":0.0, "Sm149":0.0}, 
            "pin2" : {"O16": 4.667480e-02, "U238": 2.257430e-02, "U234": 7.039170e-06, "U235": 7.560370e-04, "Pu239":0.0, "Pu241":0.0, "Gd155":0.0, "Gd157":0.0, "Xe135":0.0, "Sm149":0.0},
            "pin3": {"U235": 9.686590e-04, "U234": 9.163680e-06, "U238": 2.236200e-02, "O16": 4.667960e-02, "Pu239":0.0, "Pu241":0.0, "Gd155":0.0, "Gd157":0.0, "Xe135":0.0, "Sm149":0.0}, 
            "pin4": {"O16": 4.668150e-02, "U238": 2.227940e-02, "U234": 9.991530e-06, "U235": 1.051340e-03, "Pu239":0.0, "Pu241":0.0, "Gd155":0.0, "Gd157":0.0, "Xe135":0.0, "Sm149":0.0},
            "pin5": {"U234": 1.058330e-05, "U235": 1.110400e-03, "U238": 2.222040e-02, "O16": 4.668280e-02, "Pu239":0.0, "Pu241":0.0, "Gd155":0.0, "Gd157":0.0, "Xe135":0.0, "Sm149":0.0}, 
            "pin6":{"U234": 1.117530e-05, "U235": 1.169460e-03, "O16": 4.668410e-02, "U238": 2.216140e-02, "Pu239":0.0, "Pu241":0.0, "Gd155":0.0, "Gd157":0.0, "Xe135":0.0, "Sm149":0.0}, 
            "pin7":{"Gd160": 2.994740e-04, "Gd157": 2.143990e-04, "Gd158": 3.403000e-04, "Gd156": 2.804310e-04,"U238": 2.107540e-02, "O16": 4.621410e-02, "Gd155": 2.027540e-04, "U234": 9.451580e-06, "Gd154": 2.986510e-05, "U235": 9.945290e-04, "Pu239":0.0, "Pu241":0.0, "Xe135":0.0, "Sm149":0.0}, 
            "pin8": {"O16": 4.621230e-02, "U238": 2.115350e-02, "Gd156": 2.804310e-04, "Gd158": 3.403000e-04, "U235": 9.163120e-04, "Gd154": 2.986510e-05, "U234": 8.668470e-06, "Gd155": 2.027540e-04, "Gd157": 2.143990e-04, "Gd160": 2.994740e-04, "Pu239":0.0, "Pu241":0.0, "Xe135":0.0, "Sm149":0.0}
        }
    fission_rates = {}
    ngamma_rates={}

    for i in range(ngroups):
        for j in range(ncells):
            if f"cell{j+1}" not in fission_rates.keys():
                fission_rates[f"cell{j+1}"] = {}
            if f"cell{j+1}" not in ngamma_rates.keys():
                ngamma_rates[f"cell{j+1}"] = {}
            for k in range(ntallies):
                print(f"Energy group {i+1}, Cell {j+1}, Reaction {tally_index_to_react[k]} : {detector.detectors['_lattice_2G'].tallies[i,j,k]}")
                #print(detector.detectors['_pins_2G'].tallies[i,j,k])
                if tally_index_to_react[k] not in fission_rates[f"cell{j+1}"].keys():
                    fission_rates[f"cell{j+1}"][tally_index_to_react[k]] = []
                if tally_index_to_react[k] not in ngamma_rates[f"cell{j+1}"].keys():
                    ngamma_rates[f"cell{j+1}"][tally_index_to_react[k]] = []
                    
                if tally_index_to_react[k].endswith("ngamma"):
                    ngamma_rates[f"cell{j+1}"][tally_index_to_react[k]].append(detector.detectors['_lattice_2G'].tallies[i,j,k] * N_iso[f"pin{j+1}"][tally_index_to_react[k][:-7]] * vol)
                elif tally_index_to_react[k].endswith("fiss"):
                    fission_rates[f"cell{j+1}"][tally_index_to_react[k]].append(detector.detectors['_lattice_2G'].tallies[i,j,k] * N_iso[f"pin{j+1}"][tally_index_to_react[k][:-5]] * vol)


    summed_fission_rates_over_isos = sum_S2rates_over_iso(fission_rates)
    summed_ngamma_rates_over_isos = sum_S2rates_over_iso(ngamma_rates)
    fission_rates["TOT"] = summed_fission_rates_over_isos
    ngamma_rates["TOT"] = summed_ngamma_rates_over_isos

    n_groups_fine = detector.detectors['spectrum_295G'].tallies.shape[0]


    return keff, fission_rates, ngamma_rates


def parse_S2_pin_mat_det(name_case, XS_lib_S2, fission_isotopes, ngamma_isotopes, bu):
    """
    Serpent2 assembly detector post-treatment for Assembly reaction rates :
    use lattice detector definition
    name_case : (str) name of the Serpent2 case
    XS_lib_S2 : (str) name of the Serpent2 XS library
    fission_isotopes : (list) list of fission isotopes
    ngamma_isotopes : (list) list of (n,gamma) isotopes
    bu : (int) burnup step
    unfold_symmetry : (bool) if True, unfold the symmetry of the assembly : create a 10x10 diagonally symmetric grid

    returns : reaction rates for fission and (n,gamma) reactions for the specified isotopes
    """
    iso_code_to_num = {"U235":0, "U238": 2, "Pu239": 5, "Pu241": 7, 
                        "Xe135":18, "Sm149": 30, "Gd154":21, "Gd155":22, "Gd156":23, "Gd157":24}
    detectorFile = st.read(f"{os.environ['SERPENT_RESULTS']}/{name_case}_{XS_lib_S2}_inp_det{bu}.m")
    resultsFile = st.read(f"{os.environ['SERPENT_RESULTS']}/{name_case}_{XS_lib_S2}_inp_res.m")
    keff = resultsFile.resdata["absKeff"].T[0]
    
    # try reading depletion 
    try:
        depletionFile = st.read(f"{os.environ['SERPENT_RESULTS']}/{name_case}_{XS_lib_S2}_BU_inp_dep.m")
    except:
        depletionFile = None


    print(f"keff = {keff}")
    print(detectorFile.detectors.keys())
    
    N_iso = {"24UOX" : {"U234":  5.15910E-06, "U235":  5.67035E-04, "U238":  2.27631E-02, "O16" :  4.66705E-02, "Pu239":0.0, "Pu241":0.0, "Gd155":0.0, "Gd157":0.0, "Xe135":0.0, "Sm149":0.0}, 
        "32UOX" : {"O16": 4.667480e-02, "U238": 2.257430e-02, "U234": 7.039170e-06, "U235": 7.560370e-04, "Pu239":0.0, "Pu241":0.0, "Gd155":0.0, "Gd157":0.0, "Xe135":0.0, "Sm149":0.0},
        "42UOX": {"U235": 9.686590e-04, "U234": 9.163680e-06, "U238": 2.236200e-02, "O16": 4.667960e-02, "Pu239":0.0, "Pu241":0.0, "Gd155":0.0, "Gd157":0.0, "Xe135":0.0, "Sm149":0.0}, 
        "45UOX": {"O16": 4.668150e-02, "U238": 2.227940e-02, "U234": 9.991530e-06, "U235": 1.051340e-03, "Pu239":0.0, "Pu241":0.0, "Gd155":0.0, "Gd157":0.0, "Xe135":0.0, "Sm149":0.0},
        "48UOX": {"U234": 1.058330e-05, "U235": 1.110400e-03, "U238": 2.222040e-02, "O16": 4.668280e-02, "Pu239":0.0, "Pu241":0.0, "Gd155":0.0, "Gd157":0.0, "Xe135":0.0, "Sm149":0.0}, 
        "50UOX":{"U234": 1.117530e-05, "U235": 1.169460e-03, "O16": 4.668410e-02, "U238": 2.216140e-02, "Pu239":0.0, "Pu241":0.0, "Gd155":0.0, "Gd157":0.0, "Xe135":0.0, "Sm149":0.0}, 
        "45Gd":{"Gd160": 2.994740e-04, "Gd157": 2.143990e-04, "Gd158": 3.403000e-04, "Gd156": 2.804310e-04,"U238": 2.107540e-02, "O16": 4.621410e-02, "Gd155": 2.027540e-04, "U234": 9.451580e-06, "Gd154": 2.986510e-05, "U235": 9.945290e-04, "Pu239":0.0, "Pu241":0.0, "Xe135":0.0, "Sm149":0.0}, 
        "42Gd": {"O16": 4.621230e-02, "U238": 2.115350e-02, "Gd156": 2.804310e-04, "Gd158": 3.403000e-04, "U235": 9.163120e-04, "Gd154": 2.986510e-05, "U234": 8.668470e-06, "Gd155": 2.027540e-04, "Gd157": 2.143990e-04, "Gd160": 2.994740e-04, "Pu239":0.0, "Pu241":0.0, "Xe135":0.0, "Sm149":0.0}
    }


    # Extracting the detector response
    ngroups = detectorFile.detectors["_pin_pos_1_2G"].tallies.shape[0] # 2 = n_groups
    ntallies = detectorFile.detectors["_pin_pos_1_2G"].tallies.shape[1] # 7 reactions
    vol = np.pi * 0.4435 ** 2
    tally_index_to_react = { 0: "U234_fission", 1: "U235_fission", 2:"U236_fission", 3:"U238_fission", 4:"Pu239_fission", 5:"Pu241_fission",
                            6:"U235_ngamma", 7:"U238_ngamma", 8:"Pu239_ngamma", 9:"Pu241_ngamma", 
                            10:"Xe135_ngamma", 11:"Sm149_ngamma", 12:"Gd154_ngamma", 13:"Gd155_ngamma", 14:"Gd156_ngamma", 15:"Gd157_ngamma", 16:"Gd158_ngamma", 17:"Gd160_ngamma"}

    pin_positions = [1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 
                        11, 12, 13, 14, 15, 16, 17, 18, 19,
                            20, 21, 22, 23, 24, 25, 26, 27,
                                28, 29, 30, 31, 32, 33, 34,
                                                38, 39, 40, 
                                                43, 44, 45,
                                                47, 48, 49,
                                                50, 51, 52,
                                                    53, 54,
                                                        55]
    pos_on_diag = [1, 11, 20, 28, 50, 53, 55]
    npos = max(pin_positions)
    fission_rates = np.zeros((ngroups, npos))
    ngamma_rates = np.zeros((ngroups, npos))
    n_isotopes_per_mat = {}
    groups = {}
    pattern = re.compile(r"^(\d{2}(?:UOX|Gd))_([A-Za-z])_(\d+)$")
    # Identify rings in groups of identical 'initial' materials ie same pin position, same initial enrichment.

    for id_nb, material_name in enumerate(depletionFile):
        print(f"Depletion File Material name {material_name}")
        m = pattern.match(material_name)
        n_iso = 0
        if m:
            material = m.group(1)   # ie material code example : "24UOX"
            ring_number = m.group(2)   # ie "A", "B", "C" or "D" for UOX and + "E", "F" for Gd
            pos  = int(m.group(3))  # e.g. 1
            key = (material, pos)
            groups.setdefault(key, []).append(material_name)
            
    #print(groups)
    
    for iso in fission_isotopes:
        n_iso = 0
        for key, material_names in groups.items():
            #print(f"Processing isotope {iso} for material {key}")
            material, pos = key
            n_iso = 0
            #print(f"Material {material}, Position {pos}, Material names : {material_names}")
            for material_name in material_names:
                # Example of material_name : "24UOX_A_1" or "45Gd_D_31"
                if key not in n_isotopes_per_mat.keys():
                    n_isotopes_per_mat[key] = {}
                n_iso_ring = depletionFile[material_name]['adens'][iso_code_to_num[iso]][bu]*depletionFile[material_name]['volume'][bu]
                #print(f"Material {material_name}, Isotope {iso}, Density {depletionFile[material_name]['adens'][iso_code_to_num[iso]][bu]}, Volume {depletionFile[material_name]['volume'][bu]}")
                n_iso += n_iso_ring
            #print(f"Total number density of isotope {iso} in material {material} at position {pos} : {n_iso/vol} atoms/b-cm")
            #print(f"Comparing to pin definition : {N_iso[material][iso]} atoms/b-cm")
            n_isotopes_per_mat[key][iso] = n_iso/vol
            #print(f"Material {material_name}, Isotope {iso}, Density {n_iso}, pos {pos}")
        print(n_isotopes_per_mat)
    # extract fision rates
    for pos_idx in pin_positions:
        detector_name = f"_pin_pos_{pos_idx}_2G"
        n_groups = detectorFile.detectors[detector_name].tallies.shape[0]
        n_reactions = detectorFile.detectors[detector_name].tallies.shape[1]
        for g in range(n_groups):
            for r in range(n_reactions):
                rate = detectorFile.detectors[detector_name].tallies[g, r]
                reaction = tally_index_to_react[r]
                isotope = reaction.split("_")[0]
                # Identify material corresponding to the pin position
                for key in n_isotopes_per_mat.keys():
                    material, pos = key
                    if pos == pos_idx:
                        if isotope in fission_isotopes: #or isotope in ngamma_isotopes:
                            ### Sanity check : compare N_iso from pin definition and from depletion file
                            symmetry_factor = 1 if pos in pos_on_diag else 2
                            N_iso_pin = N_iso[material][isotope]
                            N_iso_dep = n_isotopes_per_mat[key][isotope] / symmetry_factor
                            if abs(N_iso_pin - N_iso_dep)/N_iso_pin > 0.1:
                                print(f"Warning : large discrepancy between N_iso from pin definition and from depletion file for isotope {isotope} in material {material} at position {pos_idx} : N_iso_pin = {N_iso_pin}, N_iso_dep = {N_iso_dep}")
                            print(f"Material {material}, Position {pos_idx}, Isotope {isotope}, Rate {rate}, N_iso_pin = {N_iso_pin}, N_iso_dep = {N_iso_dep}")
                            ### Now fill the rates arrays
                            if reaction.endswith("fission") and isotope in fission_isotopes:
                                fission_rates[g, pos_idx-1] += rate / symmetry_factor * N_iso_dep 
                            elif reaction.endswith("ngamma") and isotope in ngamma_isotopes:
                                ngamma_rates[g, pos_idx-1] += rate / symmetry_factor * N_iso_dep 
                            print(f"Group {g+1}, Cell {pos_idx}, Reaction {reaction}, Rate: {rate}, N_iso_dep = {N_iso_dep}")
        #for gr in range(ngroups):
        

    print(n_isotopes_per_mat)

    

    return keff, fission_rates, ngamma_rates


def sum_S2rates_over_iso(rates_dict):
    """
    rates_dict is a nested dictionary with the structure:
    {
        "mix1": {
            "iso1": [val1, val2, ...],
            "iso2": [...],
            ...
        },
        "mix2": {
            ...
        },
        ...
    }

    This function sums the values over isotopes while preserving the mix structure,
    returning:
    {
        "mix1": [sum_val1, sum_val2, ...],
        "mix2": [...],
        ...
    }
    """
    result = {}

    for mix, iso_dict in rates_dict.items():
        sum_vector = None
        for values in iso_dict.values():
            if sum_vector is None:
                sum_vector = [0.0] * len(values)
            for i, val in enumerate(values):
                sum_vector[i] += val
        result[mix] = sum_vector

    return result


