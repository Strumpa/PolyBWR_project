## Python3 file to group Serpent2 parsing functions
import os
import numpy as np
import serpentTools as st

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


def parse_S2_ASSBLY_rates_lat_det(name_case, XS_lib_S2, fission_isotopes, n_gamma_isotopes, bu, unfold_symmetry):
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
    if unfold_symmetry:
        sym_factor = 1
    else:
        sym_factor = 0.5

    # Extracting the detector response
    ngroups = detector.detectors["_pins_2G"].tallies.shape[0]
    ncells = detector.detectors["_pins_2G"].tallies.shape[1]
    ntallies = detector.detectors["_pins_2G"].tallies.shape[2]
    tally_index_to_react = { 0: "U235_ngamma", 1 : "U238_ngamma", 2 : "Pu239_ngamma", 3 : "Pu241_ngamma",
                             4 : "Gd155_ngamma", 5: "Gd157_ngamma", 6: "Xe135_ngamma", 7 : "Sm149_ngamma", 
                             8 : "U235_fiss", 9 : "U238_fiss", 10 : "Pu239_fiss", 11 : "Pu241_fiss"}
    number_of_each_mix = {"pin1":4, "pin2":8, "pin3":10, "pin4":20, "pin5":6, "pin6":27, "pin7":14, "pin8":2}
    vol = np.pi * 0.4435 ** 2
    N_iso = {"pin1" : {"U234":  5.15910E-06, "U235":  5.67035E-04, "U238":  2.27631E-02, "O16" :  4.66705E-02}, 
            "pin2" : {"O16": 4.667480e-02, "U238": 2.257430e-02, "U234": 7.039170e-06, "U235": 7.560370e-04},
            "pin3": {"U235": 9.686590e-04, "U234": 9.163680e-06, "U238": 2.236200e-02, "O16": 4.667960e-02}, 
            "pin4": {"O16": 4.668150e-02, "U238": 2.227940e-02, "U234": 9.991530e-06, "U235": 1.051340e-03},
            "pin5": {"U234": 1.058330e-05, "U235": 1.110400e-03, "U238": 2.222040e-02, "O16": 4.668280e-02}, 
            "pin6":{"U234": 1.117530e-05, "U235": 1.169460e-03, "O16": 4.668410e-02, "U238": 2.216140e-02}, 
            "pin7":{"Gd160": 2.994740e-04, "Gd157": 2.143990e-04, "Gd158": 3.403000e-04, "Gd156": 2.804310e-04,"U238": 2.107540e-02, "O16": 4.621410e-02, "Gd155": 2.027540e-04, "U234": 9.451580e-06, "Gd154": 2.986510e-05, "U235": 9.945290e-04}, 
            "pin8": {"O16": 4.621230e-02, "U238": 2.115350e-02, "Gd156": 2.804310e-04, "Gd158": 3.403000e-04, "U235": 9.163120e-04, "Gd154": 2.986510e-05, "U234": 8.668470e-06, "Gd155": 2.027540e-04, "Gd157": 2.143990e-04, "Gd160": 2.994740e-04}
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
                print(f"Energy group {i+1}, Cell {j+1}, Reaction {tally_index_to_react[k]} : {detector.detectors['_pins_2G'].tallies[i,j,k]}")
                #print(detector.detectors['_pins_2G'].tallies[i,j,k])
                if tally_index_to_react[k] not in fission_rates[f"cell{j+1}"].keys():
                    fission_rates[f"cell{j+1}"][tally_index_to_react[k]] = []
                if tally_index_to_react[k] not in ngamma_rates[f"cell{j+1}"].keys():
                    ngamma_rates[f"cell{j+1}"][tally_index_to_react[k]] = []
                    
                if tally_index_to_react[k].endswith("ngamma"):
                    ngamma_rates[f"cell{j+1}"][tally_index_to_react[k]].append(detector.detectors['_pins_2G'].tallies[i,j,k] * N_iso[f"pin{j+1}"][tally_index_to_react[k][:-7]] * vol * sym_factor / number_of_each_mix[f"pin{j+1}"])
                elif tally_index_to_react[k].endswith("fiss"):
                    fission_rates[f"cell{j+1}"][tally_index_to_react[k]].append(detector.detectors['_pins_2G'].tallies[i,j,k] * N_iso[f"pin{j+1}"][tally_index_to_react[k][:-5]] * vol * sym_factor / number_of_each_mix[f"pin{j+1}"])


    summed_fission_rates_over_isos = sum_S2rates_over_iso(fission_rates)
    summed_ngamma_rates_over_isos = sum_S2rates_over_iso(ngamma_rates)
    fission_rates["TOT"] = summed_fission_rates_over_isos
    ngamma_rates["TOT"] = summed_ngamma_rates_over_isos


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


