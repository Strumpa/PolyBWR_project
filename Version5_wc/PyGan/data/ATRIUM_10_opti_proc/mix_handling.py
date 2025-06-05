# # PyGan deck to handle mix numbering and microlib creation for DRAGON5 calculations.
# Author : R. Guasch
# Date : 2025-06-03
# createLib reads input from main PyGan procedure, creates LIBRARY according to :
#   - composition type
#   - mix numbering option selected
#   - anisotropy level selected
#   - self-shielding method selected
#   - potential transport correction
import lifo 
import cle2000
from collections import defaultdict


def createLib(mix_numbering_option, draglib_name, anisotropy_level, self_shielding_method, resonance_correlation, transport_correction, composition_option, mix_connectivity=None):
    """
    mix_numbering_option : (str)
        Option for mix numbering, can be "number_mix_families_per_enrichment" or "number_mix_families_per_region" or "number_mix_families_per_selfshielding_region".
    draglib_name : (str)
        Name of the draglib to read evaluation data from (should be accessed via .access file).
    anisotropy_level : (int)
        Level of anisotropy for the mix handling, can be 1 (isotropic), 2 (linearly anisotropic), 3 (anisotropy order P_2), or 4 (anisotropy order P_3).
    self_shielding_method : (str)
        Method to be used for self-shielding calculations, "PT" for Mathematical Probability Tables, "SUBG" for Physical Probaility tables, "RSE" for Resonant Spectrum Expansion,
        All subgroup methods that can be treated via the USS: module
    resonance_correlation : (str)
        Specify if the resonance correlation model should be applied. Only available for "RSE" and "PT".
        This will use a correlation model to treat reonances of U238, Pu240 and Gd157.
    transport_correction : (str)
        Specify the transport correction to be applied in isotropic cases,  "APOL", or "NONE" in this context.
    composition_option : (str)
        Specify which composition of mixes should be used for the LIBRARY creation. For now "AT10_void_0" is the only option available.
    """

    if mix_numbering_option not in ["number_mix_families_per_enrichment", "number_mix_families_per_region", "number_mix_families_per_selfshielding_region"]:
        raise ValueError("Invalid mix numbering option selected.")
    if anisotropy_level not in [1, 2, 3, 4]:
        raise ValueError("Invalid anisotropy level selected. Must be 1, 2, 3, or 4.")
    if self_shielding_method not in ["PT", "SUBG", "RSE"]:
        raise ValueError("Invalid self-shielding method selected. Must be 'PT', 'SUBG', or 'RSE'.")
    if transport_correction not in ["CTRA", "APOL", "NONE", ""]:
        raise ValueError("Invalid transport correction selected. Must be 'APOL', 'NONE', or ''.")
    if draglib_name == "":
        raise ValueError("Draglib name cannot be empty.")
    if not isinstance(draglib_name, str):
        raise TypeError("Draglib name must be a string.")
    

    if mix_numbering_option == "number_mix_families_per_enrichment":
        # Use Mix_ASSBLY1.c2m procedure
        print("Creating LIBRARY with mix numbering per enrichment.")
        myLifo = lifo.new()
        myLifo.pushEmpty("LIBRARY", "LCM")
        myLifo.push(draglib_name)
        myLifo.push(composition_option)
        myLifo.push(self_shielding_method)
        myLifo.push(resonance_correlation)
        myLifo.push(anisotropy_level)
        myLifo.push(transport_correction)

        # Create a cle2000 object to handle the library creation
        lib_proc = cle2000.new('MIX_ASSBLY1', myLifo, 1)
        # Execute the library creation procedure
        lib_proc.exec()

        # Recover the results from the LIFO stack
        myLifo.lib()
        lib_lcm = myLifo.node('LIBRARY')   

    elif mix_numbering_option == "number_mix_families_per_region":
        """
        Define material mixes based on the region numbering.
        Retrieve mix connectivity dict from geometry creation.
        """

        proc_name = "MIX_R_A.c2m"
        fill_lib_proc(proc_name, composition_option, mix_connectivity, resonance_correlation)
    
    return lib_lcm

def fill_generating_fuel_mix_definition_template(fuel_compositions_dict, correlation_option, generating_mixes_keys):
    """
    This function fills the mix definition template with data based on composition_option.
    Additionally if a compisition had already been defined, it will be "copied" through mix duplication.
    """
    corr_isotopes = ["U238", "Pu240", "Gd157"]
    self_shielded_isotopes = ["U234", "U235", "U238", "Pu239", "Pu240", "Pu241", "Pu242", "Gd154", "Gd155", "Gd156", "Gd157", "Gd158", "Gd160"] # This and inrs option could be passed as input for more flexibility
    definition_of_fuel_mixes = ""
    for key in generating_mixes_keys:
        """
        This loop will fill the definition of fuel mixes with the generating mixes keys.
        It will create a mix definition for each generating mix, using the fuel_compositions_dict.
        """
        _, composition_id, unique_id, ring_id = key.split("_")
        mix_compo = fuel_compositions_dict[str(composition_id)]
        definition_of_fuel_mixes += f"    MIX <<FMIX{composition_id}{unique_id}{ring_id}>> 750.0 \n" # Could pass temperature as input / depending on composition option ?
        for isotope, value in mix_compo.items():
            # Check if correlation is needed for specific isotopes
            if correlation_option == "CORR" and isotope in corr_isotopes:
                # Add correlation keyword for isotopes that require it
                definition_of_fuel_mixes += f"    {isotope} = {isotope} {value:.5E} CORR 1\n"
            elif isotope in self_shielded_isotopes:
                definition_of_fuel_mixes += f"    {isotope} = {isotope} {value:.5E} 1\n"
            else:
                definition_of_fuel_mixes += f"    {isotope} = {isotope} {value:.5E}\n"
            
    return definition_of_fuel_mixes

def duplicate_fuel_mixes(definition_of_fuel_mixes, generating_mixes_keys, mix_connectivity_dict):
    """
    Generating mixes have already been defined in the mix_composition_definition function.
    This function duplicates the fuel mixes based on the mix connectivity dictionary.

    """

def fill_lib_call(composition_option, mix_connectivity_dict, correlation_option):
    """
    This function fills the library call with the necessary data.
    Parameters:
    ----------
    composition_option : str
        The option for the mix composition definition (e.g., "AT10_void_0").
    Returns:
    ----------
    lib_call : str
        The filled library call string.
    """
    max_mix_number = max(mix_connectivity_dict.values()) if mix_connectivity_dict else 0
    generating_mixes_keys = findGeneratingMixesKeys(mix_connectivity_dict)
    print(f"Generating mixes found: {generating_mixes_keys}")

    fuel_compositions_dict, non_fuel_mixes = mix_composition_definition(composition_option)
    definition_of_fuel_mixes = fill_generating_fuel_mix_definition_template(fuel_compositions_dict, correlation_option, generating_mixes_keys)
    definition_of_fuel_mixes = duplicate_fuel_mixes(definition_of_fuel_mixes, generating_mixes_keys, mix_connectivity_dict)

    print(definition_of_fuel_mixes)
    lib_call = (
        "LIBRARY := LIB: ::\n"
        "EDIT 0\n"
        f"NMIX {max_mix_number}  ! MAXIMUM OF MATERIAL MIXTURES\n"
        "<<ssh_method>>\n" 
        "ANIS <<anis_level>>\n"
        "CTRA <<tran_correc>>\n"
        "DEPL LIB: DRAGON FIL: <<Library>>\n"
        "MIXS LIB: DRAGON FIL: <<Library>>\n"
        f"{definition_of_fuel_mixes}"
        f"{non_fuel_mixes}\n"
        ";\n"
    )
    print(lib_call)
    return

def findGeneratingMixesKeys(mix_connectivity_dict):
    """
    mix_connectivity_dict is a dictionary containing the mix connectivity information.
    the keys correspond to the unique mix identifiers, and the values are the mix numbers.
    mix indentifiers are of the form "FMIX_{composition_id}_{unique_id}_{ring_id}".
    generating mixes are extracted from the list of keys as the first mix identifier with unique composition_id
    """
    generating_mixes = []
    unique_enrich_ids = []
    for key in mix_connectivity_dict.keys():
        if key.startswith("FMIX_"):
            try:
                _, composition_id, unique_id, ring_id = key.split("_")
                if composition_id not in unique_enrich_ids:
                    unique_enrich_ids.append(composition_id)
                    generating_mixes.append(f"FMIX_{composition_id}_{unique_id}_{ring_id}")
            except ValueError:
                raise ValueError(f"Unexpected FMIX key format: {key}")
    # Remove duplicates while preserving order
            
    return generating_mixes

def fill_lib_proc(file_name, composition_option, mix_connectivity_dict, correlation_option):
    """
    This function fills the library procedure with the necessary data.
    Parameters:
    ----------
    composition_option : str
        The option for the mix composition definition (e.g., "AT10_void_0").
    mix_connectivity_dict : dict
        Dictionary containing the mix connectivity information.
    correlation_option : str
        Option to activate the resonance correlation model when using PT and RSE methods. 
        This adds the CORR keyword to U238, Pu240 and Gd157 isotopes in the LIB: call.
    """
    header = (  "* --------------------------------\n"
                "* Procedure generated by mix_handling.py\n"
                "* Author: R. Guasch\n"  
                "* --------------------------------\n"
                "*    INPUT & OUTPUT PARAMETERS\n"
                "* --------------------------------\n"
                "PARAMETER LIBRARY ::\n"
                "::: LINKED_LIST LIBRARY ; ;\n"
                "\n"
                "! Library name, composition option, ssh method\n"
                "STRING Library compos_opt ssh_method ;\n"
                ":: >>Library<<  >>compos_opt<< >>ssh_method<< ;\n" 
                "INTEGER anis_level ;\n"
                ":: >>anis_level<< ; ! Anisotropy level for sections angular depedance e.g. 1, 2, 3, 4, etc...\n"
                "STRING tran_correc ;\n"
                ":: >>tran_correc<< ; ! Transport correction option, e.g. , 'APOL', 'NONE'\n"
                "\n"
                "* -------------------------------\n"
                "*    STRUCTURES AND MODULES\n"
                "* -------------------------------\n"
                "MODULE  LIB: UTL: DELETE: END: ABORT: ;\n"
                "\n"
            )
    
    variable_declaration = cle2000_variable_declaration(mix_connectivity_dict)
    
    print(variable_declaration)
    lib_proc = (
        f"*PROCEDURE {file_name}\n"
        f"{header}\n"
        "* --------------------------------\n"
        "*    MIX NUMBERING DEFINITION\n"
        "* --------------------------------\n"
        f"{variable_declaration}\n"
        "* --------------------------------\n"
        "*    MIX COMPOSITION DEFINITION\n"
        "*    CALL TO LIB: MODULE\n"
        "* --------------------------------\n"
        f"{fill_lib_call(composition_option, mix_connectivity_dict, correlation_option)}\n"
        )
    return


def mix_composition_definition(composition_option, fuel_temperature=750.0, coolant_temperature=559.0):
    """
    This function retrieves the mix composition definition based on the specified option.
    Parameters:
    ----------
    composition_option : str
        The option for the mix composition definition (e.g., "AT10_void_0").
    Returns:
    ----------
    mix_composition_definition : dict of dict with str keys [mix][isotope] = value
    """
    N_H_H20 = 4.94546E-02  # Hydrogen in water
    N_O16_H2O = 2.47298E-02  # Oxygen in water
    
    fuel_composition = {
    "1": {"O16": 4.66705E-02, "U234": 5.15910E-06, "U235": 5.67035E-04, "U238": 2.27631E-02, "Pu239": 0.0, "Pu240": 0.0, "Pu241": 0.0, "Pu242": 0.0},
    "2": {"O16": 4.667480E-02, "U234": 7.039170E-06, "U235": 7.560370E-04, "U238": 2.257430E-02, "Pu239": 0.0, "Pu240": 0.0, "Pu241": 0.0, "Pu242": 0.0},
    "3": {"O16": 4.667960E-02, "U234": 9.163680E-06, "U235": 9.686590E-04, "U238": 2.236200E-02, "Pu239": 0.0, "Pu240": 0.0, "Pu241": 0.0, "Pu242": 0.0},
    "4": {"O16": 4.668150E-02, "U234": 9.991530E-06, "U235": 1.051340E-03, "U238": 2.227940E-02, "Pu239": 0.0, "Pu240": 0.0, "Pu241": 0.0, "Pu242": 0.0},
    "5": {"O16": 4.668280E-02, "U234": 1.058330E-05, "U235": 1.110400E-03, "U238": 2.222040E-02, "Pu239": 0.0, "Pu240": 0.0, "Pu241": 0.0, "Pu242": 0.0},
    "6": {"O16": 4.668410E-02, "U234": 1.117530E-05, "U235": 1.169460E-03, "U238": 2.216140E-02, "Pu239": 0.0, "Pu240": 0.0, "Pu241": 0.0, "Pu242": 0.0},
    "7": {"O16": 4.621410E-02, "U234": 9.451580E-06, "U235": 9.945290E-04, "U238": 2.107540E-02, "Pu239": 0.0, "Pu240": 0.0, "Pu241": 0.0,"Pu242": 0.0,
        "Gd154": 2.986510E-05, "Gd155": 2.027540E-04, "Gd156": 2.804310E-04, "Gd157": 2.143990E-04, "Gd158": 3.403000E-04, "Gd160": 2.994740E-04},
    "8":{"O16": 4.621230E-02, "U234": 8.668470E-06, "U235": 9.163120E-04, "U238": 2.115350E-02, "Pu239": 0.0, "Pu240": 0.0, "Pu241": 0.0,"Pu242": 0.0,
        "Gd154": 2.986510E-05, "Gd155": 2.027540E-04, "Gd156": 2.804310E-04, "Gd157": 2.143990E-04, "Gd158": 3.403000E-04, "Gd160": 2.994740E-04}
    }
    ### Look to pass these as input potentially, to prepare for multi-physics MULTICOMPO creation (temperatures).
    clad_temperature = 559.0  # Clad temperature in Kelvin
    box_temperature = 559.0  # Box temperature in Kelvin
    mode_temperature = 559.0  # Mode temperature in Kelvin
    fuel_temperature = 750.0  # Fuel temperature in Kelvin
    inrs_box = 2  # INRS option for box
    inrs_clad = 2  # INRS option for clad
    NOEV_definition = (
    f"MIX <<BOX>> {box_temperature} NOEV\n"
    "    O16  = O16   3.08132E-04\n"
    "    O17  = O17   1.17376E-07\n"
    "    Cr50  = Cr50   3.29613E-06\n"
    "    Cr52  = Cr52   6.35631E-05\n"
    "    Cr53  = Cr53   7.20746E-06\n"
     "   Cr54  = Cr54   1.79408E-06\n"
    "    Fe54  = Fe54   5.57337E-06\n"
    "    Fe56  = Fe56   8.74901E-05\n"
    "    Fe57  = Fe57   2.02053E-06\n"
    "    Fe58  = Fe58   2.68895E-07\n"
    "    Ni58  = Ni58   2.51627E-05\n"
    "    Ni60  = Ni60   9.69262E-06\n"
    "    Ni61  = Ni61   4.21331E-07\n" 
    "    Ni62  = Ni62   1.34344E-06\n"
    "    Ni64  = Ni64   3.42085E-07\n"
    f"    Zr90  = Zr90   2.18317E-02 {inrs_box}\n"
    f"    Zr91  = Zr91   4.76096E-03 {inrs_box}\n"
    f"    Zr92  = Zr92   7.27723E-03 {inrs_box}\n"
    f"    Zr94  = Zr94   7.37482E-03 {inrs_box}\n"
    f"    Zr96  = Zr96   1.18812E-03 {inrs_box}\n"
    "    Sn112  = Sn112   4.67343E-06\n"
    "    Sn114  = Sn114   3.17985E-06\n"
    "    Sn115  = Sn115   1.63812E-06\n"
    "    Sn116  = Sn116   7.00535E-05\n"
    "    Sn117  = Sn117   3.70020E-05\n"
    "    Sn118  = Sn118   1.16691E-04\n"
    "    Sn119  = Sn119   4.13864E-05\n"
    "    Sn120  = Sn120   1.56970E-04\n"
    "    Sn122  = Sn122   2.23073E-05\n"
    "    Sn124  = Sn124   2.78961E-05\n"
    "\n"
    f"MIX <<CLAD>> {clad_temperature} NOEV\n"
    "    O16  = O16   3.08132E-04\n"
    "    O17  = O17   1.17376E-07\n"
    "    Cr50  = Cr50   3.29613E-06\n"
    "    Cr52  = Cr52   6.35631E-05\n"
    "    Cr53  = Cr53   7.20746E-06\n"
     "   Cr54  = Cr54   1.79408E-06\n"
    "    Fe54  = Fe54   5.57337E-06\n"
    "    Fe56  = Fe56   8.74901E-05\n"
    "    Fe57  = Fe57   2.02053E-06\n"
    "    Fe58  = Fe58   2.68895E-07\n"
    "    Ni58  = Ni58   2.51627E-05\n"
    "    Ni60  = Ni60   9.69262E-06\n"
    "    Ni61  = Ni61   4.21331E-07\n" 
    "    Ni62  = Ni62   1.34344E-06\n"
    "    Ni64  = Ni64   3.42085E-07\n"
    f"    Zr90  = Zr90   2.18317E-02 {inrs_clad}\n"
    f"    Zr91  = Zr91   4.76096E-03 {inrs_clad}\n"
    f"    Zr92  = Zr92   7.27723E-03 {inrs_clad}\n"
    f"    Zr94  = Zr94   7.37482E-03 {inrs_clad}\n"
    f"    Zr96  = Zr96   1.18812E-03 {inrs_clad}\n"
    "    Sn112  = Sn112   4.67343E-06\n"
    "    Sn114  = Sn114   3.17985E-06\n"
    "    Sn115  = Sn115   1.63812E-06\n"
    "    Sn116  = Sn116   7.00535E-05\n"
    "    Sn117  = Sn117   3.70020E-05\n"
    "    Sn118  = Sn118   1.16691E-04\n"
    "    Sn119  = Sn119   4.13864E-05\n"
    "    Sn120  = Sn120   1.56970E-04\n"
    "    Sn122  = Sn122   2.23073E-05\n"
    "    Sn124  = Sn124   2.78961E-05\n"
    "\n"
    f"MIX <<GAP>> {fuel_temperature} NOEV\n"
    "    He4      = He4 1.50456E-04\n"
    "\n"
    f"MIX <<MODE>> {mode_temperature} NOEV"
    "    H1      = H1_H2O 4.94546E-02 ! Hydrogen in moderating water\n" 
    "    O16     = O16    2.47298E-02 ! Oxygen in moderating water\n"
    )
    if composition_option == "AT10_void_0":
        NOEV_definition += (
        "MIX <<COOL>> 559.0 NOEV\n"
        f"    H1      = H1_H2O {N_H_H20:.5E}  ! variable void coolant\n" 
        f"    O16     = O16    {N_O16_H2O:.5E}  ! variable void coolant\n"
        )
    elif composition_option == "AT10_void_40":
        NOEV_definition += (
        "MIX <<COOL>> 559.0 NOEV\n"
        f"    H1      = H1_H2O {N_H_H20*0.4:.5E}  ! variable void coolant\n" 
        f"    O16     = O16    {N_O16_H2O*0.4:.5E}  ! variable void coolant\n"
        )
    elif composition_option == "AT10_void_60":
        NOEV_definition += (
        "MIX <<COOL>> 559.0 NOEV\n"
        f"    H1      = H1_H2O {N_H_H20*0.6:.5E}  ! variable void coolant\n" 
        f"    O16     = O16    {N_O16_H2O*0.6:.5E}  ! variable void coolant\n"
        )
    elif composition_option == "AT10_void_80":
        NOEV_definition += (
        "MIX <<COOL>> 559.0 NOEV\n"
        f"    H1      = H1_H2O {N_H_H20*0.8:.5E}  ! variable void coolant\n" 
        f"    O16     = O16    {N_O16_H2O*0.8:.5E}  ! variable void coolant\n"
        )
    
    return fuel_composition, NOEV_definition


def cle2000_variable_declaration(mix_connectivity_dict):
    """
    Parameters:
    ----------
    mix_connectivity_dict : dict
        Dictionary containing the mix connectivity information.
    Returns:
    ----------
    variable_declaration : str
        String containing the CLE2000 variable declarations for the mix connectivity.
    """
    # Group FMIX variables by unique_id
    fmix_groups = defaultdict(list)
    other_vars = {}

    for key, value in mix_connectivity_dict.items():
        if key.startswith("FMIX"):
            try:
                _, composition_id, unique_id, ring_id = key.split("_")
                var_name = f"FMIX{composition_id}{unique_id}{ring_id}"
                fmix_groups[unique_id].append((var_name, value))
            except ValueError:
                raise ValueError(f"Unexpected FMIX key format: {key}")
        else:
            other_vars[key] = value

    variable_declaration = ""

    # Generate FMIX variable declarations
    for unique_id in sorted(fmix_groups, key=lambda x: int(x)):
        var_list = [var for var, _ in fmix_groups[unique_id]]
        value_list = [val for _, val in fmix_groups[unique_id]]
        var_line = " ".join(var_list)
        val_line = format_list_of_mix_numbers(value_list)
        variable_declaration += f"INTEGER {var_line} := {val_line} ;\n"

    # Add other variables
    for key, value in other_vars.items():
        variable_declaration += f"INTEGER {key} := {value} ;\n"

    return variable_declaration


def format_list_of_mix_numbers(list):
    """
    format a list of int to a string that can be used in a CLE2000 procedure call
    """

    if not list:
        return ""
    if len(list) == 1:
        return str(list[0])
    return " ".join(str(x) for x in list)

