## PyGan procedure to create ATRIUM-10 geometry, 
# Parameters which can be specified : 
#           - submeshing parameters for outer assembly water, bounding boxes, moderating water (in box),
#           - complexity in mix definiton fr self-shielding and rates condensation / homogenization,
import lifo 
import cle2000
import os
import numpy as np
import math

def createGeo(geo_name, split_water_in_moderator_box, split_moderator_box, split_coolant_around_moderator_box,
            split_intra_assembly_coolant, split_assembly_box, split_out_assembly_moderator,
            mix_numbering_option):
    """
    Create the geometry for the ATRIUM-10 BWR fuel assembly.
    Parameters
    ----------
    geo_name : str
        Name of the geometry.
    split_water_in_moderator_box : int
        Number of mesh splits for the water IN the moderator box (same value for X and Y splits).
    split_moderator_box : int
        Number of mesh splits for the moderator box (same value for X and Y splits).
        Effectively, this divides the zircalloy box bounding the central moderating water into sub regions.
    split_coolant_around_moderator_box : int
        Number of mesh splits for the coolant around the moderator box (same value for X and Y splits).
    split_intra_assembly_coolant : int
        Number of mesh splits for the intra-assembly water (coolant) (same value for X and Y splits).
        This represents the extra coolant between pincells lattice and assembly box.
    split_assembly_box : int
        Number of mesh splits for the assembly box (same value for X and Y splits).
        Effectively, this divides the zircalloy box bounding the fuel assembly into sub regions.
    split_out_assembly_moderator : int
        Number of mesh splits for the outer assembly water (same value for X and Y splits).
        This splits the moderator surrounding the assembly box.
    """
    
    if mix_numbering_option == "number_mix_families_per_enrichment":
        # Use the same mix number for the same composition : ie all fuel rods with the same composition will be self-shielded within the same mix
        # This is the simplest case, gives good agreement on integral quantities and total fission rates.
        # However due to the way MERG MIX works in the EDI: module, this effectively averages out rates per material MIX ie per composition in this case.
        ipLifo = lifo.new() 
        ipLifo.pushEmpty("GEOM", "LCM")
        ipLifo.pushEmpty("GEOMSSH", "LCM")
        ipLifo.push(geo_name)
        ipLifo.push(split_water_in_moderator_box)
        ipLifo.push(split_moderator_box)
        ipLifo.push(split_coolant_around_moderator_box)
        ipLifo.push(split_intra_assembly_coolant)
        ipLifo.push(split_assembly_box)
        ipLifo.push(split_out_assembly_moderator)

        # Create the geometry
        geomCreator = cle2000.new('GEO_A', ipLifo, 1)
        geomCreator.exec()
        # Recover the geometry nodes
        ipLifo.lib()
        pyGEOM = ipLifo.node("GEOM")
        pyGEOM_SSH = ipLifo.node("GEOMSSH") 

        return pyGEOM, pyGEOM_SSH, None

    elif mix_numbering_option == "number_mix_families_per_region":
        """
        In this option, the user selects to number mixes per region.
        The region numbering follows the X- to X+ and Y- to Y+ convention.
        Fuel mixes are numbered first, the remainder of the mix indices are attributed to box, cladding, moderator and coolant materials. 
        """
        lower_diag, UOX_mixes, Gd_mixes, water_box_sym_desc = getRegionToComposition(geo_name)
        # Renumber the lattice description and lower diagonal according to the numbering rule
        lattice_description, lower_diag = renumber_lattice_REGI(lower_diag, water_box_sym_desc )
        print(f"Renumbered lattice description for {geo_name}:")
        print(f"Lattice desc = {lattice_description}")
        print(f"lower diag = {lower_diag}")
        # Fill the template for CLE-2000 mix number definition
        if geo_name == "AT10_ASSBLY":
            remaining_mixes = ["GAP", "CLAD", "BOX", "MODE", "COOL"] # these are the remaining mixes, which are not defined in the lower diagonal list
            # Think of how to specify geometric dimensions for the geometry creation, for now assume template with geometrical data defined in header.
            #rfuel, rclad_in, rclad_out = 0.4435, 0.4520, 0.5140
            ## Check if the procedure already exists, if not create it
            proc_file_name = "GEO_R_A.c2m"
            if proc_file_name is not os.listdir:
                connectivity_dict = fill_geom_proc(proc_file_name, lattice_description, lower_diag, UOX_mixes, Gd_mixes, remaining_mixes)
                # change executable permission to the procedure file
                os.chmod(proc_file_name, 0o755)
            else:
                print("Geometry procedure already exists, skipping creation.")
        elif geo_name == "AT10_ASSBLY_CTRL":
            remaining_mixes = ["GAP", "CLAD", "BOX", "MODE", "COOL" "CTRL_ROD","CTRL_CRS"]
            proc_file_name = "GEO_R_CTR.c2m"
        ### Hardcoding geometrical info, should be recovered from config file and passed as argument when creating the procedure.
        moderator_box_inner_side = 3.34 # cm
        cell_pitch = 1.295 # cm
        split_water_in_moderator_box_center = int(math.ceil(split_water_in_moderator_box * (cell_pitch / moderator_box_inner_side))) # proportion of lines in center of water box.
        split_water_in_moderator_box_sides = int((split_water_in_moderator_box - split_water_in_moderator_box_center)/2) # keep the same number of splits as in previous defition, 
        # recompute in order to match the geometry definition in the procedure file.

        print(split_water_in_moderator_box_center)
        print(split_water_in_moderator_box_sides)
        # Execution of the generated geometry procedure
        ipLifo = lifo.new()
        ipLifo.pushEmpty("GEOM", "LCM")
        ipLifo.pushEmpty("GEOMSSH", "LCM")
        ipLifo.push(geo_name)
        ipLifo.push(split_water_in_moderator_box_center)
        ipLifo.push(split_water_in_moderator_box_sides)
        ipLifo.push(split_moderator_box)
        ipLifo.push(split_coolant_around_moderator_box)
        ipLifo.push(split_intra_assembly_coolant)
        ipLifo.push(split_assembly_box)
        ipLifo.push(split_out_assembly_moderator)

        # Create the geometry
        geomCreator = cle2000.new("GEO_R_A", ipLifo, 1)
        geomCreator.exec()
        # Recover the geometry nodes
        ipLifo.lib()
        pyGEOM = ipLifo.node("GEOM")
        pyGEOM_SSH = ipLifo.node("GEOMSSH")
        

        
        return pyGEOM, pyGEOM_SSH, connectivity_dict


def getRegionToComposition(name_geom):
    """
    This function returns a list of fuel enrichments, denominated by alias (ie C1/C2 ... to C8, water), for which the index in the list corresponds to the region number.
    Due to the diagonal symmetry : unique mixes only need to be defined for a given (Lower,Upper) pair and for mixes on the diagonal.

    Parameters
    ----------
    name_geom : str
        Name of the geometry. --> used to identify the lattice layout and define the region numbering. 
    numbering_rule : str
        Rule for numbering the mixes. Options are "number_mix_families_per_enrichment" or "number_mix_families_per_region".
        This defines how the mix numbers are assigned to the regions in the geometry.
        --> Sub optimal but best was for now : define the region numbering according to the lower diagonal of the lattice : renumber here depending on 
        numbering option so that the rest of procedure preparing functions can remain general
    """

    if name_geom == "AT10_ASSBLY": # ATRIUM-10 assembly, all rods, 10 by 10 lattice with 3x3 moderating water box
        lattice_description = ["C1", "C2", "C3", "C5", "C6", "C5", "C4", "C3", "C2", "C1", # BOTTOM OF LATTICE
                                "C2", "C4", "C7", "C6", "C7", "C6", "C6", "C7", "C4", "C2",
                                "C3", "C7", "C6", "C6", "C6", "C7", "C6", "C6", "C7", "C3",
                                "C5", "C6", "C6", "C6", "C6", "C6", "C7", "C6", "C5", "C4",
                                "C6", "C7", "C6", "C6", "W1", "WB", "W2", "C4", "C6" ,"C4",
                                "C5", "C6", "C7", "C6", "WL", "W0", "WR", "C3", "C7", "C4",
                                "C4", "C6", "C6", "C7", "W3", "WT", "W4", "C4", "C4", "C4",
                                "C3", "C7", "C6", "C6", "C4", "C3", "C4", "C4", "C8", "C3", 
                                "C2", "C4", "C7", "C5", "C6", "C7", "C4", "C8", "C4", "C2",
                                "C1", "C2", "C3", "C4", "C4 ","C4" ,"C4", "C3", "C2", "C1",] # TOP OF LATTICE
        
        lower_diag = ["C1", "C2", "C3", "C5", "C6", "C5", "C4", "C3", "C2", "C1", # BOTTOM OF LATTICE
                            "C4", "C7", "C6", "C7", "C6", "C6", "C7", "C4", "C2",
                                  "C6", "C6", "C6", "C7", "C6", "C6", "C7", "C3",
                                        "C6", "C6", "C6", "C7", "C6", "C5", "C4",
                                              "W1", "WB", "W2", "C4", "C6" ,"C4",
                                                    "W0", "WR", "C3", "C7", "C4",
                                                          "W4", "C4", "C4", "C4",
                                                                "C4", "C8", "C3", 
                                                                      "C4", "C2",
                                                                            "C1",] # TOP OF LATTICE
        
        UOX_mixes = ["C1", "C2", "C3", "C4", "C5", "C6"] # these require 4 mix numbers per unique position in the lower diagonal lattice
        Gd_mixes = ["C7", "C8"] # these require 6 mix numbers per unique position in the lower diagonal lattice

        water_box_desc = {"WB":"WT", "WR":"WL", "W2":"W3", "W0":"W0", "W1":"W1",  "W4":"W4"}

    return lower_diag, UOX_mixes, Gd_mixes, water_box_desc


def renumber_lattice_REGI(lower_diag, water_box_desc = {"WB":"WT", "WR":"WL", "W2":"W3", "W0":"W0", "W1":"W1",  "W4":"W4"}):
    """
    Renumber lattice according to the "number_mix_families_per_region" rule.
    
    Parameters:
    - lower_diag: list representing the lower triangular part (including diagonal) of a 10x10 lattice
    - water_box_desc: dictionnary associating the water box symmetry description
    
    Returns:
    - renumbered_lattice_description: full 10x10 list of numbers
    - renumbered_lower_diag: renumbered lower diagonal list
    """

    # Step 1: Assign unique integer labels to each composition
    renumbered_lower_diag = []

    for idx, desc in enumerate(lower_diag):
        if "W" in desc:
            renumbered_lower_diag.append(f"{desc}")
        else:
            renumbered_lower_diag.append(f"{desc}_{idx+1}")


    # Step 2: Unfold lower diagonal into full 10x10 symmetric matrix
    size = 10
    full_matrix = [["" for _ in range(size)] for _ in range(size)]
    index = 0
    ## Fill line by line starting from the bottom left corner : first element in the lower diagonal
    for i in range(size):  # extract line [row]
        for j in range(size): # iterate over columns in row
            if j >= i:  # Fill lower diagonal and diagonal elements
                if renumbered_lower_diag[index][0] == "W":
                    full_matrix[i][j] = renumbered_lower_diag[index]
                    full_matrix[j][i] = water_box_desc[renumbered_lower_diag[index]]  # Symmetric element
                else:
                    full_matrix[i][j] = renumbered_lower_diag[index]
                    full_matrix[j][i] = renumbered_lower_diag[index]  # Symmetric element
                index += 1

    # Step 3: Convert full_matrix to list of lists for output
    renumbered_lattice_description = full_matrix

    return renumbered_lattice_description, renumbered_lower_diag


def fill_template_mixes(lower_diag, UOX_mixes, Gd_mixes, remaining_mixes):
    """
    This function fills a template for CLE-2000 mix number definition.
    It is used to defined unique material mixtures for the ATRIUM-10 assembly.

    The naming scheme for CLE-2000 INTEGER variables is : INTEGER FMIX{index_in_lower_diag}{ring_number} 
    where index_in_lower_diag is the index of the mix in the lower diagonal list, associating a unique mixture to each position in the lower diagonal of the lattice,
    and ring_number is the ring number of the mix in the pincell (1 to 4 for UO2, 1 to 6 for UO2+Gd2O3).

    The remaining mixes are used for the box, cladding, moderator and coolant materials, which are not defined in the lower diagonal list.

    Parameters: 
    ----------
    lower_diag : list
        List of Cn,m materials, wher n defines the composition / initial fuel composition and m defined the individual numbering according to numbering rule.
        The position [index] in the list defines the unique mix according to its postion in the lower diagonal of the lattice.
    UOX_mixes : list
        List of UO2 mixes, which require to have 4 unique mix number per position.
    Gd_mixes : list
        List of Gd2O3 mixes, which require to have 6 unique mix number per position.
    remaining_mixes : list
        List of remaining mixes, which are not defined in the lower diagonal list.
        These are used for the box, cladding, moderator and coolant materials.

    Returns
    -------
    filled_mix_template : str
        String containing the filled CLE-2000 mix (INTEGER) definition template.
    connectivity_dict : dict
        Dictionary mapping the mix names FMIX_{composition}_{position}{ring_number} to their corresponding mix number.
    """
    filled_mix_template = ""
    connectivity_dict = {}
    last_n_mix = 0
    for index in range(len(lower_diag)):
        composition = lower_diag[index].split("_")[0] 
        composition_id = composition[-1]  # Extract the composition number from the mix description, ie C1, C2, C3, C4, C5, C6, C7, C8
        if composition in UOX_mixes:
            # For UO2 mixes, 4 rings -> 4 mix numbers
            if last_n_mix == 0:
                mix1, mix2, mix3, mix4 = 1, 2, 3, 4
            else:
                mix1, mix2, mix3, mix4 = last_mix + 1, last_mix + 2, last_mix + 3, last_mix + 4
            filled_mix_template += f"INTEGER FMIX{composition_id}{index+1}01 FMIX{composition_id}{index+1}02 FMIX{composition_id}{index+1}03 FMIX{composition_id}{index+1}04 := {mix1} {mix2} {mix3} {mix4} ;\n"
            connectivity_dict[f"FMIX_{composition_id}_{index+1}_01"] = mix1
            connectivity_dict[f"FMIX_{composition_id}_{index+1}_02"] = mix2
            connectivity_dict[f"FMIX_{composition_id}_{index+1}_03"] = mix3
            connectivity_dict[f"FMIX_{composition_id}_{index+1}_04"] = mix4
            last_n_mix = 4
            last_mix = mix4
        elif composition in Gd_mixes:
            if last_n_mix == 0:
                mix1, mix2, mix3, mix4, mix5, mix6 = 1, 2, 3, 4, 5, 6
            else:
                mix1, mix2, mix3, mix4, mix5, mix6 = last_mix + 1, last_mix + 2, last_mix + 3, last_mix + 4, last_mix + 5, last_mix + 6
            # For Gd2O3 mixes, 6 rings -> 6 mix numbers
            filled_mix_template += f"INTEGER FMIX{composition_id}{index+1}01 FMIX{composition_id}{index+1}02 FMIX{composition_id}{index+1}03 FMIX{composition_id}{index+1}04 FMIX{composition_id}{index+1}05 FMIX{composition_id}{index+1}06 := {mix1} {mix2} {mix3} {mix4} {mix5} {mix6} ;\n"
            connectivity_dict[f"FMIX_{composition_id}_{index+1}_01"] = mix1
            connectivity_dict[f"FMIX_{composition_id}_{index+1}_02"] = mix2
            connectivity_dict[f"FMIX_{composition_id}_{index+1}_03"] = mix3
            connectivity_dict[f"FMIX_{composition_id}_{index+1}_04"] = mix4
            connectivity_dict[f"FMIX_{composition_id}_{index+1}_05"] = mix5
            connectivity_dict[f"FMIX_{composition_id}_{index+1}_06"] = mix6
            last_n_mix = 6
            last_mix = mix6

    # Fill the remaining mixes remaining materials : specified as input to this function depending on the geometry
    for mix in remaining_mixes:
        if last_n_mix == 0:
            mix_number = 1
        else:
            mix_number = last_mix + 1
        filled_mix_template += f"INTEGER {mix} := {mix_number} ;\n"
        connectivity_dict[mix] = mix_number
        last_mix = mix_number

    return filled_mix_template.strip(), connectivity_dict


def format_CELL_definition(lattice_description):
    """
    Format the CELL definition for the ATRIUM-10 assembly geometry.
    This function takes the lattice description and formats it into a string that can be used in the CLE-2000 geometry definition.
    
    Parameters
    ----------
    lattice_description : list
        List of strings representing the lattice description, where each string is a mix description (e.g., "C1", "C2", etc.).
    
    Returns
    -------
    str
        Formatted string for the CELL definition.
    """
    cell_definition = ""
    for i in range(len(lattice_description)):
        cell_definition += " ".join(lattice_description[i])+ "\n" + "    "
    
    return cell_definition.strip()


def fill_CARCEL_template(index_in_lower_diag, lower_diag, UOX_mixes, Gd_mixes, sectorization, generating_cell):
    """
    Fill a generic CARCEL template for geometry definition according to a specific mix numbering.
    Assumes that sectorization is 4,n where n is the number of rings in the CARCEL. This ensure discretization of coolant only. 
    Parameters
    ----------
    index_in_lower_diag : int
        Index of the mix in the lower diagonal list, which defines the unique mix according to its position in the lower diagonal of the lattice.
    lower_diag : list
        List of Cn materials, where n defines the composition / initial fuel composition.
        The position [index] in the list defines the unique mix according to its postion in the lower diagonal of the lattice.
    UOX_mixes : list
        List of UO2 mixes, which require to have 4 unique mix numbers per position and 4 radii dividing the fuel, --> total .
    Gd_mixes : list
        List of Gd2O3 mixes, which require to have 6 unique mix numbers per position and 6 radii defining the cell.
    sectorization : boolean
        Flag to indicate whether the CARCEL should be sectorized or not. If True, sectorization is 4,n where n is the number of rings in the CARCEL.
        Sectorization must only be set to true in flux geometry, not in self-shielding geometry.
    generating_cell : str
        Name of the generating cell : first CARCEL to be defined, then used to defined subsequent cells by simply updating mixes.
        ie should be first UOX and fist Gd cells.
    """

    composition = lower_diag[index_in_lower_diag].split("_")[0]  # Extract the composition from the lower diagonal description, ie C1, C2, C3, C4, C5, C6, C7, C8
    composition_id = composition[-1]  # Extract the composition number from the mix description, ie C1, C2, C3, C4, C5, C6, C7, C8
    if composition in UOX_mixes:
        mix_numbering_list = [f"<<FMIX{composition_id}{index_in_lower_diag+1}0{i+1}>>" for i in range(4)]
        mix_numbering_str = " ".join(mix_numbering_list)
        if sectorization:
            if lower_diag[index_in_lower_diag] != generating_cell:
                CARCEL_def = (f"::: {lower_diag[index_in_lower_diag]}  :=  GEO: {generating_cell}\n" 
                            f"   MIX {mix_numbering_str} <<GAP>> <<CLAD>> \n"
                            "   <<COOL>> <<COOL>> <<COOL>> <<COOL>>\n"
                            "   <<COOL>> <<COOL>> <<COOL>> <<COOL>>\n"
                            "   <<COOL>> <<COOL>> <<COOL>> <<COOL>> ;\n"
                            )
            else:
                CARCEL_def = (f"::: {lower_diag[index_in_lower_diag]} :=  GEO: CARCEL 6\n"
                            "    SECT 4 6\n"
                            "    RADIUS 0.0 <<Rfuel1>> <<Rfuel2>> <<Rfuel3>> <<Rfuel4>> <<Rgap>> <<Rclad>>\n"
                            f"    MIX {mix_numbering_str}\n"
                            "   <<GAP>> <<CLAD>>\n" 
                            "   <<COOL>> <<COOL>> <<COOL>> <<COOL>>\n"
                            "   <<COOL>> <<COOL>> <<COOL>> <<COOL>>\n"
                            "   <<COOL>> <<COOL>> <<COOL>> <<COOL>>\n"
                            "   MESHX 0.0 <<pitch>>\n" 
                            "   MESHY 0.0 <<pitch>> ;\n"
                            "\n"
                            )
            return CARCEL_def.strip()
        else:
            if lower_diag[index_in_lower_diag] != generating_cell:
                CARCEL_def = (f"::: {lower_diag[index_in_lower_diag]}  :=  GEO: {generating_cell}\n" 
                            f"   MIX {mix_numbering_str}\n"
                            "   <<GAP>> <<CLAD>> <<COOL>> ;\n"
                            )
            else:
                CARCEL_def = (f"::: {lower_diag[index_in_lower_diag]}  :=  GEO: CARCEL 6\n"
                            "    RADIUS 0.0 <<Rfuel1>> <<Rfuel2>> <<Rfuel3>> <<Rfuel4>> <<Rgap>> <<Rclad>>\n"
                            f"    MIX {mix_numbering_str}\n"
                            "   <<GAP>> <<CLAD>> <<COOL>>\n" 
                            "   MESHX 0.0 <<pitch>>\n" 
                            "   MESHY 0.0 <<pitch>> ;\n"
                            "\n"
                            )
            return CARCEL_def.strip()
    elif composition in Gd_mixes:
        mix_numbering_list = [f"<<FMIX{composition_id}{index_in_lower_diag+1}0{i+1}>>" for i in range(6)]
        mix_numbering_str = " ".join(mix_numbering_list)
        if sectorization:
            if lower_diag[index_in_lower_diag] != generating_cell:
                CARCEL_def = (f"::: {lower_diag[index_in_lower_diag]}  :=  GEO: {generating_cell}\n" 
                            f"   MIX {mix_numbering_str}\n"
                            "  <<GAP>> <<CLAD>> \n"
                            "   <<COOL>> <<COOL>> <<COOL>> <<COOL>>\n"
                            "   <<COOL>> <<COOL>> <<COOL>> <<COOL>>\n"
                            "   <<COOL>> <<COOL>> <<COOL>> <<COOL>> ;\n"
                            )
            else:
                CARCEL_def = (f"::: {lower_diag[index_in_lower_diag]}  :=  GEO: CARCEL 8\n"
                            "    SECT 4 8\n"
                            "    RADIUS 0.0 <<RfuelGd1>> <<RfuelGd2>> <<RfuelGd3>> <<RfuelGd4>> <<RfuelGd5>> <<RfuelGd6>> <<Rgap>> <<Rclad>>\n"
                            f"    MIX {mix_numbering_str}\n"
                            "   <<GAP>> <<CLAD>>\n" 
                            "   <<COOL>> <<COOL>> <<COOL>> <<COOL>>\n"
                            "   <<COOL>> <<COOL>> <<COOL>> <<COOL>>\n"
                            "   <<COOL>> <<COOL>> <<COOL>> <<COOL>>\n"
                            "   MESHX 0.0 <<pitch>>\n" 
                            "   MESHY 0.0 <<pitch>> ;\n "
                            "\n"
                            ) 
            return CARCEL_def.strip()
        
        else:
            if lower_diag[index_in_lower_diag] != generating_cell:
                CARCEL_def = (f"::: {lower_diag[index_in_lower_diag]}  :=  GEO: {generating_cell}\n" 
                              f"MIX {mix_numbering_str} \n"
                              "<<GAP>> <<CLAD>> <<COOL>> ;\n"
                            )
            else:
                CARCEL_def = (f"::: {lower_diag[index_in_lower_diag]}  :=  GEO: CARCEL 8\n"
                            "    RADIUS 0.0 <<RfuelGd1>> <<RfuelGd2>> <<RfuelGd3>> <<RfuelGd4>> <<RfuelGd5>> <<RfuelGd6>> <<Rgap>> <<Rclad>>\n"
                            f"    MIX {mix_numbering_str}\n"
                            "   <<GAP>> <<CLAD>> <<COOL>>\n" 
                            "   MESHX 0.0 <<pitch>>\n" 
                            "   MESHY 0.0 <<pitch>> ;\n"
                            "\n"
                            ) 
            return CARCEL_def.strip()


def create_carcel_definitions(lower_diag, UOX_mixes, Gd_mixes, sectorization):
    """
    Create the CARCEL definitions for the ATRIUM-10 assembly geometry.
    """
    carcel_definitions = ""
    first_uox_cell = None
    first_gd_cell = None
    for index_in_lower_diag in range(len(lower_diag)):
        composition = lower_diag[index_in_lower_diag].split("_")[0] # Extract the composition from the lower diagonal description, ie C1, C2, C3, C4, C5, C6, C7, C8
        if composition in UOX_mixes and first_uox_cell is None:
            first_uox_cell = lower_diag[index_in_lower_diag]
        if composition in Gd_mixes and first_gd_cell is None:
            first_gd_cell = lower_diag[index_in_lower_diag]
        if composition not in UOX_mixes and composition not in Gd_mixes:
            print(f"Warning: Mix {composition} not found in UOX or Gd mixes, skipping CARCEL definition.")
            continue
        else:
            if composition in UOX_mixes:
                generating_cell = first_uox_cell
            elif composition in Gd_mixes:
                generating_cell = first_gd_cell
            carcel_definitions += fill_CARCEL_template(index_in_lower_diag, lower_diag, UOX_mixes, Gd_mixes, sectorization, generating_cell) + "\n"
    return carcel_definitions.strip()


def fill_self_shielding_geometry_template(lattice_description, lower_diag, UOX_mixes, Gd_mixes):
    """
    function used to fill the ATRIUM-10 geometry template with the appropriate fuel mix numbering
    geometrical dimensions are specificed in the geometrical_dimensions_template function.
    Need to handle control cross definition, for now stick to regular ATRIUM-10 geometry without control cross.
    """
    sectorization = False # Set to True if sectorization is required, False otherwise.

    ATRIUM_10_geo_template = (
        "GEOMSSH := GEO: :: CAR2D 3 3\n"  
        "   EDIT 1\n"
        "   X- DIAG X+ REFL\n"
        "   Y- REFL Y+ DIAG\n"
        "   CELL\n"
        "   BotCL BMidW BotCR\n"
        "           LAT RMidW\n"
        "               TopCR\n"
        "   MESHX 0.0 <<X1>> <<X4>> <<Pitch_A>>\n"
        "   MESHY 0.0 <<Y1>> <<Y4>> <<Pitch_A>>\n"

        "   ::: BotCL :=  GEO: CAR2D 3 3\n" # This depends on the presence of the control cross : have an option to handle both with and without control cross
        "       MESHX 0.0 <<W_gap>> <<x_step>> <<X1>>\n"
        "       MESHY 0.0 <<W_gap>> <<x_step>> <<Y1>>\n"
        "       MIX\n" 
        "           <<MODE>> <<MODE>> <<MODE>> \n" 
        "           <<MODE>>  <<BOX>>  <<BOX>> \n"
        "           <<MODE>>  <<BOX>> <<COOL>> ;\n"  

        "   ::: BotCR := GEO: CAR2D 3 3\n"
        "       MESHX 0.0 <<XtraCool>> <<y_step>> <<X1>>\n"
        "       MESHY 0.0 <<W_gap>> <<x_step>> <<X1>>\n"
        "       MIX\n"
        "           <<MODE>> <<MODE>> <<MODE>>\n"
        "           <<BOX>>  <<BOX>> <<MODE>> \n"
        "           <<COOL>> <<BOX>> <<MODE>> ;\n"
        
        "   ::: TopCR := GEO: CAR2D 3 3\n"
        "       MESHX 0.0 <<XtraCool>> <<y_step>> <<X1>>\n"
        "       MESHY 0.0 <<XtraCool>> <<y_step>> <<Y1>>\n"
        "       MIX\n"
        "           <<COOL>> <<BOX>> <<MODE>>  \n"
        "           <<BOX>>  <<BOX>> <<MODE>>  \n"
        "           <<MODE>> <<MODE>> <<MODE>> ;\n"

        "   ::: RMidW := GEO: CAR2D 3 1\n"
        "       MESHX 0.0 <<XtraCool>> <<y_step>> <<X1>>\n"
        "       MESHY 0.0 <<LLat>>\n"
        "       MIX <<COOL>> <<BOX>> <<MODE>> ;\n"

        "   ::: BMidW := GEO:  CAR2D 1 3\n" # This depends on the presence of the control cross : have an option to handle both
        "       MESHX 0.0 <<LLat>>\n"
        "       MESHY 0.0 <<W_gap>> <<x_step>> <<Y1>>\n"
        "       MIX <<MODE>> <<BOX>> <<COOL>> ;\n" # This depends on the presence of the control cross : have an option to handle both, in case of control cross need to call a cell definition for 3rd level
    # Lattice definition ; where the newly implemented material numbering comes into play, could consider breaking down the lattice definition into a separate function        
        "   ::: LAT :=  GEO: CAR2D 10 10\n"
        "       MESHX 0.0 <<Pitch_C>> <<2_Pitch_C>> <<3_Pitch_C>> <<4_Pitch_C>> <<5_Pitch_C>>\n"
        "                 <<6_Pitch_C>> <<7_Pitch_C>> <<8_Pitch_C>> <<9_Pitch_C>> <<LLat>>\n"
        "       MESHY 0.0 <<Pitch_C>> <<2_Pitch_C>> <<3_Pitch_C>> <<4_Pitch_C>> <<5_Pitch_C>>\n"
        "                 <<6_Pitch_C>> <<7_Pitch_C>> <<8_Pitch_C>> <<9_Pitch_C>> <<LLat>>\n"
        "       CELL\n"
        f"  {format_CELL_definition(lattice_description)}\n" # unfold geometry to fill the whole 10x10 square, creating a unique cell for each position in the lower diagonal, the same cell is used for the symmetric position in the upper diag.
        f"  {create_carcel_definitions(lower_diag, UOX_mixes, Gd_mixes, sectorization)}\n"
        "   ::: W1 := GEO: CAR2D 3 3\n"
        "       MESHX 0.0 <<XCHNL1>> <<XCHNL2>> <<Pitch_C>>\n" 
        "       MESHY 0.0 <<XCHNL1>> <<XCHNL2>> <<Pitch_C>>\n"
        "       MIX <<COOL>> <<COOL>> <<COOL>>\n"
        "           <<COOL>> <<BOX>> <<BOX>>\n"
        "           <<COOL>> <<BOX>> <<MODE>> ;\n"

        "   ::: WB := GEO: CAR2D 1 3\n"
        "       MESHX 0.0 <<Pitch_C>>\n"
        "       MESHY 0.0 <<XCHNL1>> <<XCHNL2>> <<Pitch_C>>\n"
        "       MIX <<COOL>> <<BOX>> <<MODE>> ;\n"
        
        "   ::: W2 := GEO: CAR2D 3 3\n"
        "       MESHX 0.0 <<X1sym>> <<X2sym>> <<Pitch_C>>\n"
        "       MESHY 0.0 <<XCHNL1>> <<XCHNL2>> <<Pitch_C>>\n"
        "       MIX <<COOL>> <<COOL>> <<COOL>>\n"
        "           <<BOX>>  <<BOX>>  <<COOL>>\n"
        "            <<MODE>> <<BOX>> <<COOL>> ;\n" 

        "   ::: WL := GEO: CAR2D 3 1\n" 
        "       MESHX 0.0 <<XCHNL1>> <<XCHNL2>> <<Pitch_C>>\n"
        "       MESHY 0.0 <<Pitch_C>>\n"
        "       MIX <<COOL>> <<BOX>> <<MODE>> ;\n"

        "   ::: W0 := GEO: CAR2D 1 1\n"
        "       MESHX 0.0 <<Pitch_C>>\n"
        "       MESHY 0.0 <<Pitch_C>>\n"
        "        MIX <<MODE>> ;\n"

        "   ::: WR := GEO: CAR2D 3 1\n"
        "       MESHX 0.0 <<X1sym>> <<X2sym>> <<Pitch_C>>\n"
        "       MESHY 0.0 <<Pitch_C>>\n"
        "       MIX <<MODE>> <<BOX>> <<COOL>> ;\n"

        "   ::: W3 := GEO: CAR2D 3 3\n" 
        "       MESHX 0.0 <<XCHNL1>> <<XCHNL2>> <<Pitch_C>>\n"
        "       MESHY 0.0 <<X1sym>> <<X2sym>> <<Pitch_C>>\n"
        "       MIX <<COOL>> <<BOX>> <<MODE>>\n"
        "           <<COOL>> <<BOX>> <<BOX>>\n"
        "           <<COOL>> <<COOL>> <<COOL>> ;\n"
            
        "   ::: WT := GEO: CAR2D 1 3\n" 
        "       MESHX 0.0 <<Pitch_C>>\n"
        "       MESHY 0.0 <<X1sym>> <<X2sym>> <<Pitch_C>>\n"
        "       MIX <<MODE>> <<BOX>> <<COOL>> ;\n"
        
        "   ::: W4 := GEO: CAR2D 3 3\n"
        "        MESHX 0.0 <<X1sym>> <<X2sym>> <<Pitch_C>>\n"
        "        MESHY 0.0 <<X1sym>> <<X2sym>> <<Pitch_C>>\n"
        "        MIX <<MODE>> <<BOX>> <<COOL>>\n"
        "            <<BOX>> <<BOX>> <<COOL>>\n"
        "            <<COOL>> <<COOL>> <<COOL>> ;\n"
        "   ;\n"
        ";\n"
    )
    return ATRIUM_10_geo_template.strip()


def fill_flux_geometry_template(lattice_description, lower_diag, UOX_mixes, Gd_mixes):
    """
    function used to fill the ATRIUM-10 geometry template with the appropriate fuel mix numbering
    geometrical dimensions are specificed in the geometrical_dimensions_template function.
    Need to handle control cross definition, for now stick to regular ATRIUM-10 geometry without control cross.
    """
    sectorization = True # Set to True if sectorization is required, False otherwise.

    ATRIUM_10_geo_template = (
        "GEOM := GEO: :: CAR2D 3 3\n"  
        "   EDIT 1\n"
        "   X- DIAG X+ REFL\n"
        "   Y- REFL Y+ DIAG\n"
        "   CELL\n"
        "   BotCL BMidW BotCR\n"
        "           LAT RMidW\n"
        "               TopCR\n"
        "   MESHX 0.0 <<X1>> <<X4>> <<Pitch_A>>\n"
        "   MESHY 0.0 <<Y1>> <<Y4>> <<Pitch_A>>\n"

        "   ::: BotCL :=  GEO: CAR2D 3 3\n" # This depends on the presence of the control cross : have an option to handle both with and without control cross
        "       MESHX 0.0 <<W_gap>> <<x_step>> <<X1>>\n"
        "       MESHY 0.0 <<W_gap>> <<x_step>> <<Y1>>\n"
        "       SPLITX <<sp_out_mod>> <<sp_as_box>> <<sp_in_as_co>>\n"
        "       SPLITY <<sp_out_mod>> <<sp_as_box>> <<sp_in_as_co>>\n"
        "       MIX\n" 
        "           <<MODE>> <<MODE>> <<MODE>> \n" 
        "           <<MODE>>  <<BOX>>  <<BOX>> \n"
        "           <<MODE>>  <<BOX>> <<COOL>> ;\n"  

        "   ::: BotCR := GEO: CAR2D 3 3\n"
        "       MESHX 0.0 <<XtraCool>> <<y_step>> <<X1>>\n"
        "       MESHY 0.0 <<W_gap>> <<x_step>> <<X1>>\n"
        "       SPLITX <<sp_in_as_co>> <<sp_as_box>> <<sp_out_mod>>\n"
        "       SPLITY <<sp_out_mod>> <<sp_as_box>> <<sp_in_as_co>>\n"
        "       MIX\n"
        "           <<MODE>> <<MODE>> <<MODE>>\n"
        "           <<BOX>>  <<BOX>> <<MODE>> \n"
        "           <<COOL>> <<BOX>> <<MODE>> ;\n"
        
        "   ::: TopCR := GEO: CAR2D 3 3\n"
        "       MESHX 0.0 <<XtraCool>> <<y_step>> <<X1>>\n"
        "       MESHY 0.0 <<XtraCool>> <<y_step>> <<Y1>>\n"
        "       SPLITX <<sp_in_as_co>> <<sp_as_box>> <<sp_out_mod>>\n"
        "       SPLITY <<sp_in_as_co>> <<sp_as_box>> <<sp_out_mod>>\n"
        "       MIX\n"
        "           <<COOL>> <<BOX>> <<MODE>>\n"
        "           <<BOX>>  <<BOX>> <<MODE>> \n"
        "           <<MODE>> <<MODE>> <<MODE>> ;\n"

        "   ::: RMidW := GEO: CAR2D 3 1\n"
        "       MESHX 0.0 <<XtraCool>> <<y_step>> <<X1>>\n"
        "       MESHY 0.0 <<LLat>>\n"
        "       SPLITX <<sp_in_as_co>> <<sp_as_box>> <<sp_out_mod>>\n"
        "       SPLITY <<sp_out_mod>>\n"
        "       MIX <<COOL>> <<BOX>> <<MODE>> ;\n"

        "   ::: BMidW := GEO:  CAR2D 1 3\n" # This depends on the presence of the control cross : have an option to handle both
        "       MESHX 0.0 <<LLat>>\n"
        "       MESHY 0.0 <<W_gap>> <<x_step>> <<Y1>>\n"
        "       SPLITX <<sp_out_mod>>\n"
        "       SPLITY <<sp_out_mod>> <<sp_as_box>> <<sp_in_as_co>>\n"
        "       MIX <<MODE>> <<BOX>> <<COOL>> ;\n" # This depends on the presence of the control cross : have an option to handle both, in case of control cross need to call a cell definition for 3rd level
    # Lattice definition ; where the newly implemented material numbering comes into play, could consider breaking down the lattice definition into a separate function        
        "   ::: LAT :=  GEO: CAR2D 10 10\n"
        "       MESHX 0.0 <<Pitch_C>> <<2_Pitch_C>> <<3_Pitch_C>> <<4_Pitch_C>> <<5_Pitch_C>>\n"
        "                 <<6_Pitch_C>> <<7_Pitch_C>> <<8_Pitch_C>> <<9_Pitch_C>> <<LLat>>\n"
        "       MESHY 0.0 <<Pitch_C>> <<2_Pitch_C>> <<3_Pitch_C>> <<4_Pitch_C>> <<5_Pitch_C>>\n"
        "                 <<6_Pitch_C>> <<7_Pitch_C>> <<8_Pitch_C>> <<9_Pitch_C>> <<LLat>>\n"
        "       CELL\n"
        f"  {format_CELL_definition(lattice_description)}\n" # unfold geometry to fill the whole 10x10 square, creating a unique cell for each position in the lower diagonal, the same cell is used for the symmetric position in the upper diag.
        f"  {create_carcel_definitions(lower_diag, UOX_mixes, Gd_mixes, sectorization)}\n"
        "   ::: W1 := GEO: CAR2D 3 3\n"
        "       MESHX 0.0 <<XCHNL1>> <<XCHNL2>> <<Pitch_C>>\n" 
        "       MESHY 0.0 <<XCHNL1>> <<XCHNL2>> <<Pitch_C>>\n"
        "       SPLITX <<sp_co_aroud>> <<sp_md_box>> <<sp_mod_s>>\n"
        "       SPLITY <<sp_co_aroud>> <<sp_md_box>> <<sp_mod_s>>\n"
        "       MIX <<COOL>> <<COOL>> <<COOL>>\n"
        "           <<COOL>> <<BOX>> <<BOX>>\n"
        "           <<COOL>> <<BOX>> <<MODE>> ;\n"

        "   ::: WB := GEO: CAR2D 1 3\n"
        "       MESHX 0.0 <<Pitch_C>>\n"
        "       MESHY 0.0 <<XCHNL1>> <<XCHNL2>> <<Pitch_C>>\n"
        "       SPLITX <<sp_mod_s>> \n"
        "       SPLITY <<sp_co_aroud>> <<sp_md_box>> <<sp_mod_s>>\n"
        "       MIX <<COOL>> <<BOX>> <<MODE>> ;\n"
        
        "   ::: W2 := GEO: CAR2D 3 3\n"
        "       MESHX 0.0 <<X1sym>> <<X2sym>> <<Pitch_C>>\n"
        "       MESHY 0.0 <<XCHNL1>> <<XCHNL2>> <<Pitch_C>>\n"
        "       SPLITX <<sp_mod_s>> <<sp_md_box>> <<sp_co_aroud>>\n"
        "       SPLITY <<sp_co_aroud>> <<sp_md_box>> <<sp_mod_s>>\n"
        "       MIX <<COOL>> <<COOL>> <<COOL>>\n"
        "           <<BOX>>  <<BOX>>  <<COOL>>\n"
        "            <<MODE>> <<BOX>> <<COOL>> ;\n" 

        "   ::: WL := GEO: CAR2D 3 1\n" 
        "       MESHX 0.0 <<XCHNL1>> <<XCHNL2>> <<Pitch_C>>\n"
        "       MESHY 0.0 <<Pitch_C>>\n"
        "       SPLITX <<sp_co_aroud>> <<sp_md_box>> <<sp_mod_s>>\n"
        "       SPLITY <<sp_mod_s>>\n"
        "       MIX <<COOL>> <<BOX>> <<MODE>> ;\n"

        "   ::: W0 := GEO: CAR2D 1 1\n"
        "       MESHX 0.0 <<Pitch_C>>\n"
        "       MESHY 0.0 <<Pitch_C>>\n"
        "       SPLITX <<sp_mod_c>>\n"
        "       SPLITY <<sp_mod_c>>\n"
        "       MIX <<MODE>> ;\n"

        "   ::: WR := GEO: CAR2D 3 1\n"
        "       MESHX 0.0 <<X1sym>> <<X2sym>> <<Pitch_C>>\n"
        "       MESHY 0.0 <<Pitch_C>>\n"
        "       SPLITX <<sp_mod_s>> <<sp_md_box>> <<sp_co_aroud>>\n"
        "       SPLITY <<sp_mod_s>>\n"
        "       MIX <<MODE>> <<BOX>> <<COOL>> ;\n"

        "   ::: W3 := GEO: CAR2D 3 3\n" 
        "       MESHX 0.0 <<XCHNL1>> <<XCHNL2>> <<Pitch_C>>\n"
        "       MESHY 0.0 <<X1sym>> <<X2sym>> <<Pitch_C>>\n"
        "       SPLITX <<sp_co_aroud>> <<sp_md_box>> <<sp_mod_s>>\n"
        "       SPLITY <<sp_mod_s>> <<sp_md_box>> <<sp_co_aroud>>\n"
        "       MIX <<COOL>> <<BOX>> <<MODE>>\n"
        "           <<COOL>> <<BOX>> <<BOX>>\n"
        "           <<COOL>> <<COOL>> <<COOL>> ;\n"
            
        "   ::: WT := GEO: CAR2D 1 3\n" 
        "       MESHX 0.0 <<Pitch_C>>\n"
        "       MESHY 0.0 <<X1sym>> <<X2sym>> <<Pitch_C>>\n"
        "       SPLITX <<sp_mod_c>>\n"
        "       SPLITY <<sp_mod_c>> <<sp_md_box>> <<sp_co_aroud>>\n"
        "       MIX <<MODE>> <<BOX>> <<COOL>> ;\n"
        
        "   ::: W4 := GEO: CAR2D 3 3\n"
        "        MESHX 0.0 <<X1sym>> <<X2sym>> <<Pitch_C>>\n"
        "        MESHY 0.0 <<X1sym>> <<X2sym>> <<Pitch_C>>\n"
        "        SPLITX <<sp_mod_s>> <<sp_md_box>> <<sp_co_aroud>>\n"
        "        SPLITY <<sp_mod_s>> <<sp_md_box>> <<sp_co_aroud>>\n"
        "        MIX <<MODE>> <<BOX>> <<COOL>>\n"
        "            <<BOX>> <<BOX>> <<COOL>>\n"
        "            <<COOL>> <<COOL>> <<COOL>> ;\n"
        "   ;\n"
        ";\n"
    )
    return ATRIUM_10_geo_template.strip()


def fill_geom_proc(file_name, lattice_description, lower_diag, UOX_mixes, Gd_mixes, remaining_mixes):

    """
    Function to fill/create the CLE2000 procedure file with the ATRIUM-10 geometry definition.
    This function is used to create the geometry for the ATRIUM-10 assembly, which is a 10x10 lattice with a 3x3 moderating water box.
    
    Parameters
    ----------
    file_name : str
        Name of the file to write the geometry definition to.
    
    Returns
    -------
    connectivity_dict : dict
    """
    header = (  "* --------------------------------\n"
                "*    GEOMETRY DEFINITION PROCEDURE\n"
                "*    ATRIUM-10 ASSEMBLY\n"
                "* Procedure generated by create_geo.py\n"
                "* Author: R. Guasch\n"
                "* --------------------------------\n"
                "*    INPUT & OUTPUT PARAMETERS\n"
                "* --------------------------------\n"
                "PARAMETER GEOM GEOMSSH ::\n"
                "::: LINKED_LIST GEOM ;\n"
                "::: LINKED_LIST GEOMSSH ; ;\n"
                "\n"
                "STRING name_geom ;\n"
                ":: >>name_geom<< ;\n"
                "INTEGER sp_mod_c sp_mod_s sp_md_box sp_co_aroud sp_in_as_co sp_as_box sp_out_mod ;\n"
                ":: >>sp_mod_c<< >>sp_mod_s<< >>sp_md_box<< >>sp_co_aroud<< >>sp_in_as_co<< >>sp_as_box<< >>sp_out_mod<< ;\n"
                "\n"
                "* -------------------------------\n"
                "*    STRUCTURES AND MODULES\n"
                "* -------------------------------\n"
                "\n"
                "MODULE  GEO: END: ;\n"
                "\n"
            )
    

    # Fill the geometrical dimensions template
    geometrical_dimensions = geometrical_dimensions_template()
    # Mix numbers definition
    mix_numbers_definition, connectivity_dict = fill_template_mixes(lower_diag, UOX_mixes, Gd_mixes, remaining_mixes)

    
    # consrtuct GEO_REGI_A.c2m procedure
    
    geo_proc = (
        f"*PROCEDURE {file_name}\n"
        f"{header}\n"
        "* --------------------------------\n"
        "*    DIMENSIONS DEFINITION\n"
        "* --------------------------------\n"
        f"{geometrical_dimensions}\n"
        "*Mix numbering definition\n"
        f"{mix_numbers_definition}\n"
        "* -----------------------------------------\n"
        "*    SELF-SHIELDING GEOMETRY DEFINITION\n"
        "* -----------------------------------------\n"

        f"{fill_self_shielding_geometry_template(lattice_description, lower_diag, UOX_mixes, Gd_mixes)}\n"
        
        "* -----------------------------------------\n"
        "*         FLUX GEOMETRY DEFINITION\n"
        "* -----------------------------------------\n"
        f"{fill_flux_geometry_template(lattice_description, lower_diag, UOX_mixes, Gd_mixes)}\n"

        "* -----------------------------------------\n"
        "*         END OF GEOMETRY DEFINITION\n"
        "* -----------------------------------------\n"
        "END: ;\n"
        "QUIT ."
    )

    # Write the procedure file
    with open(file_name, "w") as file:
        file.write(geo_proc)
    file.close()
    print(f"Geometry procedure file {file_name} created successfully.")
    return connectivity_dict


def geometrical_dimensions_template():
    """
    To be generalized to be able to define more general geometries,
    for now "hard code" header to CLE-2000 procedure definition.
    
    Parameters
    ----------
    None

    Returns
    -------
        geometrical_dimensions_def : str
            String containing the geometrical dimensions definition / CLE-2000 variable definition for geometrical data.
    """

    geometrical_dimensions_def = """
!Assembly scale data

REAL Pitch_A W_gap := 15.24 0.75 ; ! Assembly pitch and water gap thickness
REAL Box_out Box_in Box_thi := 13.74 13.4 0.17 ; ! Outer box outer, inner sides and thickness
REAL Chan_out Chan_in Chan_thi := 3.5 3.34 0.08 ; ! Channel box outer, inner sides and thickness


! UOX Pincell scale data

REAL Pitch_C := 1.29500 ;
REAL Rgap Rclad pitch := 0.4520 0.5140 1.29500 ;
REAL Rfuel1 Rfuel2 Rfuel3 Rfuel4 ;
EVALUATE Rfuel4 := 0.4435 ;
EVALUATE Rfuel1 := 0.313602 ;
EVALUATE Rfuel2 := 0.396678 ; 
EVALUATE Rfuel3 := 0.43227 ;

REAL RfuelGd1 RfuelGd2 RfuelGd3 RfuelGd4 RfuelGd5 RfuelGd6 ;
EVALUATE RfuelGd1 RfuelGd2 RfuelGd3 RfuelGd4 RfuelGd5 RfuelGd6 
        := 0.19834 0.28049 0.34353 0.39668 0.43227 0.4435 ;

! Evaluating dependencies

! Multiples of cell pitch : used to define 2nd level geometry mesh 
REAL 2_Pitch_C 3_Pitch_C 4_Pitch_C 5_Pitch_C 6_Pitch_C 7_Pitch_C 8_Pitch_C 9_Pitch_C ;
EVALUATE 2_Pitch_C := Pitch_C 2.0 * ;
EVALUATE 3_Pitch_C := Pitch_C 3.0 * ;
EVALUATE 4_Pitch_C := Pitch_C 4.0 * ;
EVALUATE 5_Pitch_C := Pitch_C 5.0 * ;
EVALUATE 6_Pitch_C := Pitch_C 6.0 * ;
EVALUATE 7_Pitch_C := Pitch_C 7.0 * ;
EVALUATE 8_Pitch_C := Pitch_C 8.0 * ;
EVALUATE 9_Pitch_C := Pitch_C 9.0 * ;
REAL LLat := Pitch_C 10.0 * ;


! X1 = Water gap + outer box thickness +extra water : in box water not a multiple of Pitch_C
REAL XtraCool ;
EVALUATE XtraCool := Box_in 10.0 Pitch_C * - 2.0 / ;

!Meshing points for macro geometry :
REAL X1 X2 X3 X4  ;
! X2 = X1 + 4*Pitch_C
! X3 = X2 + 3*Pitch_C
! X4 = X3 + 3*Pitch_C
EVALUATE X1 := W_gap Box_thi + XtraCool + ;
EVALUATE X2 := 4.0 Pitch_C * X1 + ;
EVALUATE X3 := 3.0 Pitch_C * X2 + ;
EVALUATE X4 := 3.0 Pitch_C * X3 + ;  

REAL Y1 Y2 Y3 Y4 ;
EVALUATE Y1 := X1 ;
EVALUATE Y2 := 4.0 Pitch_C * Y1 + ;
EVALUATE Y3 := 3.0 Pitch_C * Y2 + ;
EVALUATE Y4 := 3.0 Pitch_C * Y3 + ;

REAL L1 L2 L3 ;
EVALUATE L1 := X2 X1 - ;
EVALUATE L2 := X3 X2 - ;
EVALUATE L3 := X4 X3 - ;

REAL H1 H2 H3 ;
EVALUATE H1 := Y2 Y1 - ;
EVALUATE H2 := Y3 Y2 - ;
EVALUATE H3 := Y4 Y3 - ;

REAL x_step y_step ; ! allows for submeshing of in-assembly box water
! in-assembly box water included in geom level with outer assembly water and box
EVALUATE x_step := W_gap Box_thi + ;
EVALUATE y_step := XtraCool Box_thi + ;

! inner CHaNneL box related parameters, used to define 2nd level CHNL Geometry
! CHNL = channel box + surroudning water to fit with pin CARCEL edges 
REAL XCHNL1 XCHNL2 XCHNL3 XCHNL4 XCHNL5 ;
EVALUATE XCHNL1 := L2 Chan_out - 2.0 / ;
EVALUATE XCHNL2 := XCHNL1 Chan_thi + ;
EVALUATE XCHNL3 := XCHNL2 Chan_in + ;
EVALUATE XCHNL4 := XCHNL3 Chan_thi + ;
EVALUATE XCHNL5 := L2 ;

REAL X1sym := Pitch_C XCHNL2 - ;
REAL X2sym := Pitch_C XCHNL1 - ;
    """
    return geometrical_dimensions_def.strip()