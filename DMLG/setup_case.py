## Main setup_case.py file to be used to set up a neutronics BWR case
# Idea is to set up a case for given a geometry and material composition
# Options could involve setting up a Serpent2 input file, a DRAGON5 case (collection of c2m files) 
# Options could allow for different types of calculations.

# General philosophy : yaml files for configuration, python for logic and setting up the case.
import yaml
import numpy as np
from pathlib import Path
from DMLG_composition_handling.material_mixture import Material_Mixture
from template.Serpent2.geometry_definitions import Lattice_Definition, Pin_Universe_Definition
#from DMLG_geometry_handling.lattice import BWR_lattice



class DMLG_case:
    def __init__(self, name:str, output_codes:list, nuclear_data_evaluation:str = "endfb8r1"):
        """
        name (str) : used to identify the input and output files associated with this case
        output_codes (list) : list of codes used for the output files, e.g. ['Serpent2', 'DRAGON5']
        """
        self.name = name
        self.output_codes = output_codes
        self.nuclear_data_evaluation = nuclear_data_evaluation
        self.case_dir = Path(f"./output_cases/{self.name}")
        self.case_dir.mkdir(parents=True, exist_ok=True)
        self.materials = {}
        self.geometry = None
        self.settings = {}
        self.load_materials(f"case_input/{self.name}_mix_def.yml")

    def load_materials(self, materials_file: str):
        """
        Load materials from a YAML file.
        materials_file (str): Path to the YAML file containing material definitions.
        """
        with open(materials_file, 'r') as file:
            materials_data = yaml.safe_load(file)
        compositions_info = materials_data["MIX_COMPOSITIONS"]
        for material in compositions_info:
            material_name = material['name']
            isotopic_composition = material['isotopic_composition']
            temperature_points = material['temperature_points'] if 'temperature_points' in material else None
            is_depletable = material['depletable'] if 'depletable' in material else False
            print(f"Loading material: {material_name}")
            for temp in temperature_points:
                self.materials[f"{material_name} {temp}"] = Material_Mixture(material_name, isotopic_composition, temp, self.nuclear_data_evaluation, is_depletable)

    def print_materials_to_D5(self):
        """
        Print the material definitions to a DRAGON5 compatible format.
        This will create a .c2m file with a LIBRARY definition throught the LIB: module of DRAGON5.
        """
        c2m_file_path = self.case_dir / f"MIX_DEF.c2m"
        LIB_definition = ['LIBRARY := LIB: ::\n']    
                          
        if self.settings["DRAGON5"]["Self-Shiedling Method"]:
            LIB_definition.append(f'{self.settings["DRAGON5"]["Self-Shiedling Method"]}\n')
        if self.settings["DRAGON5"]["Order Anisotropic Scattering"]:
            LIB_definition.append(f'ANIS {self.settings["DRAGON5"]["Order Anisotropic Scattering"]+1}\n')
        if self.settings["DRAGON5"]["Number of mixtures"]:
            LIB_definition.append(f'NMIX {self.settings["DRAGON5"]["Number of mixtures"]}\n')
        if self.settings["DRAGON5"]["DRAGLIB"]:
            LIB_definition.append(f'MIXS LIB: DRAGON FIL: {self.settings["DRAGON5"]["DRAGLIB"]}\n')
        if self.settings["DRAGON5"]["Depletion"] and self.settings["DRAGON5"]["DRAGLIB"]:
            LIB_definition.append(f'DEPL LIB: DRAGON FIL: {self.settings["DRAGON5"]["DRAGLIB"]}\n')

        for material_name, material in self.materials.items():
            if material.isdepletable:
                LIB_definition.append(f'MIX {material_name.split(" ")[0]} {material_name.split(" ")[-1]}\n')
            else:
                LIB_definition.append(f'MIX {material_name.split(" ")[0]} {material_name.split(" ")[-1]} NOEV\n')
            LIB_definition.append(material.print_to_D5_format())

        with open(c2m_file_path, 'w') as c2m_file:
            for line in LIB_definition:
                print(line)
                c2m_file.write(line)

        return

    def set_lattice_description(self, lattice_desc:list, pitch:float):
        """setting the lattice desciption associated to treated case. 
        assume square cartesian lattice without symmeties (for now).
        Args:
            lattice_desc (list): list of pincell identifiers making up the lattice in x-increasing y-increasing order.
                    --> Only cells interior to outer box should be included in this description.
            
        """
        ## assuming square lattice : 
        nx = int(np.sqrt(len(lattice_desc)))
        ny = nx
        


        if "Serpent2" in self.output_codes:
            # assume lattice is centered at (0.0, 0.0)
            self.lattice_description = Lattice_Definition("BWR_lattice", pitch, 0.0, 0.0, nx, ny, pincell_universes)
            
    
    def set_pincell_id_to_mat_id(self, pincell_to_material_association):
        """ 
        set a dictionnary associating pincell id to material id.
        each pair of ends of getting a unique id 
        """
        self.pincell_nb_to_mat_name = pincell_to_material_association

if __name__ == "__main__":
    # Test usage on OECD PHASE IIIB BWR assembly benchmark case :
    case = DMLG_case(name="OECD_NEA_PHASE_IIIB", output_codes=["Serpent2", "DRAGON5"], nuclear_data_evaluation="endfb8r1")
    case.settings = {"DRAGON5": {"DRAGLIB": "endfb8r1_295",
                                "Self-Shiedling Method": "PT",
                                "Order Anisotropic Scattering": 3,
                                "Number of mixtures": 10,
                                "Depletion": True}}
    case.print_materials_to_D5()

    # Usage on ATRIUM-10 BWR fuel bundle.

    ATRIUM10 = DMLG_case(name="ATRIUM10", output_codes = ["Serpent2", "DRAGON5"], nuclear_data_evaluation="endfb8r1")
    ATRIUM10.settings = {"DRAGON5": {"DRAGLIB": "endfb8r1_295",
                                "Self-Shiedling Method": "PT",
                                "Order Anisotropic Scattering": 3,
                                "Number of mixtures": 10,
                                "Depletion": True}}
    ATRIUM10.set_pincell_id_to_mat_id({1:"24UOX", 2:"32UOX", 3:"42UOX", 4:"45UOX", 5:"48UOX", 6:"50UOX", 7:"45Gd", 8:"42Gd"})
    ATRIUM10.set_lattice_description([  1, 2, 3, 4, 4, 4, 4, 3, 2, 1,
                                        2, 4, 7, 5, 6, 7, 4, 8, 4, 2,
                                        3, 7, 6, 6, 4, 3, 4, 4, 8, 3,
                                        4, 6, 6, 7, 9, 9, 9, 4, 4, 4,
                                        5, 6, 7, 6, 9, 9, 9, 3, 7, 4,
                                        6, 7, 6, 6, 9, 9, 9, 4, 6, 4,
                                        5, 6, 6, 6, 6, 6, 7, 6, 5, 4,
                                        3, 7, 6, 6, 6, 7, 6, 6, 7, 3,
                                        2, 4, 7, 6, 7, 6, 6, 7, 4, 2,	
                                        1, 2, 3, 5, 6, 5, 4, 3, 2, 1
                                    ])