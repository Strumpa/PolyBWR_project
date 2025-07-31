## Main setup_case.py file to be used to set up a neutronics BWR case
# Idea is to set up a case for given a geometry and material composition
# Options could involve setting up a Serpent2 input file, a DRAGON5 case (collection of c2m files) 
# Options could allow for different types of calculations.

# General philosophy : yaml files for configuration, python for logic and setting up the case.
import yaml
from pathlib import Path
from DMLG_composition_handling.material_mixture import Material_Mixture



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
        LIB_definition = (
            "LIBRARY := LIB: ::\n",       
                          )
        if self.settings["DRAGON5"]["Self-Shiedling Method"]:
            LIB_definition += (f'{self.settings["DRAGON5"]["Self-Shiedling Method"]}',)
        if self.settings["DRAGON5"]["Order Anisotropic Scattering"]:
            LIB_definition += (f'ANIS {self.settings["DRAGON5"]["Order Anisotropic Scattering"]}',)
        if self.settings["DRAGON5"]["Number of mixtures"]:
            LIB_definition += (f'NMIX {self.settings["DRAGON5"]["Number of mixtures"]}',)
        if self.settings["DRAGON5"]["DRAGLIB"]:
            LIB_definition += (f'MIXS LIB: DRAGON FIL: {self.settings["DRAGON5"]["DRAGLIB"]}',)
        if self.settings["DRAGON5"]["Depletion"] and self.settings["DRAGON5"]["DRAGLIB"]:
            LIB_definition += (f'DEPL LIB: DRAGON FIL: {self.settings["DRAGON5"]["DRAGLIB"]}',)

        for material_name, material in self.materials.items():
            if material.isdepletable:
                LIB_definition += (f'MIX {material_name.split(" ")[0]} {material_name.split(" ")[-1]}',)
            else: 
                LIB_definition += (f'MIX {material_name.split(" ")[0]} {material_name.split(" ")[-1]} NOEV',)
            LIB_definition += (material.print_to_D5_format(),)

        with open(c2m_file_path, 'w') as c2m_file:
            for line in LIB_definition:
                c2m_file.write(line)

        return
            
            
    

if __name__ == "__main__":
    # Test usage on OECD PHASE IIIB BWR assembly benchmark case :
    case = DMLG_case(name="OECD_NEA_PHASE_IIIB", output_codes=["Serpent2", "DRAGON5"], nuclear_data_evaluation="endfb8r1")
    case.settings = {"DRAGON5": {"DRAGLIB": "endfb8r1_295",
                                "Self-Shiedling Method": "PT",
                                "Order Anisotropic Scattering": 3,
                                "Number of mixtures": 10,
                                "Depletion": True}}
    case.print_materials_to_D5()