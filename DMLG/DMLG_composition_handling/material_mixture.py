## Python3 class to define an abstract material mixture. 
# This class is used to handle the isotopic composition of a material mixture.
# Author : R. Guasch
# Date : 2025/07/22

import mendeleev
import re
import yaml


class Material_Mixture:
    """
    Abstract class to define a material mixture.
    This class is used to handle the isotopic composition of a material mixture.
    """

    def __init__(self, name: str, composition: dict, temperature:float, evaluation_name: str = None, isdepletable: bool = False):
        """
        Initialize the Material_Mixture object.

        :param name: Name of the material mixture.
        :param composition: Dictionary containing the isotopic composition of the material mixture.
        """
        self.therm_isotopes = {"H1": "H1_H2O"}
        self.HM_isotopes = ["U234", "U235", "U236", "U238", "Pu239", "Pu240", "Pu241", "Pu242", "Am241", "Am243", "Cm244", "Cm245", "Cm246"]
        self.name = name
        self.composition = composition
        self.set_mixture_temperature(temperature)
        if evaluation_name is not None:
            self.set_evalutation(evaluation_name)
        if isdepletable:
            self.isdepletable = True
        else:
            self.isdepletable = False

        self.convert_nat_to_iso()  # Convert natural isotopic composition to isotopic composition if needed

    def get_composition(self):
        """Return the isotopic composition of the material mixture."""
        return self.composition
    
    def get_name(self):
        """Return the name of the material mixture."""
        return self.name
    
    def set_composition(self, composition: dict):
        """
        Set the isotopic composition of the material mixture.

        :param composition: Dictionary containing the isotopic composition of the material mixture.
        """
        self.composition = composition

    def convert_nat_to_iso(self):
        """
        Convert the natural isotopic composition to isotopic composition.
        Updates the composition attribute with the isotopic values.
        """

        if not hasattr(self, 'composition'):
            raise ValueError("Composition is not defined for this material mixture.")
        
        # Check if any nuclides are present in the composition
        if not self.composition:
            raise ValueError("Composition is empty. Cannot convert natural to isotopic composition.")
        
        # Check if any nuclides are in "natural" element format :
        nat_nuclides = []
        iso_composition = {}
        for nuclide in self.composition.keys():
            ## Identify if the nuclide is in natural format : ends in 0 or in "_nat" or in "_NAT" or doesnt have a nucleon number and just atomic symbol.
            if nuclide.endswith('_nat') or nuclide.endswith('_NAT') or nuclide.isalpha():
                # try converting to isotopic composition
                nuclide_symbol = nuclide.replace('_nat', '').replace('_NAT', '').replace('0', '')
                nat_nuclides.append(nuclide_symbol)
                try:
                    ## use mendeleev to get the isotopes of element
                    isotopes = mendeleev.element(nuclide_symbol).isotopes
                    if not isotopes:
                        raise ValueError(f"No isotopes found for element {nuclide_symbol}.")
                    # create a dictionary to hold the isotopic composition
                    for iso in isotopes:
                        if iso.mass_number is None or iso.abundance is None:
                            continue  # skip isotopes without mass number or abundance
                        # print the conversion for debugging
                        print(f"Converting nuclide {nuclide} to isotopic composition: {iso.mass_number}, abundance: {iso.abundance}")
                        print(iso.mass_number, iso.abundance)
                        print(f"compostion - {self.composition}")
                        print(f"compute isotopic density for {nuclide_symbol}{iso.mass_number} = {self.composition[nuclide]*iso.abundance/100}")
                        iso_composition[f"{nuclide_symbol}{iso.mass_number}"] = self.composition[nuclide] * iso.abundance/100 # convert natural element density to densities of natural isotopic abundances
                
                    # remove the natural element from the composition
                except Exception as e:
                    raise ValueError(f"Error converting nuclide {nuclide_symbol} to isotopic composition: {e}")
                
        # remove the natural elements from the composition
        for nuclide in nat_nuclides:
            if nuclide in self.composition:
                del self.composition[nuclide]
        # update the composition with the isotopic composition
        for iso in iso_composition.keys():
            self.composition[iso] = iso_composition[iso]
                

    def format_zaid(self, iso: str):
        # special treatment for thermal scattering data
        if iso == "H1_H2O":
            zaid = 1001
            print(f"Treatment of {iso} as {zaid}, H1_H2O")
            return
        elif iso == "H2_D2O":
            zaid = 1002
            print(f"Treatment of {iso} as {zaid}, H2_D2O")
            return
        elif iso == "H1_CH2":
            zaid = 1001
            print(f"Treatment of {iso} as {zaid}, H1_CH2")
            return
        elif iso == "H1_ZRH":
            zaid = 1001
            print(f"Treatment of {iso} as {zaid}, H1_ZRH")
            return
        elif iso == "Zr90_ZrH":
            zaid = 40090
            print(f"Treatment of {iso} as {zaid}, Zr90_ZrH")
            return
        elif iso == "C12_GR":
            zaid = 6012
            print(f"Treatment of {iso} as {zaid}, C12_GR")
            return
        elif iso == "Be9":
            zaid = 4009
            print(f"Treatment of {iso} as {zaid}, Be9")
            return
        else:
            # special treatment for meta-stable isotopes
            if iso[-1] == "m":
                isMetaStable = True
                iso = iso[:-1]
                print(f"Treatment of {iso} as metastable")
            else:
                isMetaStable = False
            match = re.match(r"([A-Za-z]+)(\d+)", iso)
            if match:
                element = match.group(1)
                Z = mendeleev.element(element).atomic_number
                nucleons = int(match.group(2))
                #return {"element": element, "nucleons": nucleons}
            else:
                raise ValueError(f"Invalid isotope format: {iso}")
            if isMetaStable:
                zaid = 1000 * Z + 300 + nucleons%100
            else:
                zaid = 1000 * Z + nucleons
            #print(f"Treatment of {iso} as {zaid}")
            return zaid
        
    def set_evalutation(self, evaluation_name):
        """
        set the nuclear data evaluation to be used.
        """
        self.evaluation_name = evaluation_name
        print(f"Evaluation set to {self.evaluation_name} for material {self.name}")
        return
    
    # Get suffix for the closest inferior temperature
    def get_suffix(self, evaluation_name, temperature):
        """
        get suffix corresonding to the closest inferior temperature for a given evaluation.
        this allows to load the proper acefiles in Serpent2
        :param evaluation_name: Name of the evaluation (e.g., "endfb8r1").
        :param temperature: Temperature in Kelvin.
        :return: Suffix corresponding to the closest inferior temperature.
        :raises ValueError: If the evaluation or temperature is not found.
        """
        with open("eval_temp_suffix.yml", "r") as file:
            data = yaml.safe_load(file)
        for eval in data["Evaluation"]:
            if eval["name"] == evaluation_name:
                temps = eval["temperatures"]
                suffixes = eval["suffixes"]

                # Find all temperatures less than or equal to the given temperature
                valid_temps = [(t, suffixes[i]) for i, t in enumerate(temps) if t <= temperature]

                if not valid_temps:
                    raise ValueError(f"No temperature ≤ {temperature} found for {evaluation_name}")

                # Get the one with the maximum value (i.e., closest inferior)
                closest_temp, suffix = max(valid_temps, key=lambda x: x[0])
                return suffix

        raise ValueError(f"Evaluation '{evaluation_name}' not found")
    
    def get_evaluation(self):
        """
        get the nuclear data evaluation to be used.
        """
        if hasattr(self, 'evaluation_name'):
            return self.evaluation_name
        else:
            raise ValueError("Evaluation is not set for this material mixture.")
        
    def set_mixture_temperature(self, temperature):
        """
        Set the temperature of the material mixture.
        
        :param temperature: Temperature of the material mixture in Kelvin.
        """
        self.temperature = temperature
        return
        
    def print_to_S2_format(self):
        """
        Print the material mixture composition in S2 format.
        """
        print(f"Material: {self.name}")
        for iso, density in self.composition.items():
            zaid = self.format_zaid(iso)
            # if temperature is set, get the suffix for the evaluation
            if hasattr(self, 'temperature'):
                suffix = self.get_suffix(self.get_evaluation(), self.temperature)
                print(f"{zaid}{suffix} {density:.6E}")
            else:
                print(f"{zaid} {density:.6E}") 

    def print_to_D5_format(self):
        """
        Print the material mixture definition to DRAGON5 format.
        """
        print(f"Material: {self.name}")
        dens_definitions = ()
        for iso, density in self.composition.items():
            if iso in self.therm_isotopes:
                # Special treatment for thermal scattering data
                dens_definitions += (f"{iso}    = {self.therm_isotopes[iso]}  {density:.6E}",)
            elif iso in self.HM_isotopes:
                # Special treatment for elf shielding of heavy metal isotopes
                dens_definitions += (f"{iso}    = {iso}  {density:.6E}  1",)
            else:
                dens_definitions += (f"{iso}    = {iso}  {density:.6E}",)

        return dens_definitions



        

if __name__ == "__main__":

    # Test / Example usage :
    # Clad / box from OECD/NEA/Phase IIIB BWR benchmark :
    clad = Material_Mixture("clad", {})
    clad.set_composition({"Cr": 7.5891E-05, "Fe":1.4838E-04, "Zr": 4.2982E-02})
    clad.convert_nat_to_iso()
    clad.set_evalutation("endfb8r1")
    clad.set_mixture_temperature(559.0)  # example temperature in Kelvin
    #print(clad.get_composition())
    clad.print_to_S2_format()