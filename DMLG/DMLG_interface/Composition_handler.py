## Python3 Script Composition_handler tool
# Author : R. Guasch
# Purpose : dealing with Serpent2 isotope codes and converting to Dragon5 LIB: MIX format,
# This version is an extension of ProcesComposition.py, with additional functionalities to homogenize geometric features.

class Compo_Processor:
    def __init__(self, mix_name, filename): 
        """
        filename : text file to read from to convert Serpent2 composition to Dragon5 format
        more code to isotopes to be added in order to allow for processing of all compos
        """
        self.mix_name = mix_name
        self.association_dict = {
                        "1001" : "H1", "1002" : "H2", "1003" : "H3",
                        "2003" : "He3", "2004" : "He4",
                        "8016" : "O16", "8017" : "O17", 
                        "24050": "Cr50", "24052": "Cr52", "24053": "Cr53", "24054": "Cr54", 
                        "26054": "Fe54", "26056": "Fe56", "26057": "Fe57", "26058": "Fe58",
                        "28058": "Ni58", "28060": "Ni60", "28061": "Ni61", "28062": "Ni62", "28064": "Ni64",
                        "40090": "Zr90", "40091": "Zr91", "40092": "Zr92", "40094": "Zr94", "40096": "Zr96",
                        "50112": "Sn112", "50114": "Sn114", "50115": "Sn115", "50116": "Sn116","50117": "Sn117",  "50118": "Sn118", "50119": "Sn119", "50120": "Sn120", "50122": "Sn122", "50124": "Sn124",
                        "64154": "Gd154", "64155": "Gd155", "64156": "Gd156", "64157": "Gd157", "64158": "Gd158", "64160": "Gd160",
                        "92235": "U235", "92238": "U238", "92234": "U234", "92236": "U236"}
   
        lines = open(filename, "r")
        isos_code_adens = {}
        for line in lines:
            line = " ".join(line.split())
            if "--" not in line and "Nuclide" not in line:
                isos_code_adens[line.split(" ")[0].split(".")[0]] = float(line.split(" ")[3])
        isos_code_adens.pop("sum")
        self.compo = isos_code_adens
    
    def check_consistency(self):
        incompatible_isos = []
        for iso_code in self.compo.keys():
            if iso_code not in self.association_dict.keys():
                print("Isotope code "+ iso_code + " is in compo but not in association dict") 
                incompatible_isos.append(iso_code)
        if len(incompatible_isos) == 0 :
            print("No issues with isotope compatibility detected, process may continue.")
        return
    
    def iso_code_to_nuclide(self):
        dict_Adens = {}
        for iso_code in self.compo.keys():
            if iso_code in self.association_dict.keys():
                dict_Adens[self.association_dict[iso_code]] = self.compo[iso_code]
        self.Ndens_isos = dict_Adens
        return

    def print_to_Dragon5_format(self):
        """
        rewrite the info in Dragon iso = iso [atomic density] format
        isotopes as keys and atomic density as values

        condition to add autop options ? To implement ?
        """
        for iso in self.Ndens_isos.keys():
            print(str(iso)+"  = "+str(iso)+"   "+str((self.Ndens_isos[iso])))


    def get_enrichment(self,U_isotopic_compo):
        """
        U_isotopic_compo (dict) = isotopic vector with values [N_i] = isotopic number of atoms / (b*cm)
                                                with keys = isotope name
        Assuming enrichment is calculated as N_U235/(sum N_Us)
        """
        Ntot_U = 0
        for iso in U_isotopic_compo.keys():
            Ntot_U += U_isotopic_compo[iso]
        return U_isotopic_compo["U235"]*100/(Ntot_U)

# Idea of class for library generation : write a class that takes a list of Compo_Processor objects and writes a Dragon5 LIB: structure definition or a Serpent2 material definition 




