## Python3 Script ProcessSerpCompo
# Author : R. Guasch
# Purpose : dealing with Serpent2 isotope codes and converting to Dragon5 LIB: MIX format

class Compo_Processor:
    def __init__(self, filename): 
        """
        filename : text file to read from to convert Serpent2 composition to Dragon5 format
        
        
        more code to isotopes to be added in order to allow for processing of all compos
        """
        self.association_dict = {
                        "8016" : "O16", "8017" : "O17", 
                        "24050": "Cr50", "24052": "Cr52", "24053": "Cr53", "24054": "Cr54", 
                        "26054": "Fe54", "26056": "Fe56", "26057": "Fe57", "26058": "Fe58",
                        "28058": "Ni58", "28060": "Ni60", "28061": "Ni61", "28062": "Ni62", "28064": "Ni64",
                        "40090": "Zr90", "40091": "Zr91", "40092": "Zr92", "40094": "Zr94", "40096": "Zr96",
                        "50112": "Sn112", "50114": "Sn114", "50115": "Sn115", "50116": "Sn116","50117": "Sn117",  "50118": "Sn118", "50119": "Sn119", "50120": "Sn120", "50122": "Sn122", "50124": "Sn124"
                                 }
   
        lines = open(filename, "r")
        isos_code_adens = {}
        for line in lines:
            if "--" not in line:
                isos_code_adens[line.split("  ")[0].split(".")[0][1:]] = line.split("  ")[3]
                print(line.split("  ")[0].split(".")[0])
                print(line.split("  ")[3])
        isos_code_adens.pop("Nuclide")
        isos_code_adens.pop("")
        print(len(isos_code_adens))
        print(len(self.association_dict))
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
    


Processor = Compo_Processor("BOX_compo.txt")
Processor.check_consistency()
Processor.iso_code_to_nuclide()
Processor.print_to_Dragon5_format()

compo_fuel1 = {"U234": 5.15910E-06, "U235": 5.67035E-04, "U238": 2.27631E-02}
compo_fuel2 = {"U234": 7.039170e-06, "U235": 7.560370e-04, "U238": 2.257430e-02}
compo_fuel3 = {"U234": 9.163680e-06, "U235": 9.686590e-04, "U238": 2.236200e-02}
compo_fuel4 = {"U234": 9.991530e-06, "U235": 1.051340e-03, "U238": 2.227940e-02}
compo_fuel5 = {"U234": 1.058330e-05, "U235": 1.110400e-03, "U238": 2.222040e-02}
compo_fuel6 = {"U234": 1.117530e-05, "U235": 1.169460e-03, "U238": 2.216140e-02}
compo_fuel7 = {"U234": 9.451580e-06, "U235": 9.945290e-04, "U238": 2.107540e-02}
compo_fuel8 = {"U234": 8.668470E-06, "U235": 9.163120E-04, "U238": 2.115350E-02}
percentages = [Processor.get_enrichment(compo_fuel1), Processor.get_enrichment(compo_fuel2), Processor.get_enrichment(compo_fuel3),
               Processor.get_enrichment(compo_fuel4), Processor.get_enrichment(compo_fuel5), Processor.get_enrichment(compo_fuel6),
               Processor.get_enrichment(compo_fuel7), Processor.get_enrichment(compo_fuel8)]
print(percentages)




