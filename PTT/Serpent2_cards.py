# Serpent2 card class structures in Python3
#  
# Instances of this class are initialized in the DMLG_Interface class when parsing Serpent2 cards
# Author: R. Guasch
# Part of the DMLG code Package.
#
# Contribution to BWR geometry handling for the Version5 environment.
# 
# Purpose : retrieve geometric info from parsed data and create Serpent2 cards object
# 

class S2_mat_card:
    """
    Definition of Serpent2 output geometry object.
    io_mode = output read or input read
    data = dict with key = material name 
                     values (Atom density, Mass density, Volume, Mass, + composition pairs (iso, a. dens))
    """

    def __init__(self, name, data, io_mode):
        self.io_mode = io_mode
        #print(data)
        self.mat_name = name
        self.atom_dens = float(data.split(" ")[0])
        self.mass_dens = float(data.split(" ")[1])
        self.volume =  float(data.split(" ")[2])
        self.mass = float(data.split(" ")[3])
        # now retrieve compositions !
        print(f"Processing material name : {self.mat_name} with volume {self.volume}")
