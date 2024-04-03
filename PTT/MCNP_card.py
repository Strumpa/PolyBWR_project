# MCNP card class structures in Python3
#  
# Instances of this class are initialized in the DMLG_Interface class when parsing MCNP cards
# Author: R. Guasch
# Part of the DMLG code Package.
#
# Contribution to BWR geometry handling for the Version5 environment.
# 
# Structure of the Cell card object obtained from "MCNP® USER’S MANUAL Code Version 6.2 October 27, 2017"
# 

import numpy as np
import re
from sympy import symbols, Eq
from sympy.geometry import Plane
#from sympy.geometry import Cylinder

class MCNP_Cell_Card:
    """
    definition of cell cards from MCNP format
    """
    def __init__(self, format, cell_group_name, data):
        """
        format = input format , only format  is supported
        cell name (str) (optional) : to be edited in case cell cards do not procide comments with cells' names
        data : list of lists parsed in "parseMCNP_deck", data will be accessed and used to intialize de cell cards' attributes
        """
        if format == 1:
            self.format = format
        else:
            print("Input cell format to supported yet")
            return
        self.cells_group = cell_group_name
        self.surfaces = []
        self.neutron_importance=np.zeros(len(data))
        data=data.replace("  ", " ").replace("   ", " ")
        data = data.split(" ")
        cleaned_data=[]
        for entry in data:
            if entry != '' and data:
                cleaned_data.append(entry)
        self.cell_number = int(cleaned_data[0])
        self.material_number = int(cleaned_data[1])
        self.material_density = float(cleaned_data[2])
        self.surfaces.append(cleaned_data[3:-4]) # for ith material in cell : append list of surfaces
        self.neutron_importance=int(cleaned_data[-4][-1])

class MCNP_Surface_Card:
    """
    Defintion of surface card object according to MCNP format
    general format : J N A list
    J = surface number given by user, if * or + opertor present before surface number : specular or reflective BC to be applied to surface
    N = optional entry specifying a coordinate transformation with TRn card or periodic translation. Not implemented here.
    A = surface type (plane, cylinder, etc)
    list = corresponding coefficients defining surface equation.
    printlvl = boolean to chose to print verbose
    """
    def __init__(self, surface_group, data,printlvl):
        self.surfaces_group = surface_group # just a comment from INP file, carried in for reference. Might need to get rid of later ?
        # need to determine surface type. Should be given by second entry of surface card, assuming N entry not used. 
        self.printlvl=printlvl
        self.surface_data=[]
        data=data.replace("  ", " ").replace("   ", " ")
        data = data.split(" ")
        cleaned_data=[]
        for entry in data:
            if entry != '':
                cleaned_data.append(entry) 
        if self.printlvl:
            print(f"cleaned up surface type is : {cleaned_data[0]}")
        surface_number=cleaned_data[0]
        if "*" in cleaned_data[0]:
            self.BC_type="Reflective"
            self.surface_number=int(surface_number.replace("*",""))
        elif "+" in cleaned_data[0]:
            self.BC_type="White"
            self.surface_number=int(surface_number.replace("+",""))
        else:
            self.BC_type="Not a surface at a boundary"
            self.surface_number=int(surface_number)
        self.surface_type=cleaned_data[1]
        self.surface_data=cleaned_data[2:]
        if self.printlvl:
            print(f"cleaned up surface data is {cleaned_data[2:]}")
        self.setSurfaceEquation()
        
    def setSurfaceEquation(self):
        """
        Retrieve ceofficients from recovered surface data and set the right form for the equation        
        """
        if self.surface_type == "p":
            if self.printlvl:
                print(f"processing plane surface with equation given by Ax+By+Cz+D=0: A={self.surface_data[0]}, B={self.surface_data[1]}, C={self.surface_data[2]}, D={self.surface_data[3]}")
        elif self.surface_type == "px":
            if self.printlvl:
                print(f"processing plane surface with equation given by x-D=0: D={self.surface_data[0]}")
        elif self.surface_type == "py":
            if self.printlvl:
                print(f"processing plane surface with equation given by y-D=0: D={self.surface_data[0]}")
        elif self.surface_type == "pz":
            if self.printlvl:
                print(f"processing plane surface with equation given by z-D=0: D={self.surface_data[0]}")




