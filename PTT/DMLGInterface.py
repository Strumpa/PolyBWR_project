## Python3 script for parsing MCNP input file
# Author : R. Guasch
# Idea : convert MCNP geom to Serpent2
# Try to keep it general : idea could be to use functions in centralized package for geometry handling
# Goals :
# Step (1) : Convert MCNP input file for ATRIUM10 with control blade to Serpent2 format
# Step (2) : define a general Geometry object that could be written in Serpent2 or Dragon5 native geom
# Step (3) : allow for surfacic geometry ?

import numpy as np
import math
from GeometricTools import geom_PIN
from MCNP_card import *

# Class DMLG_Interface()

class DMLG_Interface:
    def __init__(self, input_deck, type, mode):
        """
        loading geometry object to be parsed according to :
        input file (str) : path to file to be parsed
        type (str) : "MCNP", "Native", "Serpent2"
        mode (str) : 'input' or 'output'
        """
        self.type = type
        self.input_deck = input_deck
        self.mode = mode
        if self.type == "MCNP" and self.mode == 'input' :
            self.parseMCNP_input()
            self.createMCNPcard_objects()
        elif self.type == "Native_Dragon" and self.mode == 'input' :
            self.parseNative_input()
        elif self.type == "Serpent2" and self.mode == 'input' :
            self.parseSerpent2_input()
        elif self.type == "Serpent2" and self.mode == 'output' : 
            self.parseSerpent2_output()
        else:
            print("Error: invalid Geometry type to be parsed")
        
        
    def parseMCNP_input(self):
        """
        parsing the MCNP input file
        Its main structure contains a defintion of "Cell Cards", "Surface Cards" and "Material Cards" 
        Each section has a list associated in order to store the info about each card.
        The lists are :
        for Cell Cards : group of cells with associated data
        for Surface Cards : group of surfaces with asociated data
        for Material Cards : material name with associated data

        Classes defining the proper data structure for each card/ group of cards are called in this function.
        """
        print(f"Parsing {self.type} {self.mode} file, from input file : {self.input_deck}")
        # Define flags to track when to start/stop parsing each section
        parsing_cell_cards = False
        parsing_surface_cards = False
        parsing_material_cards = False

        # Initialize dictionaries to store information from each section
        cell_cards = []
        surface_cards = []
        material_cards = []
        cell_titles = []
        surface_titles = []
        material_titles = []

        # Initialize variables to store current card title and contents
        current_title = None

        # Open the input file
        with open(self.input_deck, 'r') as file:
            for line in file:
                line = line.strip()  # Remove leading/trailing whitespace

                # Check if the line indicates the start of a section
                if "Cell  Cards" in line:
                    parsing_cell_cards = True
                    parsing_surface_cards = False
                    parsing_material_cards = False
                elif "Surface Cards" in line:
                    parsing_cell_cards = False
                    parsing_surface_cards = True
                    parsing_material_cards = False
                elif "Material  Cards" in line:
                    parsing_cell_cards = False
                    parsing_surface_cards = False
                    parsing_material_cards = True
                elif "Tally Cards" in line:
                    parsing_cell_cards = False
                    parsing_surface_cards = False
                    parsing_material_cards = False
                
                elif parsing_cell_cards or parsing_surface_cards:
                    if line.startswith("c"):
                        if line.lstrip("c").strip():
                            current_title=line.lstrip("c").strip()

                # Parse and store information from each section
                elif parsing_cell_cards or parsing_surface_cards and line:
                    # Store previous card's contents (if any) in the corresponding dictionary
                    if parsing_cell_cards:
                        cell_cards.append(line.replace("  ", " ").replace("   "," "))
                        cell_titles.append(current_title)
                    elif parsing_surface_cards:
                        surface_cards.append(line.replace("  ", " ").replace("   "," "))
                        surface_titles.append(current_title)
                elif parsing_material_cards and line:
                    if line.startswith("m"):
                        line=line.replace("  ", " ").replace("   "," ")
                        current_title=line.split(" ")[0]
                        material_titles.append(current_title)
                    if not line.startswith("c"):
                        line=line.replace("  ", " ").replace("   "," ")
                        material_cards.append(line)
        file.close()              

 

        self.Cell_Cards = self.clean_list(cell_cards)
        self.Surface_Cards = self.clean_list(surface_cards)
        self.cell_titles = self.clean_list(cell_titles)
        self.surface_titles = self.clean_list(surface_titles)
        
        self.Material_Cards = self.combine_list(material_cards) 
        
    def clean_list(self,list_):
        cleaned_list=[]
        for entry in list_:
            if entry != '' or entry !="\n":
                cleaned_list.append(entry.strip())
        return cleaned_list
    def combine_list(self, list_):
        combined_list=[]
        current_list=[]
        for sub_list in list_:
            sub_list = sub_list.split(" ")
            if sub_list[0][0]!="m" and current_list:
                current_list.extend(sub_list)
            else:
                if current_list:
                    combined_list.append(current_list)
                current_list=sub_list
        if current_list:
            combined_list.append(current_list)
        return combined_list

    def createMCNPcard_objects(self):
        for i in range(len(self.Cell_Cards)):
            self.Cell_Cards[i] = MCNP_Cell_Card(1, self.cell_titles[i], self.Cell_Cards[i])
        for i in range(len(self.Surface_Cards)):
            self.Surface_Cards[i] = MCNP_Surface_Card(self.surface_titles[i], self.Surface_Cards[i],printlvl=False)
        for i in range(len(self.Material_Cards)):
            self.Material_Cards[i] = MCNP_Material_Card(self.Material_Cards[i], printlvl=True)
    
    def getMCNP_card_data(self, print_cells, print_surfaces, print_materials):
        """
        print_cells/surfaces/materials = boolean to specify print level
        """
        if print_cells:
            print(f"There are a total of {len(self.Cell_Cards)} cell cards")
            for cellcard in self.Cell_Cards:
                print(f"Cell card attrtibutes are : group of cells name = {cellcard.cells_group}, \n user defined cell numbers = {cellcard.cell_number} \n material numbers : {cellcard.material_number} \n and material densities : {cellcard.material_density} ")
                print(f"neutron importance is : {cellcard.neutron_importance}")
        if print_surfaces:
            print(f"There are a total of {len(self.Surface_Cards)} surface cards")
            for surface_card in self.Surface_Cards:
                print(f"Surface card attributes are : group of surfaces name : = {surface_card.surfaces_group}")
                print(f"The corresonding surface equation is : {surface_card.surfaceEquation}")
        if print_materials:
            print(f"There are a total of {len(self.Material_Cards)} Material cards")
        
        return





    def parseNative_input(self):
        print("Native geom type not suported yet")
        return

    def parseSerpent2_input(self):
        print("Serpent2 geom type not supported yet")
        return
    def parseSerpent2_output(self):
        """
        parsing Serpent2 output file

       
        """
        lines=[]
        materials=[]
        material_data_={}
        current_material=None
        parsing_new_mat = False
        parsing_same_mat = False
        print(f"Parsing {self.type} {self.mode} file, from input file : {self.input_deck}")
        with open(self.input_deck, 'r') as file:
            for line in file:
                line = line.strip()  # Remove leading/trailing whitespace
                #print(line)
                if "Material" in line and ":" in line:
                    current_material=line.split(" ")[1].strip('"').strip(":").strip('"')
                    materials.append(current_material)
                    parsing_new_mat=True
                #print(current_material)
                #print(current_material in material_data_.keys())
                elif parsing_new_mat:
                    #print("Parsing new line")
                    #print(line)
                    material_data_[current_material] = line
                    parsing_new_mat=False
                    parsing_same_mat=True
                elif parsing_same_mat:
                    if "density" in line :
                        data = material_data_[current_material] + line.split(' ')[3]+" "
                    elif "Volume" in line or "Mass" in line :
                        data = material_data_[current_material] + line.split(' ')[2]+" "
                    elif ".06c" in line or ".05c" in line or ".03c" in line:
                        data = material_data_[current_material] + line.split('  ')[0]+" "+line.split('  ')[4]+" "
                    else:
                        data = material_data_[current_material]
                    #print("combinin new line with preivous data")
                    #print(data)
                    material_data_.pop(current_material)
                    material_data_[current_material]=data
                """
                if current_material in material_data_.keys():
                    data=material_data_[current_material]+line
                    material_data_.pop(current_material)
                    material_data_[current_material]=data
                elif line:
                    material_data_[current_material]=line
                """
        file.close()
        print(materials)
        print(material_data_)

        return

"""
input_f = "MCNP_AT10_sanitized.inp"
AT10_CRTL_MCNP = DMLG_Interface(input_f, type="MCNP",mode="input")
#print(AT10_CRTL_MCNP.Cell_Cards)
AT10_CRTL_MCNP.getMCNP_card_data(print_cells=False, print_surfaces=False, print_materials=True)
#print(AT10_CRTL_MCNP.Cell_Cards)
#print(AT10_CRTL_MCNP.Cell_Cards['water box centered at ( 8.267, 6.973)'].material_densities[2])
"""
cell_output = "Serpent2/data_AT10_24UOX_try/AT10_24UOX_mc.out"
AT10_24UOX_cell = DMLG_Interface(cell_output, type="Serpent2",mode="output")