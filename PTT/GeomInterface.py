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
    def __init__(self, input_deck, type):
        """
        loading geometry object to be parsed according to :
        input file (str) : path to file to be parsed
        type (str) : "MCNP", "Native", "Serpent2"
        """
        self.type = type
        self.input_deck = input_deck
        if type == "MCNP":
            self.parseMCNP_deck()
            self.createMCNPcard_objects()
        elif type == "Native_Dragon":
            self.parseNative_deck()
        elif type == "Serpent2":
            self.parseSerpent2_deck()
        else:
            print("Error: invalid Geometry type to be parsed")
        
        
    def parseMCNP_deck(self):
        """
        parsing the MCNP input file
        Its main structure contains a defintion of "Cell Cards", "Surface Cards" and "Material Cards" 
        Each section has a dictionnary associated in order to store the info about each card.
        The dictionnarys are :
        for Cell Cards : group of cells with associated data
        for Surface Cards : group of surfaces with asociated data
        for Material Cards : material name with associated data

        Classes defining the proper data structure for each card/ group of cards are called in this function.
        """
        print(f"Parsing geometry of type : {self.type}, from input file : {self.input_deck}")
        # Define flags to track when to start/stop parsing each section
        parsing_cell_cards = False
        parsing_surface_cards = False
        parsing_material_cards = False

        # Initialize dictionaries to store information from each section
        cell_cards = {}
        surface_cards = {}
        material_cards = {}

        # Initialize variables to store current card title and contents
        current_title = None
        current_contents = []

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

                # Parse and store information from each section
                elif parsing_cell_cards or parsing_surface_cards or parsing_material_cards:
                    if line:  # Skip empty lines
                        # Check if the line indicates a new card title
                        if line.startswith("c"):
                            # Store previous card's contents (if any) in the corresponding dictionary
                            if current_title and current_contents:
                                if parsing_cell_cards:
                                    cell_cards[current_title] = current_contents.copy()
                                elif parsing_surface_cards:
                                    surface_cards[current_title] = current_contents.copy()
                                elif parsing_material_cards:
                                    material_cards[current_title] = current_contents.copy()

                            # Update current card title and reset contents list
                            current_title = line.split("    ")[-1]
                            current_contents = []
                        else:
                            # Add line to current card's contents
                            current_contents.append(line.replace("  ", " ").replace("   ", " ").replace("    ", " ").replace("     ", " ").replace("      "," "))

            # Store contents of the last group of cards in each section
            if current_title and current_contents:
                if parsing_cell_cards:
                    cell_cards[current_title] = current_contents.copy()
                elif parsing_surface_cards:
                    surface_cards[current_title] = current_contents.copy()
                elif parsing_material_cards:
                    material_cards[current_title] = current_contents.copy()

        self.Cell_Cards = cell_cards
        self.Surface_Cards = surface_cards
        self.Material_Cards = material_cards
        #print(self.Surface_Cards)

    def createMCNPcard_objects(self):
        for cell_group_name in self.Cell_Cards.keys():
            self.Cell_Cards[cell_group_name] = MCNP_Cell_Card(1, cell_group_name, self.Cell_Cards[cell_group_name])
        for surface_group_name in self.Surface_Cards.keys():
            self.Surface_Cards[surface_group_name] = MCNP_Surface_Card(surface_group_name, self.Surface_Cards[surface_group_name],printlvl=True)
    
    def getMCNP_card_data(self, print_cells, print_surfaces, print_materials):
        """
        print_cells/surfaces/materials = boolean to specify print level
        """
        if print_cells:
            for cellcard in self.Cell_Cards.values():
                print(f"Cell card attrtibutes are : group of cells name = {cellcard.cells_group}, \n user defined cell numbers = {cellcard.cell_numbers}, \n material numbers : {cellcard.material_numbers} \n and material densities : {cellcard.material_densities} ")
                print(f"neutron importance is : {cellcard.neutron_importance}")
        if print_surfaces:
            for surfcard in self.Surface_Cards.values():
                print(f"Surface card attributes are : group of surfaces name : = {surfcard.surfaces_group}")
        if print_materials:
            print("Materials not implemented yet")





    def parseNative_deck(self):
        print("Native geom type not suported yet")
        return

    def parseSerpent2_deck(self):
        print("Serpent2 geom type not supported yet")
        return


input_f = "MCNP_AT10_sanitized.inp"
AT10_CRTL_MCNP = DMLG_Interface(input_f, type="MCNP")
#print(AT10_CRTL_MCNP.Cell_Cards)
AT10_CRTL_MCNP.getMCNP_card_data(print_cells=False, print_surfaces=False, print_materials=False)
print(AT10_CRTL_MCNP.Cell_Cards['water box centered at ( 8.267, 6.973)'].surfaces)
#print(AT10_CRTL_MCNP.Cell_Cards['water box centered at ( 8.267, 6.973)'].material_densities[2])