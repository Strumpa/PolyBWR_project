# Serpent2 card class structures in Python3
#  
# Instances of this class are initialized in the DMLG_Interface class when parsing Serpent2 cards
# Author: R. Guasch
# Part of the DMLG code Package.
#
# Contribution to BWR geometry handling for the Version5 environment.
# 
# Purpose : retrieve geometric info from parsed data and create Serpent2 cards objects, material properties objects and geometry objects.
# 

import os
import numpy as np

class S2_mat_output:
    def __init__(self, name, data):
        """
        Definition of Serpent2 output geometry object.
        data = dict with key = material name 
        values (Atom density, Mass density, Volume, Mass, + composition pairs (iso, a. dens))
        """
        #print(data)
        self.mat_name = name
        self.atom_dens = float(data.split(" ")[0])
        self.mass_dens = float(data.split(" ")[1])
        self.volume =  float(data.split(" ")[2])
        self.mass = float(data.split(" ")[3])
        # now retrieve compositions !
        print(f"Processing material name : {self.mat_name} with volume {self.volume}")



class S2_mat_properties:
    def __init__(self, name, S2_materials, nbDim):
        """
        name = name of the geometry, user defined
        S2_materials = list of S2_material_output objects contained in the geometry
        """
        
        self.Geo_name = name
        self.S2materials = S2_materials
        self.nbDim=nbDim
        if self.nbDim == 2:
            self.units = "cm2"
        elif self.nbDim ==3 :
            self.units = "cm3"

        print(f"The Serpent2 geometry being processed is {self.nbDim}-Dimensional, the 'volume' units associated are {self.units}")
        self.ComputeTotalFuelMass()
        self.ComputeTotalFuelVolume()
        self.ComputeTotalFuelMassDens()
        

    def getMaterialsandVolumes(self):
        material_vol_dict = {}
        for material in self.S2materials:
            material_vol_dict[material.mat_name]=material.volume
        return material_vol_dict
    
    
    def getMaterialsandMass(self):
        material_Mass_dict = []
        for material in self.S2materials:
            material_Mass_dict[material.mat_name]=material.mass
        return material_Mass_dict
    
    def getMaterialsandMassDens(self):
        material_MassDens_dict = []
        for material in self.S2materials:
            material_MassDens_dict[material.mat_name]=material.mass_dens 
        return material_MassDens_dict


    # Compute and get total fuel volume for a given geometry
    def ComputeTotalFuelVolume(self):
        fuel_vol = 0
        for material in self.S2materials:
            if "UOx" in material.mat_name or "fuel" in material.mat_name or "FUEL" in material.mat_name:
                fuel_vol+=material.volume
        self.TotalFuelVol = fuel_vol
        return
    
    def getFuelVolume(self):
        return self.TotalFuelVol
    

    # Compute and get total fuel mass for the given geometry
    def ComputeTotalFuelMass(self):
        fuel_mass = 0
        for material in self.S2materials:
            if "UOx" in material.mat_name or "fuel" in material.mat_name or "FUEL" in material.mat_name:
                fuel_mass+=material.mass
        self.TotalFuelMass=fuel_mass
        return
    
    def getTotalFuelMass(self):
        return self.TotalFuelMass

    # Compute and get total fuel mass density for the given geometry
    def ComputeTotalFuelMassDens(self):
        self.TotalFuelMassDens = self.TotalFuelMass/self.TotalFuelVol
        return 
    
    def getTotalFuelMassDens(self):
        return self.TotalFuelMassDens
    
class S2_Material_Vol:
    def __init__(self, name, data, nbDim):
        self.mode = "MVol"
        self.mat_vol_name = name
        self.nbDim=nbDim
        if self.nbDim == 2:
            self.units = "cm2"
        elif self.nbDim ==3 :
            self.units = "cm3"
        self.volume = data
        print(f"{self.mat_vol_name} with volume {self.volume} {self.units} processed")
    
class S2_geom:
    def __init__(self, name, data, nbDim, height):
        """
        Minimal implementation of the S2 geometry :
        at this stage only material volumes are taken to be included in data
        """
        self.geo_name = name
        self.nbDim = nbDim
        self.height = height
        material_volumes = {}
        #print(data)
        #for mat in data.keys():
        for elem in data:
            if elem.mode == "MVol":
                material_volumes[elem.mat_vol_name] = float(elem.volume)
        self.material_volumes = dict(sorted(material_volumes.items(), key=lambda item: item[1]))
        #print(len(self.material_volumes))
    
    def getOrderedMaterialVols(self):
        return self.material_volumes


class S2_case:
    def __init__(self, name, input_data_type, input_case, mode):
        """
        name (str) : name of the S2_case to be output
        input_data = input data parsed from MCNP (or Dragon) input file to create equivalent Serpent2 case. Either MCNP_case object or Dragon5_case object as created from DMLG_Interface.
        input_data_type = type of input data (MCNP or Dragon)
        """
        self.name = name
        self.input_data_type = input_data_type
        self.mode = mode

        self.fuel_pins = []
        self.cell_cards = []
        self.surface_cards = []
        self.material_cards = []
        #self.lattice_card = []

        print(f"Processing output data for {self.name}")
        if input_data_type == "MCNP" and self.mode == "output":
            self.cell_cards, self.fuel_pins, self.surface_cards, self.material_cards, self.lattice_card = self.convert_MCNP_to_S2(input_case)
        # attributes should be added to the class to allow for the creation of an output serpent2 case and write its equivalent text file.


        
    
    def convert_MCNP_to_S2(self, input_case):
        """
        This function converts an MCNP input case to a Serpent2 input case.
        """
        cell_cards = []
        fuel_pins = []
        cells_that_go_in_fuel_pins = []
        surface_cards = []
        material_cards = []
        lattice_card = []

        for cell in input_case.cell_cards:
            if "pin" in cell.cell_name:
                cells_that_go_in_fuel_pins.append(cell)
            else:
                cell_cards.append(S2_cell(cell.cell_name, 0, cell.material_number, cell.surfaces))
        for surface in input_case.surface_cards:
            surface_cards.append(S2_surface(surface.surface_number, surface.surface_type, surface.surface_data))
        for material in input_case.material_cards:
            if "t" in material.material_name:
                material_cards.append(S2_material(material.material_name, 0, [], material.iso_codes, material.therm_lib))
            else:
                material_cards.append(S2_material(material.material_name, "sum", [], material.iso_codes, material.iso_densities))
        
        # now create the fuel pins
        fuel_pins = self.group_cells_in_fuel_pins(cells_that_go_in_fuel_pins)
        


        return cell_cards, fuel_pins, surface_cards, material_cards, lattice_card

    def group_cells_in_fuel_pins(self, cells_that_go_in_fuel_pins):
        """
        This function groups the cells that go in fuel pins together to create fuel pins.
        """
        fuel_pin_groups = {}
        for cell in cells_that_go_in_fuel_pins:
            if cell.cell_name not in fuel_pin_groups.keys():
                fuel_pin_groups[cell.cell_name] = []
            fuel_pin_groups[cell.cell_name].append(cell)
        pins = []
        pin_number = 1
        for pin_name, pin_list in fuel_pin_groups.items():
            print(f"Processing fuel pin {pin_name}")
            materials = set()
            radii = set()
            for cell in pin_list:
                #print(f"Processing cell {cell.cell_name}")
                #print(f"Cell surfaces are : {cell.surfaces}")
                materials.add(cell.material_number)
                radii.add(cell.surfaces[0][-1])
            pins.append(S2_fuel_pin(pin_number, list(materials), list(radii)))
            pin_number += 1
        return pins
    

    def print_cells_info(self, output_path, output_name):
        """
        this function writes to a text file the cell cards information for the Serpent2 case.
        """

        if not os.path.exists(output_path):
            os.mkdir(output_path)
        output = open(output_path+output_name, "w")

        for pin in self.fuel_pins:
            print(f"Fuel pin {pin.pin_id} with materials {pin.materials} and radii {pin.radii}")
            output.write(f"pin {pin.pin_id}\n")
            if len(pin.materials) == len(pin.radii):
                for i in range(len(pin.materials)):
                    output.write(f"{pin.materials[i]} {pin.radii[i]}\n")
            else:
                for i in range(len(pin.radii)):
                    output.write(f"{pin.materials[i]} {pin.radii[i]}\n")
                output.write(f"{pin.materials[-1]}")                
        for cell in self.cell_cards:
            print(f"Cell name : {cell.cell_name} with universe number {cell.universe_nb} and material {cell.material}")
            print(f"Cell surfaces are : {cell.surfaces}")
            output.write(f"{cell.cell_type} {cell.cell_name} {cell.universe_nb} {cell.material} {cell.surfaces}\n")


class S2_cell:
    def __init__(self, name, universe_nb, material, surfaces_list):
        """
        Definition of Serpent2 cell object.
        Attributes are :
        cell type, cell name, universe number, material, surfaces list
        """
        #self.cell_type = cell_type
        self.cell_name = name
        self.universe_nb = universe_nb
        self.material = material
        self.surfaces = surfaces_list
        
        print(f"Processing cell name : {self.cell_name}")
    
class S2_fuel_pin:
    def __init__(self, pin_id, materials, radii):
        """
        initializing fuel pin object.
        attributes are pin_id = universe number, materials = list of materials, radii = list of radii defining the fuel pin
        """

        self.pin_id = pin_id
        self.materials = materials
        self.radii = radii
        print(f"Processing fuel pin {self.pin_id}")
        return
class S2_surface:
    def __init__(self, surface_id, type, parameters):
        """
        Definition of Serpent2 surface object.
        Attributes are :
        surface id, type of surface, parameters defining the surface
        """
        self.surface_id = surface_id
        self.surface_type = type #type can be px, py, pz, sph, cylz or cyl, cylx, cyly, sqc, cube, cuboid, cross.
        self.surface_parameters = parameters # the parameter array definies the surface according to mathematical definition in "Serpent user manual"
        print(f"Processing surface {self.surface_id} of type {self.surface_type}")
        return

class S2_material:
    def __init__(self, name, dens, options, isotopes, data):
        """
        Definition of Serpent2 material object.
        Attributes are : material name, density, options, isotopes, data
        unit is an attribute that is deduced from the data
        """
        self.material_name = name
        print(f"In S2 material definition : Processing material : {self.material_name}")
        print(f"Associated isotopes are = {isotopes}")
        print(f"Associated received data is = {data}")
        if "t" in self.material_name:
            self.mat_type = "thermalizer"
            self.dens = 0
            self.options = options
            self.isotopes = isotopes
            self.data = data
        else:    
            if dens == "sum" or dens == 0:
                self.density = "sum"
                if data[0] > 0 :
                    self.units = "atomic density units (10^24/cm3)"
                if data[0] < 0:
                    self.units = "mass density units (g/cm3)"
            elif dens < 0:
                self.density = dens
                self.units = "mass fractions"
            elif dens > 0:
                self.density = dens
                self.units = "atom fractions"

            self.options = options
            self.isotopes = isotopes
            self.data = data

        # check that data has the same length as isotopes
        if len(data) != len(isotopes):
            print("Data and isotopes do not have the same length")
            return
        
        return