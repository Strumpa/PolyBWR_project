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

        self.pins = []
        self.cell_cards = []
        self.surface_cards = []
        self.material_cards = []
        #self.lattice_card = []
        print("$$$ DMLG: Serpent2 case class initialized")
        print("$$ Generating Serpent2 case object")
        print(f"Processing output data for {self.name}")
        if input_data_type == "MCNP" and self.mode == "output":
            self.cell_cards, self.pins, self.surface_cards, self.material_cards, self.lattice_cards = self.convert_from_MCNP(input_case)
        # attributes should be added to the class to allow for the creation of an output serpent2 case and write its equivalent text file.

        return

        
    
    def convert_from_MCNP(self, input_case):
        """
        This function converts an MCNP input case to a Serpent2 input case.
        """
        cell_cards = []
        fuel_pins = []
        cells_to_pins = []
        self.surface_cards = []
        material_cards = []
        lattice_cards = []

        for cell in input_case.cell_cards:
            if "pin" in cell.cell_name:
                cells_to_pins.append(cell)
            else:
                cell_cards.append(S2_cell(cell.cell_name, 0, cell.material_number, cell.surfaces))
        for surface in input_case.surface_cards:
            self.surface_cards.append(S2_surface(surface.surfaces_group, surface.surface_number, surface.surface_type, surface.surface_data))
        for material in input_case.material_cards:
            if "t" in material.material_name:
                material_cards.append(S2_material(material.material_name, 0, [], material.iso_codes, material.therm_lib))
            else:
                material_cards.append(S2_material(material.material_name, "sum", [], material.iso_codes, material.iso_densities))

        # Analyse cells to create fuel pins
        fuel_pins = self.group_cells_in_pins(cells_to_pins)
        merged_pins, mat_to_pins_dict = self.mergePins(fuel_pins)
        print(f"Material to correpsonding dictionnary before merge is {mat_to_pins_dict}")
        for pin in merged_pins:
            print(f"Fuel pin {pin.pin_id} with materials {pin.materials} and radii {pin.radii}")

        # Analyse cartesian geometry to create the bounding surfaces and boxes  
        self.cell_cards = self.AnalyseGeometry(cell_cards)
        for cell in self.cell_cards:
            print("$$$ Geometry analysis")
            print(f"Cell name : {cell.cell_name} with universe number {cell.universe_nb} and material {cell.material}")
            print(f"Cell surfaces are : {cell.surfaces}")



        # Combine info on fuel pins and surfaces to create the pin lattice.
        fuel_lattice_card = self.createLatticeFromPins(fuel_pins)
        print(f"Fuel lattice is {fuel_lattice_card}")

        return self.cell_cards, merged_pins, self.surface_cards, material_cards, lattice_cards
    
    def AnalyseGeometry(self, cell_cards):
        """
        Combine information about cells and surfaces to create bounding surfaces and boxes.
        """
        for cell in cell_cards:
            print(cell.surface_ids)
            for surf in cell.surface_ids:
                print(surf)
                cell.setSurfaces(self.getSurfacesFromIds(cell.surface_ids))
                #print(f"Surface {surf.surface_id} of type {surf.surface_type} in group {surf.surface_group}")
        # Keep working on this function to create bounding surfaces and boxes.
                
        return cell_cards


    def createLatticeFromPins(self, pins, type="square",diagonal_symmetry=True, box = "none"):
        """
        This function creates a lattice card from the fuel pins.
        assuming square BWR lattice.
        Need to figure out box indices.
        """
        pins_lat = self.sortPins(pins)
        return pins_lat


    def sortPins(self, pins):
        """
        Sort pins according to increasing x postion and increasing y position
        """
        centers = []
        for pin in pins:
            centers.append(pin.center)
        sorted_centers = sorted(centers, key=lambda center: (center[0],center[1]))
        
        """
        # Extract unique x and y coordinates
        unique_x = sorted(set(point[0] for point in sorted_centers))
        unique_y = sorted(set(point[1] for point in sorted_centers))
        
        # Create a 2D array with None as placeholders
        grid = [[None for _ in range(len(unique_y))] for _ in range(len(unique_x))]
        
        # Populate the 2D array with the sorted points
        for x, y in sorted_centers:
            i = unique_x.index(x)
            j = unique_y.index(y)
            for pin in pins:
                if pin.center == (x,y):
                    grid[i][j] = pin.pin_id
        return grid 
        """
        sorted_pins = [pin for center in sorted_centers for pin in pins if pin.center == center]
        return sorted_pins
    
    def group_cells_in_pins(self, cells_to_pins):
        """
        This function groups the cells that go in fuel pins together to create fuel pins.
        """
        fuel_pin_groups = {}
        for cell in cells_to_pins:
            if cell.cell_name not in fuel_pin_groups.keys():
                fuel_pin_groups[cell.cell_name] = []
            fuel_pin_groups[cell.cell_name].append(cell)
        pins = []
        pin_number = 1
        for pin_name, pin_list in fuel_pin_groups.items():
            print(f"Processing pin : {pin_name}")
            materials = []
            surfaceIds = []
            list_radii = []
            for cell in pin_list:
                #print(f"Processing cell {cell.cell_name}")
                print(f"Cell surfaces are : {cell.surfaces}")
                materials.append(cell.material_number)
                surfaceIds = cell.surfaces
                pin_surfaces = self.getSurfacesFromIds(surfaceIds)
                pin_radii,center = self.getRadiiandCenterFromSurfaces(pin_surfaces)
                for rad in pin_radii:
                    if rad not in list_radii:
                        list_radii.append(rad)
            print(f"radii are {list_radii}")
            pins.append(S2_pin(pin_number, list(materials), list(list_radii), center))
            pin_number += 1
        return pins
    
    def mergePins(self, fuel_pins):
        """
        This function merges the fuel pins together to create a single pin per material aka fuel type.
        """
        indep_mat_list=[]
        indep_pins = []
        correspondance_dict = {}
        for pin in fuel_pins:
            print(f"Fuel pin {pin.pin_id} with materials {pin.materials} and radii {pin.radii}")
            if pin.central_hole:
                fuel_mat_idx = 1
            else:
                fuel_mat_idx = 0
            if pin.materials[fuel_mat_idx] not in indep_mat_list:
                indep_mat_list.append(pin.materials[fuel_mat_idx])
                indep_pins.append(pin)
                correspondance_dict[pin.materials[fuel_mat_idx]] = [pin.pin_id]
            else:
                correspondance_dict[pin.materials[fuel_mat_idx]].append(pin.pin_id)
        return indep_pins, correspondance_dict

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
    
    def getSurfacesFromIds(self, surface_ids):
        """
        This function retrieves the surfaces from the cells of the fuel pin.
        """
        surfaces = []
        for surf in surface_ids:
            for surface in self.surface_cards:
                if surface.surface_id == np.abs(int(surf)):
                    surfaces.append(surface)
        return surfaces
    
    
    def getRadiiandCenterFromSurfaces(self, surfaces):
        """
        This function retrieves the radii and center from the surfaces defining a cylinder.
        """
        radii = []
        for surface in surfaces:
            if "Cylinder" in surface.surface_type:
                radii.append(surface.radius)
                center = surface.center # assuming concentric cylinders
                print(f"Center of the cylinder is {center}")
        return radii, center


class S2_lattice:
    def __init__(self, name, universe, type, nx, ny, pitch, lattice):
        """
        Definition of Serpent2 lattice object. Case I from https://serpent.vtt.fi/mediawiki/index.php/Input_syntax_manual#lat_.28regular_lattice_definition.29
        Assuming finite 2D lattice in (x,y) plane and infinite in z direction.
        Attributes are :
        lattice name, universe, type, (nx,ny) = number of elements in x,y direction, pitch, data specifying the lattice
        data is pre-ordered with elements in 1) increasing x, 2) increasing direction to fill the lattice universe.
        """
        self.lattice_name = name
        self.universe = universe
        self.lattice_type = type
        self.nx = nx
        self.ny = ny
        self.pitch = pitch
        self.lattice = lattice
        print(f"Processing lattice {self.lattice_name}")
        #check conformity of lattice data
        if len(self.lattice) != self.nx*self.ny:
            print("Invalid lattice data : does not match nx*ny")
        return
class S2_cell:
    def __init__(self, name, universe_nb, material, surfaces_id_list):
        """
        Definition of Serpent2 cell object.
        Attributes are :
        cell type, cell name, universe number, material, surfaces list
        """
        #self.cell_type = cell_type
        self.cell_name = name
        self.universe_nb = universe_nb
        self.material = material
        self.surface_ids = surfaces_id_list
        
        print(f"Processing cell name : {self.cell_name}")
    def setSurfaces(self, surfaces):
        self.surfaces = surfaces
    
class S2_pin:
    def __init__(self, pin_id, materials, radii, center):
        """
        initializing fuel pin object.
        attributes are pin_id = universe number, materials = list of materials, radii = list of radii defining the fuel pin, pin center
        """
        self.central_hole = False # boolean to check if the fuel pin has a central hole, by default it does not as the case treated is a BWR pin lattice.
        self.pin_id = pin_id
        self.materials = materials
        self.radii = radii
        self.center = center 
        print(f"Processing fuel pin {self.pin_id}")
        return
class S2_surface:
    def __init__(self, surface_group, surface_id, type, surfaceData):
        """
        Definition of Serpent2 surface object.
        Attributes are :
        Not needed for Serpent2 surface definition : surface group attribute,
        Needed to specify Serpent2 surface : surface id, type of surface, parameters defining the surface
        """
        self.surface_group = surface_group
        self.surface_id = surface_id
        self.surface_type = type #type can be px, py, pz, sph, cylz or cyl, cylx, cyly, sqc, cube, cuboid, cross.
        if "Cylinder" in self.surface_type:
            #print(f"surface equation is = {surfaceEq}")
            self.center = (float(surfaceData[0]), float(surfaceData[1]))
            self.radius = float(surfaceData[2])
        else:
            self.surface_data = surfaceData
            
        #self.surface_parameters = parameters # the parameter array definies the surface according to mathematical definition in "Serpent user manual"
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