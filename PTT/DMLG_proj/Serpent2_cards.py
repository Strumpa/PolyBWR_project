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
    def __init__(self, name, lattice_type, input_data_type, input_case, mode, printlvl=0):
        """
        name (str) : name of the S2_case to be output
        input_data = input data parsed from MCNP (or Dragon) input file to create equivalent Serpent2 case. Either MCNP_case object or Dragon5_case object as created from DMLG_Interface.
        input_data_type = type of input data (MCNP or Dragon)
        lattice_type corresponds to the type of lattice described in the MCNP case
        """
        self.lattice_type = lattice_type

        self.name = name
        self.input_data_type = input_data_type
        self.mode = mode
        self.printlvl = printlvl

        self.pins = []
        self.cell_cards = []
        self.surface_cards = []
        self.material_cards = []
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
        
        # cards that require conversion
        input_cell_cards = []
        self.input_surface_cards = []

        # creation of serpent2 specific cards
        fuel_pins = []
        cells_to_pins = []
        lattice_cards = []

        # material cards require minor modification to be converted
        material_cards = []

        material_density_dict = {}


        for cell in input_case.cell_cards:
            if "pin" in cell.cell_name:
                cells_to_pins.append(cell)
                material_density_dict[cell.material_number] = cell.material_density
            else:
                input_cell_cards.append(S2_cell(cell.cell_name, 0, cell.material_number, cell.surface_ids))
                material_density_dict[cell.material_number] = cell.material_density
        for surface in input_case.surface_cards:
            self.input_surface_cards.append(S2_surface(surface.surfaces_group, surface.surface_number, surface.surface_type, surface.surface_data))

        for material in input_case.material_cards:
            if "t" in material.material_name:
                material_cards.append(S2_material(material.material_name, 0, [], material.iso_codes, material.therm_lib))
            else:
                material_number = material.material_name.split("m")[1]
                if material_number[0]=='0':
                    material_number = material_number[1]
                if self.printlvl > 1:
                    print(f"Processing material {material.material_name} with density {material_density_dict[int(material_number)]}")
                material_cards.append(S2_material(material.material_name, material_density_dict[int(material_number)], [], material.iso_codes, material.iso_densities))


        # Analyse cells to create fuel pins
        fuel_pins = self.group_cells_in_pins(cells_to_pins) # create Serpent2 fuel pins (obj) from MCNP cells,
        merged_pins, mat_to_pins_dict = self.mergePins(fuel_pins) # merge fuel pins with same material to create a single pin per material
    
        if self.printlvl > 1:
            print(f"Material to correpsonding dictionnary before merge is {mat_to_pins_dict}")
            for pin in merged_pins:
                print(f"Fuel pin {pin.pin_id} with materials {pin.materials} and radii {pin.radii}")
        
        merged_pins = self.renumberMergedPins(merged_pins)

        for pin in merged_pins:
            print(f"$$$$ Fuel pin {pin.pin_id} with materials {pin.materials} and radii {pin.radii}")
        # Analyse cartesian geometry to create the bounding surfaces and boxes  
        cell_cards, surface_cards = self.AnalyseMCNPGeometry(input_cell_cards)
        if self.printlvl > 1:
            for cell in cell_cards:
                print("$$$ Geometry analysis")
                print(f"Cell name : {cell.cell_name} with universe number {cell.universe_nb} and material {cell.material}")
                print(f"Cell surfaces are : {cell.surfaces}")



        # Combine info on fuel pins and surfaces to create the pin lattice.
        fuel_lattice_card = self.createLatticeFromPins(merged_pins)
        if self.printlvl > 1:
            print(f"Fuel lattice is {fuel_lattice_card}")

        return cell_cards, merged_pins, surface_cards, material_cards, lattice_cards
    

    def getCellDataFromSurfaces(self, surfaces):
        """
        Analyse surfaces to obtain information from geometry.
        """
        radii = []
        planes_x = []
        planes_y = []
        planes_z = []
        planes_ = []
        center = None
        for surface in surfaces:
            print(surface.surface_id)
            if "Cylinder" in surface.surface_type:
                radii.append(surface.radius)
                center = surface.center
            elif "px" == surface.surface_type:
                planes_x.append(float(surface.surface_data[0]))
            elif "py" == surface.surface_type:
                planes_y.append(float(surface.surface_data[0]))
            elif "pz" == surface.surface_type:
                planes_z.append(float(surface.surface_data[0]))
            elif "p" == surface.surface_type:
                planes_=[float(surface.surface_data[0]), float(surface.surface_data[1]), float(surface.surface_data[2]), float(surface.surface_data[3])]

        if radii and planes_x and planes_y and planes_z:
            if self.checkInclusion((center[0], center[1], radii[0]), (planes_x[0], planes_x[1], planes_y[0], planes_y[1], planes_z[0], planes_z[1])):
                print("Cylinder is included in the rectangular/square cell \n")
                mesh_data = (center[0], center[1], radii[0], planes_x, planes_y, planes_z)
                return "cylinder in cell", mesh_data

            elif self.checkIntersection((center[0], center[1], radii[0]), (planes_x[0], planes_x[1], planes_y[0], planes_y[1], planes_z[0], planes_z[1])):
        
                print("Cylinder intersects with the rectangular/square cell \n")
                mesh_data = (center[0], center[1], radii[0], planes_x, planes_y, planes_z)
                return "truncated cylinder", mesh_data
            else:
                print("Cylinder is outside the rectangular/square cell \n")
                mesh_data = (center[0], center[1], radii[0], planes_x, planes_y, planes_z)
                return "union or complement of cylinder and cell", mesh_data
        elif radii and planes_x:
            if len(planes_x)==1:
                print("Cylinder truncated by x plane \n")
                mesh_data = (center[0], center[1], radii[0], planes_x, planes_z)
                return "cylinder truncated by x plane", mesh_data
        
        elif planes_x and planes_y and planes_z and not radii:
            print(f"surface number is {surface.surface_id} in rectangle")
            if len(planes_x) == 2 and len(planes_y) == 2 and len(planes_z) == 2:
                mesh_data = (planes_x, planes_y, planes_z)
                return "rectangle", mesh_data
            
            
    
    def checkIntersection(self, cylinder_data, planes):
        """
        check if a cylinder intersects with a cuboid
        """

        center = cylinder_data[0], cylinder_data[1]
        radius = cylinder_data[2]
        x0, x1, y0, y1, z0, z1 = planes

        # center is included in the rectangle (cuboid projection on xy plane)
        if center[0] >= x0 and center[0] <= x1 and center[1] >= y0 and center[1] <= y1: 
            print("center is included in the rectangle")
            return True
        # center is not included in the rectangle but one of the rectangle's edges has a point in the cylinder
        elif min(abs(center[0]-x0), abs(center[0]-x1)) <= radius and min(abs(center[1]-y0), abs(center[1]-y1)) <= radius:
            print(f"abs(center[0]-x0) = {abs(center[0]-x0)}")
            print(f"abs(center[0]-x1) = {abs(center[0]-x1)}")
            print(f"abs(center[1]-y0) = {abs(center[1]-y0)}")
            print(f"abs(center[1]-y1) = {abs(center[1]-y1)}")
            print(f"radius = {radius}")
            return True
        else:
            return False
        
    def checkInclusion(self, cylinder_data, planes):
        """
        check inclusion of a cylinder in a cuboid
        """
        center = cylinder_data[0], cylinder_data[1]
        radius = cylinder_data[2]
        x0, x1, y0, y1, z0, z1 = planes

        # center is included in the rectangle (cuboid projection on xy plane)
        if center[0] >= x0 and center[0] <= x1 and center[1] >= y0 and center[1] <= y1:
            if center[0]+radius <= x1 and center[0]-radius >= x0 and center[1]+radius <= y1 and center[1]-radius >= y0:
                return True
            else: 
                return False
        else:
            return False



    def AnalyseMCNPGeometry(self, cell_cards):
        """
        Combine information about cells and surfaces to create bounding surfaces and boxes.
        """
        for cell in cell_cards:
            #print(cell.surface_ids)
            cell.setSurfaces(self.getSurfacesFromIds(cell.surface_ids, self.input_surface_cards))
            surfaces = cell.getSurfaces()
            if self.printlvl > 0:
                print(f"$$ New cell of name {cell.cell_name} is bounded by : \n") 
                for surf in surfaces:
                    print(f"surfaces of type {surf.surface_type}") 
            Cell_type, mesh_data = self.getCellDataFromSurfaces(surfaces)


            #print(f"Surface {surf.surface_id} of type {surf.surface_type} in group {surf.surface_group}")
        # Keep working on this function to create bounding surfaces and boxes.
        surface_cards = self.input_surface_cards
        return cell_cards, surface_cards



### ------------------------------------ Lattice handling functions

    def createLatticeFromPins(self, pins, type="square", diagonal_symmetry=True, box = "none"):
        """
        This function creates a lattice card from the fuel pins.
        assuming square BWR lattice.
        Need to figure out box indices.
        """
        for pin in pins:
            print(f"$ Creating lattice : Fuel pin {pin.pin_id} with materials {pin.materials} and radii {pin.radii}")
        pins_lat = self.sortPinsOnLattice(pins)
        pitch = pins[0].pitch
        print(len(pins_lat))

        lattice = S2_lattice("fuel", 0, type, 10, 10, pitch, pins_lat)
        lattice.plot_lattice()
        return lattice


    def sortPinsOnLattice(self, pins):
        """
        Sort pins according to increasing x postion and increasing y position
        """
        centers = []
        for pin in pins:
            centers.append(pin.center)
        sorted_centers = sorted(centers, key=lambda center: (center[0],center[1]))
        
        
        # Extract unique x and y coordinates
        unique_x = sorted(set(point[0] for point in sorted_centers))
        unique_y = sorted(set(point[1] for point in sorted_centers))
        
        # Create a 2D array with 0 as placeholders
        grid = np.array([[0 for _ in range(len(unique_y))] for _ in range(len(unique_x))])
        
        # Populate the 2D array with the sorted points
        for x, y in sorted_centers:
            i = unique_x.index(x)
            j = unique_y.index(y)
            for pin in pins:
                if pin.center == (x,y):
                    grid[i][j] = pin.pin_id
        #self.plot_lattice_test(grid, saveFlag=False)
        transposed_grid = grid[::-1,::-1].T

        #self.plot_lattice_test(transposed_grid, saveFlag=False)
        full_grid = np.zeros((len(transposed_grid), len(grid)))
        for i in range(len(grid)):
            for j in range(len(grid)):
                if i+j == len(grid)-1:
                    full_grid[i][j] = int(grid[i][j])
                else:
                    full_grid[i][j] = int(transposed_grid[i][j]+grid[i][j])

        return full_grid
    


### ------------------------------------- Fuel pin handling functions

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
            if self.printlvl > 5:
                print(f"Processing pin : {pin_name}")
            materials = []
            surfaceIds = []
            for cell in pin_list:
                if self.printlvl > 5:
                    print(f"Processing cell {cell.cell_name}")
                    print(f"Cell surface identifiers are : {cell.surface_ids}")
                materials.append(cell.material_number)
                surfaceIds.extend(cell.surface_ids)
            pin_surfaces = self.getSurfacesFromIds(surfaceIds, self.input_surface_cards)
            pin_radii, center, pitch, height, isOnSymmetryAxis  = self.getPinDataFromSurfaces(pin_surfaces)

            if self.printlvl > 5:
                print(f"materials are {materials}")
                print(f"radii are {pin_radii}")
                print(f"pin surface Ids are {surfaceIds}")
            pins.append(S2_pin(pin_number, list(materials), pin_radii, center, pitch, height, isOnSymmetryAxis))
            pin_number += 1
        return pins

    
    def getPinDataFromSurfaces(self, surfaces):
        """
        This function retrieves the radii and center from the surfaces defining a cylinder.
        """
        radii = []
        planes_x = []
        planes_y = []
        planes_z = []
        planes_ = []
        for surface in surfaces:
            if "Cylinder" in surface.surface_type:
                radii.append(surface.radius)
                center = surface.center # assuming concentric cylinders
                if self.printlvl > 5:
                    print(f"Center of the cylinder is {center}")
            elif "px" == surface.surface_type:
                planes_x.append(float(surface.surface_data[0]))
                if self.printlvl > 5:
                    print(f"surface data is = {surface.surface_data}")
            elif "py" == surface.surface_type:
                planes_y.append(float(surface.surface_data[0]))
                if self.printlvl > 5:
                    print(f"surface data is = {surface.surface_data}")   
            elif "pz" == surface.surface_type:
                if float(surface.surface_data[0]) not in planes_z:
                    planes_z.append(float(surface.surface_data[0]))
                    if self.printlvl > 5:
                        print(f"surface data is = {surface.surface_data}")
            elif "p" == surface.surface_type: 
                planes_=[float(surface.surface_data[0]), float(surface.surface_data[1]), float(surface.surface_data[2]), float(surface.surface_data[3])]
                if self.printlvl > 5:
                        print(f"surface data is = {surface.surface_data}")
        print(f"pincell planes are : {planes_x, planes_y, planes_z, planes_}")
        if planes_x:
            # check that there are only two planes
            if len(planes_x) == 2:
                pitch_x = abs(planes_x[1]-planes_x[0])
            else:
                print(f"Invalid number of planes in x direction for {self.lattice_type} pincell definition")
        if planes_y:
            # check that there are only two planes
            if len(planes_y) == 2:
                pitch_y = abs(planes_y[1]-planes_y[0])
            else:
                print(f"Invalid number of planes in y direction for {self.lattice_type} pincell definition")

        print(planes_z)
        if planes_z:
            # check that there are only two planes
            if len(planes_z) == 2:
                height = abs(planes_z[1]-planes_z[0])
            else:
                print(f"Invalid number of planes in z direction for {self.lattice_type} pincell definition")
        if planes_:
            isOnSymmetryAxis = True
        else:
            isOnSymmetryAxis = False
        print(pitch_x, pitch_y)
        if np.abs(pitch_x-pitch_y)<=1e-8 :
            pitch = pitch_x
            print(f"Pitch is {pitch_x} cm")
        else: 
            print("Pitch in x and y directions are different, pincell is not a square")

        return radii, center, pitch, height, isOnSymmetryAxis
    


    def mergePins(self, fuel_pins):
        """
        This function merges the fuel pins together to create a single pin per material aka fuel type.
        """
        print("$$$ in mergePins")
        indep_mat_list=[]
        indep_pins = []
        correspondance_dict = {}
        print(f"Fuelpins length is {len(fuel_pins)}")
        i=0
        for pin in fuel_pins:
            if self.printlvl > 5:
                print(f"Fuel pin {pin.pin_id} with materials {pin.materials} and radii {pin.radii}, i={i}")
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
                pin.setPinId(correspondance_dict[pin.materials[fuel_mat_idx]][0])
                indep_pins.append(pin)
            i+=1
        return indep_pins, correspondance_dict
    


    def renumberMergedPins(self, mergedPins):
        """
        This functions renumerates the merged fuel pins in order to number fuel rods in conscutive pin_ids.
        """
        print("$$$ in renumberMergedPins")
        print(f"mergedPins length is {len(mergedPins)}")
        print(f"mergedPins are {mergedPins}")
        pin_numbers=[]
        renumbered_ids= []
        for pin in mergedPins:
            if pin.pin_id not in pin_numbers:
                pin_numbers.append(pin.pin_id)
                renumbered_ids.append(pin.pin_id)
                pin_numbers.sort()
                renumbered_ids.sort()
                if len(pin_numbers)>1:
                    if pin_numbers[-1] > pin_numbers[-2]:
                        renumbered_ids[-1] = renumbered_ids[-2]+1
        
            #pin_numbers.append(pin.pin_id)
            print(f"pin_numbers are {renumbered_ids}")
            for i in range(len(renumbered_ids)):
                if pin.pin_id == pin_numbers[i]:
                    pin.setPinId(renumbered_ids[i])

        return mergedPins
    
    
### ------------------------------------- Surface handling functions

    def getSurfacesFromIds(self, surface_ids, surface_cards):
        """
        This function retrieves the surfaces from the cells of the fuel pin.
        """
        surfaces = []
        for surf in surface_ids:
            for surface in surface_cards:
                if surface.surface_id == np.abs(int(surf)):
                    if surface not in surfaces:
                        surfaces.append(surface)
        return surfaces
    


### ------------------------------------- Output generating functions

    def print_cells_info(self, output_path, output_name):
        """
        this function writes to a text file the cell cards information for the Serpent2 case.
        """

        if not os.path.exists(output_path):
            os.mkdir(output_path)
        output = open(output_path+output_name, "w")

        for pin in self.fuel_pins:
            if self.printlvl > 5:
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
            if self.printlvl > 5:
                print(f"Cell name : {cell.cell_name} with universe number {cell.universe_nb} and material {cell.material}")
                print(f"Cell surfaces are : {cell.surfaces}")
            output.write(f"{cell.cell_type} {cell.cell_name} {cell.universe_nb} {cell.material} {cell.surfaces}\n")



    def plot_lattice_test(self, grid, saveFlag=False):
        """
        write a function that plots the lattice with the fuel pins. Each pin has an associated color,
        the pins are represented by a square. they have to be ordered in increating x and y coordinate.
        """
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        fig, ax = plt.subplots()
        for i in range(len(grid)):
            for j in range(len(grid[i])):
                if grid[i][j] != 0:
                    rect = patches.Rectangle((i,j), 1, 1, linewidth=1, edgecolor='r', facecolor='none')
                    ax.add_patch(rect)
                    ax.text(i+0.5, j+0.5, str(grid[i][j]), ha='center', va='center')
        ax.set_xlim(0, len(grid))
        ax.set_ylim(0, len(grid[0]))
        ax.set_aspect('equal')
        if saveFlag:
            plt.savefig(f"lattice_.png")
        plt.show()
        plt.close()
        return


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
        self.lattice_grid = lattice
        print(f"Processing lattice {self.lattice_name}")
        #check conformity of lattice data
        if len(self.lattice_grid[0]*self.lattice_grid[1]) != self.nx*self.ny:
            print("Invalid lattice data : does not match nx*ny")
        return
    def plot_lattice(self):
        """
        write a function that plots the lattice with the fuel pins. Each pin has an associated color,
        the pins are represented by a square. they have to be ordered in increating x and y coordinate.
        """
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        fig, ax = plt.subplots()
        for i in range(len(self.lattice_grid)):
            for j in range(len(self.lattice_grid[i])):
                if self.lattice_grid[i][j] != 0:
                    rect = patches.Rectangle((i,j), 1, 1, linewidth=1, edgecolor='r', facecolor='none')
                    ax.add_patch(rect)
                    ax.text(i+0.5, j+0.5, int(self.lattice_grid[i][j]), ha='center', va='center')
        ax.set_xlim(0, len(self.lattice_grid))
        ax.set_ylim(0, len(self.lattice_grid[0]))
        ax.set_aspect('equal')
        plt.savefig(f"lattice_{self.lattice_name}.png")
        plt.show()
        plt.close()

        return
        
class S2_cell:
    def __init__(self, name, universe_nb, material, surfaces_id_list):
        """
        Definition of Serpent2 cell object.
        Attributes are :
        cell name, universe number, material, surfaces list
        """
        self.cell_name = name
        self.universe_nb = universe_nb
        self.material = material
        self.surface_ids = surfaces_id_list
        
        #print(f"Processing cell name : {self.cell_name}")
    def setSurfaces(self, surfaces):
        self.surfaces = surfaces
        return
    def getSurfaces(self):
        return self.surfaces
    
class S2_pin:
    def __init__(self, pin_id, materials, radii, center, pitch, height, isOnSymmetryAxis=False):
        """
        initializing fuel pin object.
        attributes are pin_id = universe number, materials = list of materials, radii = list of radii defining the fuel pin, pin center
        """
        self.central_hole = False # boolean to check if the fuel pin has a central hole, by default it does not as the case treated is a BWR pin lattice.
        self.pin_id = pin_id
        self.materials = materials
        self.radii = radii
        self.center = center
        self.pitch = pitch
        self.height = height
        self.isOnSymmetryAxis = isOnSymmetryAxis 
        #print(f"Processing fuel pin {self.pin_id}")
        return
    def setPinId(self, pin_id):
        self.pin_id = pin_id
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
        #print(f"Processing surface {self.surface_id} of type {self.surface_type}")
        return
    

    def setS2Surface_type(self, S2surf_type):
        """
        Intialized cell type attribute. For now only used for sqc, cuboid and cyl cells.
        """
        self.surf_type = S2surf_type
        return

    def setS2Surface_mesh(self, mesh_data):
        """
        function used to set the mesh data to be written to Serpent2 surface definition
        """
        if self.cell_type == "sqc":
            """
            2D square, assumed inifnite in the z direction
            mesh data must take the form of a tuple (x0, y0, d) where x0 and y0 are the coordinates of the square's center and d is the half width of the square.
            """
            self.x0 = mesh_data[0] # square's center x coordinate
            self.y0 = mesh_data[1] # square's center y coordinate
            self.d = mesh_data[2] # square's half width
        elif self.cell_type == "cuboid":
            """
            3D paralleliped
            mesh data must take the form (x0, x1, y0, y1, z0, z1) where x=x0, x=x1, y=y0, y=y1, z=z0, z=z1 are the equations definiting the planes bounding the 3D cuboid
            """
            self.x0 = mesh_data[0]
            self.x1 = mesh_data[1]
            self.y0 = mesh_data[2]
            self.y1 = mesh_data[3]
            self.z0 = mesh_data[4]
            self.z1 = mesh_data[5]
        elif self.cell_type == "cyl" or self.cell_type == "cylz":
            """
            infinite cylinder centered in (x0,y0) of radius r, parallel to the z axis, truncated/bounded by 2 planes (z0, z1).
            """
            self.center = (mesh_data[0],mesh_data[1])
            self.radius = mesh_data[2]
            self.z0 = mesh_data[3]
            self.z1 = mesh_data[4]
        return
    
    def setIsPincellSurface(self, isPincellSurface):
        self.isPinCellSurface = isPincellSurface
        return

class S2_material:
    def __init__(self, name, dens, options, isotopes, data):
        """
        Definition of Serpent2 material object.
        Attributes are : material name, density, options, isotopes, data
        unit is an attribute that is deduced from the data
        """
        self.material_name = name
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
            print(f"For material {self.material_name} data and isotopes do not have the same length")
            print("Data and isotopes do not have the same length")
            return
        
        return
    def setMaterialDensity(self, dens):
        self.density = dens
        return