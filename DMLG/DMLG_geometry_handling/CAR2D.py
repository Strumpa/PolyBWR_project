# Python3 class to represent a DRAGON5 CAR2D geometry object
# Author: R. Guasch
# Date : 24 July 2024
# Dependencies : GEOM.py
#

from GEOM import GEO
from CARCEL import CARCEL
import numpy as np
from matplotlib import pyplot as plt


class CAR2D(GEO):
    def __init__(self, name, level, nx, ny, nz, meshx, meshy, meshz):
        """
        Definition of a 2D cartesian geometry according to the Dragon5 definition
        :param name: name of the geometry
        :param lx: number of cells in x direction
        :param ly: number of cells in y direction
        :param meshx: mesh in x direction : array of nx+1 values, bounds to the cells
        :param meshy: mesh in y direction : array of ny+1 values, bounds to the cells
        :param level: level of the geometry in the hierarchy
        """
        super().__init__(name, level)
        self.type = "CAR2D"
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.meshx = meshx
        self.meshy = meshy
        self.meshz = meshz

        self.xmin = self.meshx[0]
        self.xmax = self.meshx[-1]
        self.ymin = self.meshy[0]
        self.ymax = self.meshy[-1]
        self.zmin = self.meshz[0]
        self.zmax = self.meshz[-1]

        self.origin = [self.xmin, self.ymin] # Origin of the geometry, lower left corner, assumung z invariance

        self.number_of_regions = self.nx * self.ny

        if self.nx != len(self.meshx) - 1:
            print("$$- CAR2D: Error in meshx definition")
        if self.ny != len(self.meshy) - 1:
            print("$$- CAR2D: Error in meshy definition")

        if self.level == 1 :
            print("$$- CAR2D: Initializing boundary conditions to reflective for level 1")
            self.BoundaryConditions = {"X- ":"RELF", "X+ ":"RELF", "Y- ":"RELF", "Y+ ":"RELF"}
            
        self.number_regions()
        self.number_bounding_surfaces()
        self.number_all_surfaces()
        print(f"$$- CAR2D: All surfaces: {self.all_surfaces}")
        self.associateSurfacesToRegions()

    def number_regions(self):
        """
        In a CAR2D geometry, the regions are numbered from the lower left corner to the upper right corner, from 1 to lx*ly
        """
        self.region_numbering = []
        self.volumes_of_regions = []
        for i in range(1,self.nx+1):
            for j in range(1,self.ny+1):
                for k in range(1,self.nz+1):
                    l = i+self.nx*(j-1+self.ny*(k-1))
                    self.region_numbering.append([i, j, k, l]) # Equation 3.1 from NXT guide
                    self.volumes_of_regions.append([l,(self.meshx[i]-self.meshx[i-1])*(self.meshy[j]-self.meshy[j-1])*(self.meshz[k]-self.meshz[k-1])]) # Equation 3.2 from NXT guide
        self.region_numbering = np.array(self.region_numbering)
        self.volumes_of_regions = np.array(self.volumes_of_regions)
        self.filled_regions = [False] * self.number_of_regions
        return
    
    def number_bounding_surfaces(self):
        """
        Order for surface numbering in CAR2D geometry :
            # 1 : surfaces at xmin bound
            # 2 : surfaces at xmax bound
            # 3 : surfaces at ymin bound
            # 4 : surfaces at ymax bound
        These represent all surfaces bounding the geometry.
        """
        # at the xmin/max bounds :
        self.bounding_surfaces =[]
        self.nxs = self.ny*self.nz # number of surfaces in x direction
        self.surface_ids_xmin = []
        self.surface_ids_xmax = []

        for j in range(1,self.ny+1):
            for k in range(1,self.nz+1):
                l = -(j+self.ny*(k-1))
                self.surface_ids_xmin.append(l)
                self.surface_ids_xmax.append(l-self.nxs)
                self.bounding_surfaces.append([self.xmin, j, k, l,(self.meshy[j]-self.meshy[j-1])*(self.meshz[k]-self.meshz[k-1])])
                self.bounding_surfaces.append([self.xmax, j, k, l-self.nxs,(self.meshy[j]-self.meshy[j-1])*(self.meshz[k]-self.meshz[k-1])])
        
        print(f"$$- CAR2D: Surface identifiers at xmin bound: {self.surface_ids_xmin}")
        print(f"$$- CAR2D: Surface identifiers at xmax bound: {self.surface_ids_xmax}")
        
        # at the ymin/max bounds :
        self.nys = self.nx*self.nz # number of surfaces in y direction
        self.surface_ids_ymin = []
        self.surface_ids_ymax = []

        for i in range(1,self.nx+1):
            for k in range(1,self.nz+1):
                l = -(k+self.nz*(i-1))-2*self.nxs
                self.surface_ids_ymin.append(l)
                self.surface_ids_ymax.append(l-self.nys)
                self.bounding_surfaces.append([i, self.ymin, k, l,(self.meshx[i]-self.meshx[i-1])*(self.meshz[k]-self.meshz[k-1])])
                self.bounding_surfaces.append([i, self.ymax, k, l-self.nys,(self.meshx[i]-self.meshx[i-1])*(self.meshz[k]-self.meshz[k-1])])

        print(f"$$- CAR2D: Surface identifiers at ymin bound: {self.surface_ids_ymin}")
        print(f"$$- CAR2D: Surface identifiers at ymax bound: {self.surface_ids_ymax}")

        self.bounding_surfaces = np.array(self.bounding_surfaces)
        self.local_surface_numbering = np.array(self.bounding_surfaces)
        self.number_of_bounding_surfaces = len(self.bounding_surfaces)
        return
    
    def number_all_surfaces(self):
        """
        Number all surfaces in the geometry, 
        These include the bounding surfaces, as well as the inner surfaces (interfaces between regions) and potential bounds to the sub-geometries
        The surfaces are numbered following the same order as the bounding surfaces :
        # 1 : surfaces at x bounds, then x1, x2, x3, ... xlx,
        # 2 : surfaces at y bounds, then y1, y2, y3, ... yly,
        """

        self.all_surfaces = []
        self.all_surface_ids = []
        self.number_of_x_oriented_surfaces = (self.nx+1)*(self.ny)
        self.number_of_y_oriented_surfaces = (self.ny+1)*(self.nx)
        print(f"$$- CAR2D: {self.nx},{self.ny}: Number of x oriented surfaces: {self.number_of_x_oriented_surfaces}")
        print(f"$$- CAR2D: {self.nx},{self.ny}: Number of y oriented surfaces: {self.number_of_y_oriented_surfaces}")

        # x bounds
        print(f"$$- CAR2D: nxs = {self.nxs}, nys = {self.nys}")
        for i in range(self.nx+1):
            for j in range(1,self.ny+1):
                for k in range(1,self.nz+1):
                    l = (j+self.ny*(k-1))+i*self.nxs
                    self.all_surfaces.append([i, j, k, l])
                    self.all_surface_ids.append(l)
        # y bounds
        # Identified an issue with this as this way of numering inner surfaces produces a non-unique numbering, 
        # x and y oriented surfaces touching the symmetry axis are numbered twice. Find a way to avoid this ?
        for j in range(self.ny+1):
            for i in range(1,self.nx+1):
                for k in range(1,self.nz+1):
                    l = (k+self.nz*(i-1))+j*self.nys + self.number_of_x_oriented_surfaces
                    self.all_surfaces.append([i, j, k, l]) 
                    self.all_surface_ids.append(l)

        self.all_surfaces = np.array(self.all_surfaces)
        print(f"$$- CAR2D: Number of all surfaces: {len(self.all_surfaces)}")


        return
    
    def add_CARCEL(self, fuel_identifier, host_region, radii, meshx, meshy):
        """
        Add a CARCEL object to the geometry object
        :param fuel_identifier (str): fuel identifier of the CARCEL object
        :param host_region: region where the CARCEL object is added
        """

        # Sanity check : does the CARCEL fit in the region ?
        # get the bounds of the region
        position = self.getPositionFromRegion(host_region)
        region_bounds = self.getBoundsFromPosition(position[0], position[1], position[2])
        delta_x = region_bounds[1] - region_bounds[0]
        delta_y = region_bounds[3] - region_bounds[2]
        print(f"$$- CAR2D: Region {host_region} bounds: {region_bounds}")
        if max(meshx)-min(meshx) > delta_x or max(meshy)-min(meshy) > delta_y:
            print(f"$$- CAR2D: Error : CARCEL object {fuel_identifier} does not fit in region {host_region}")
            print(f"$$- CAR2D: CARCEL object {fuel_identifier} bounds: {min(meshx), max(meshx), min(meshy), max(meshy)}")
            print(f"$$- CAR2D: Region {host_region} bounds: {region_bounds}")
            print(f"Error encountered defining the geometry, please check the mesh and the CARCEL object bounds")
            return
        # Check if the CARCEL object is already in the geometry
        if self.filled_regions[host_region-1]:
            print(f"$$- CAR2D: Error : region {host_region} already filled with a sub-geometry")
            return
        if host_region in self.region_numbering[:,3]:
            # Create a new CARCEL object to fill region
            print(f"$$- CAR2D: Adding CARCEL object {fuel_identifier} to region {host_region}")
            cell = CARCEL(name=fuel_identifier, level=self.level+1, nr=len(radii)-1, radii=radii, meshx=meshx, meshy=meshy)
            cell.setHostGeometry(self.name)
            cell.setHostRegion(host_region)
            position = self.getPositionFromRegion(host_region)
            print(f"$$- CAR2D: Position of the region: {position}")
            bounds = self.getBoundsFromPosition(position[0], position[1], position[2])
            for correspondance_reg_to_surf in self.region_to_surfaces:
                if correspondance_reg_to_surf[0] == host_region:
                    surfaces = correspondance_reg_to_surf[1]
                    print(f"$$- CAR2D: Surfaces bounding region {host_region}: {surfaces}")
                    cell.setHostBoundingSurfaces(surfaces)
                    break
            self.sub_geometries.append(cell)
            self.filled_regions[host_region-1] = True
            print(f"$$- CAR2D: Cell {cell.name} added to region {host_region} with bounds {bounds}")
        else:
            print(f"$$- CAR2D: Error : host region {host_region} not found in the geometry")

        return
        

    def add_cell(self, cell, host_region):
        """
        Add a cell to the geometry object
        :param cell: cell object to be added
        :param host_region: region where the cell is added
        """
        
        if cell.level != self.level+1:
            print("$$- CAR2D: Error : cell level, must be 1 above the host geometry level")
            return
        if host_region in self.region_numbering[:,3]:
            cell.setHostGeometry(self.name)
            # Check if the cell is already in the geometry
            if cell in self.sub_geometries:
                self.counter += 1
                print(f"$$- CAR2D: Cell {cell.name} already in the geometry")
                print(f"Attemping to create a new cell object with updated cell name")
                cell = cell.__class__(f"{cell.name}{self.counter}", level=cell.level, nr=len(cell.radii)-1, radii=cell.radii, meshx=cell.meshx, meshy=cell.meshy)
            else:
                print(f"$$- CAR2D: Cell {cell.name} not in the list of sub-geometries")
                self.counter = 0
            print(f"$$- CAR2D: Adding cell {cell.name} to region {host_region}")
            cell.setHostRegion(host_region)
            position = self.getPositionFromRegion(host_region)
            print(f"$$- CAR2D: Position of the region: {position}")
            bounds = self.getBoundsFromPosition(position[0], position[1], position[2])
            for correspondance_reg_to_surf in self.region_to_surfaces:
                if correspondance_reg_to_surf[0] == host_region:
                    surfaces = correspondance_reg_to_surf[1]
                    print(f"$$- CAR2D: Surfaces bounding region {host_region}: {surfaces}")
                    cell.setHostBoundingSurfaces(surfaces)
                    break
            self.sub_geometries.append(cell)
            print(f"appending cell {cell.name} to the list of sub-geometries")
            print(f"$$- CAR2D: Cell {cell.name} added to region {host_region} with bounds {bounds}")

        else:
            print(f"$$- CAR2D: Error : host region {host_region} not found in the geometry")
        return
    

    def getSurfacefromNumbering(self, surface_number):
        """
        Get the surface from its number
        :param surface_number: number of the surface
        :return: surface
        """
        for i in range(len(self.all_surfaces)):
            if self.all_surfaces[i][3] == surface_number:
                return self.all_surfaces[i]
        return

    def associateSurfacesToRegions(self):
        """
        Each region in the nth level geometry is associated to a set of bounding surfaces. These can be the outer surfaces of the nth level geom
        or the inner surfaces connecting cells. In the presence of a nested geometry, these surfaces can be used to build a connectivity map between the
        regions of the n^th level geometry and the cells of the (n+1)^th level geometry.
        """
        print(f"$$- CAR2D: Associating surfaces to regions")
        print(f"$- CAR2D: self.region_numbering: {self.region_numbering}")
        self.region_to_surfaces = [] # [0: region id, ids of its bounding surfaces: 1:xmin, 2:xmax, 3:ymin, 4:ymax]
        for region_info in self.region_numbering:
            region = region_info[3]
            i = region_info[0]
            j = region_info[1]
            k = region_info[2]
            bounds = self.getBoundsFromPosition(i, j, k)
            print(f"$$- CAR2D: Region {region} bounds: {bounds}")
            
            temp = [] # temporary list to store the surfaces bounding the region
            # Get the surfaces bounding the region
            for surface in self.all_surfaces:
                if (surface[0] == i and surface[1] == j and surface[2] == k): # Problem with the surface numbering, do this to prototype but might need to change
                    print(f"$$- CAR2D: Surface {surface} bounding region {region}")
                    if surface[3]<=self.number_of_x_oriented_surfaces:
                        matching_surfaces = [surface[3]-self.nxs, surface[3]]
                    else:
                        matching_surfaces = [surface[3]-self.nys, surface[3]]
                    temp.extend(matching_surfaces)
            self.region_to_surfaces.append([region, np.array(temp)])
        
        #self.region_to_surfaces = np.array(self.region_to_surfaces)
        print(f"$$- CAR2D: Region to surfaces: {self.region_to_surfaces}")
        return
    
    def getSurfaceGeometricalData(self, surface_number):
        """
        Get the geometrical data of the surface with the given number
        :param surface_number: number of the surface
        :return: geometrical data of the surface, combing the all_surfaces array and the mesh arrays
        """
        for i in range(len(self.all_surfaces)):
            #print(f"Total number of surfaces: {len(self.all_surfaces)}")
            #print("Total number of x oriented surfaces: ", self.number_of_x_oriented_surfaces)
            #print("Total number of y oriented surfaces: ", self.number_of_y_oriented_surfaces)
            if self.all_surfaces[i][3] == surface_number:
                surface = self.all_surfaces[i]
                print(f"$$- CAR2D: Surface {surface_number} found for surface {surface}")
                if surface[3] <= self.number_of_x_oriented_surfaces:
                    is_x_oriented = True
                    x = self.meshx[surface[0]]
                    y1 = self.meshy[surface[1]]
                    y2 = self.meshy[surface[1]-1]
                    print(f"$$- CAR2D: x-oriented surface {surface_number} with geometrical data: x = {x} from y1 = {y1} to y2 = {y2}")
                    return is_x_oriented, x, y1, y2
                else:
                    is_x_oriented = False
                    x1 = self.meshx[surface[0]-1]
                    x2 = self.meshx[surface[0]]
                    y = self.meshy[surface[1]]
                    print(f"$$- CAR2D: y-oriented surface {surface_number} with geometrical data: y = {y} from x1 = {x1} to x2 = {x2}")
                    return is_x_oriented, y, x1, x2


    def describeGeo(self):
        print(f"$$- CAR2D: Created geometry (GEO: CAR2D) object of name {self.name} and level {self.level}, with {self.nx} cells in x direction and {self.ny} cells in y direction")
        print(f"$$- CAR2D: Region numbering from lower left corner to upper right corner: {self.region_numbering}, {self.number_of_regions} regions in total")
        print(f"$$- CAR2D: Region volumes: {self.volumes_of_regions}")
        print(f"$$- CAR2D: Outer surfaces: {self.bounding_surfaces}")
        for i in range(len(self.bounding_surfaces)):
            if i % 2 == 0:
                print(f"$$- CAR2D: Surface i={self.bounding_surfaces[i][0]}, j={self.bounding_surfaces[i][1]}, k={self.bounding_surfaces[i][2]}, with label {self.bounding_surfaces[i][3]}, with area {self.bounding_surfaces[i][4]}")
            else:
                print(f"$$- CAR2D: Surface i={self.bounding_surfaces[i][0]}, j={self.bounding_surfaces[i][1]}, k={self.bounding_surfaces[i][2]}, with label {self.bounding_surfaces[i][3]}, with area {self.bounding_surfaces[i][4]}")  
        print(f"$$- CAR2D: mesh in x direction: {self.meshx}")
        print(f"$$- CAR2D: mesh in y direction: {self.meshy}")
        if self.level == 1:
            print(f"$$- CAR2D: Boundary conditions: {self.BoundaryConditions}")

        print("$$- CAR2D: Inner labelling of surfaces :")
        for i in range(len(self.all_surfaces)):
            if i < self.nxs*(self.nys+1):
                print(f"$$- CAR2D: x surface i={self.all_surfaces[i][0]}, j={self.all_surfaces[i][1]}, k={self.all_surfaces[i][2]}, with label {self.all_surfaces[i][3]}")
            else:
                print(f"$$- CAR2D: y surface i={self.all_surfaces[i][0]}, j={self.all_surfaces[i][1]}, k={self.all_surfaces[i][2]}, with label {self.all_surfaces[i][3]}")

        if self.sub_geometries:
            print("$$- CAR2D: Sub-geometries:")
            for i in range(len(self.sub_geometries)):
                print(f"$$- CAR2D: Creation of a level {self.sub_geometries[i].level} sub-geometry {self.sub_geometries[i].name} added to region {self.sub_geometries[i].host_region}")

        return
