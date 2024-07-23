# Python3 script to define the geoemtry object for the SYBILT prototype
# Author: R. Guasch
# Date : 18 July 2024
# Purpose : allow for geometry creation and manipulation, to be analysed by the SYBILT prototype
# Usage : import GEOM

import numpy as np
import matplotlib.pyplot as plt


class GEO:
    def __init__(self, name, level, host_region = None):
        self.name = name
        self.level = level
        self.sub_geometries = []

    def describeGeo(self):
        print(f"$$- GEO: Created geometry (GEO) object of name {self.name} and level {self.level}")

    def setHostRegion(self, host_region):
        self.host_region = host_region
        return
    
    def setBoundingSurfaces(self, surfaces):
        self.bounding_surfaces = surfaces
        return
    
    def setLevel(self, level):
        self.level = level
        return
    
    def buildConnectivityMap(self):
        """
        Build a connectivity map between the regions of the n^th level geometry and the cells of the (n+1)^th level geometry.
        """
        self.connectivity_map = []
        print(f"$$- GEO: Building connectivity map between regions and cells")
        print(f"$$- GEO: Region numbering: {self.region_numbering}")
        print(f"$$- GEO: Sub-geometries: {self.sub_geometries}")
        print(f"$$- GEO: level of sub-geometries: {self.sub_geometries[0].level}")
        if self.sub_geometries[0].sub_geometries:
            print(f"$$- GEO: Sub-geometries of sub-geometries: {self.sub_geometries[0].sub_geometries}")
        for i in range(len(self.region_numbering)):
            region = self.region_numbering[i][3]
            print(f"$$- GEO: Geometry {self.name} has bounding surfaces {self.bounding_surfaces}")
            tmp = []
            for cell in self.sub_geometries:
                if cell.host_region == region:
                    print(f"$$- GEO: Checking cell {cell.name} with host region {cell.host_region}")
                    tmp.append([self.level, region, cell.bounding_surfaces, cell.local_surface_numbering])
                    if cell.sub_geometries:
                        for sub_cell in cell.sub_geometries:
                            tmp.append([cell.level, cell.host_region, sub_cell.bounding_surfaces, sub_cell.local_surface_numbering])

            self.connectivity_map.append(tmp)

        print(f"$$- GEO: Connectivity map between regions and cells: {self.connectivity_map}")
        return
    
    def plotGeo(self):
        """
        Plot the geometry
        """
        fig, ax = plt.subplots()
        ax.set_xlim(self.xmin-0.5, self.xmax+0.5)
        ax.set_ylim(self.ymin-0.5, self.ymax+0.5)
        ax.set_aspect('equal')
        ax.set_title(f"Geometry {self.name}")
        # Plot the bounding surfaces
        print(f"$$- GEO: Plotting bounding surfaces: {self.bounding_surfaces}")

        for i in range(len(self.bounding_surfaces)):
            if i >=2:
                ax.plot([self.bounding_surfaces[0], self.bounding_surfaces[1]], [self.bounding_surfaces[i], self.bounding_surfaces[i]], 'r-')
            else:
                ax.plot([self.bounding_surfaces[i], self.bounding_surfaces[i]], [self.bounding_surfaces[2], self.bounding_surfaces[3]], 'r-')
        
        # Plot the inner surfaces
        print(f"$$- GEO: Plotting inner surfaces: {self.inner_surfaces}")
        for i in range(len(self.inner_surfaces)):
            circle = plt.Circle((self.inner_surfaces[i][2], self.inner_surfaces[i][3]), self.inner_surfaces[i][1], fill=False)
            ax.add_artist(circle)
        plt.show()
        plt.savefig(f"GEO_{self.name}.png")
        return
    




class CAR2D(GEO):
    def __init__(self, name, level, lx, ly, lz, meshx, meshy, meshz):
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
        self.lx = lx
        self.ly = ly
        self.lz = lz
        self.meshx = meshx
        self.meshy = meshy
        self.meshz = meshz

        self.xmin = self.meshx[0]
        self.xmax = self.meshx[-1]
        self.ymin = self.meshy[0]
        self.ymax = self.meshy[-1]
        self.zmin = self.meshz[0]
        self.zmax = self.meshz[-1]

        self.number_of_regions = self.lx * self.ly

        if self.lx != len(self.meshx) - 1:
            print("$$- CAR2D: Error in meshx definition")
        if self.ly != len(self.meshy) - 1:
            print("$$- CAR2D: Error in meshy definition")

        if self.level == 1 :
            print("$$- CAR2D: Initializing boundary conditions to reflective for level 1")
            self.BoundaryConditions = {"X- ":"RELF", "X+ ":"RELF", "Y- ":"RELF", "Y+ ":"RELF"}\
            
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
        for i in range(1,self.lx+1):
            for j in range(1,self.ly+1):
                for k in range(1,self.lz+1):
                    l = i+self.lx*(j-1+self.ly*(k-1))
                    self.region_numbering.append([i, j, k, l]) # Equation 3.1 from NXT guide
                    self.volumes_of_regions.append([l,(self.meshx[i]-self.meshx[i-1])*(self.meshy[j]-self.meshy[j-1])*(self.meshz[k]-self.meshz[k-1])]) # Equation 3.2 from NXT guide
        self.region_numbering = np.array(self.region_numbering)
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
        self.nxs = self.ly*self.lz # number of surfaces in x direction
        self.surface_ids_xmin = []
        self.surface_ids_xmax = []

        for j in range(1,self.ly+1):
            for k in range(1,self.lz+1):
                l = -(j+self.ly*(k-1))
                self.surface_ids_xmin.append(l)
                self.surface_ids_xmax.append(l-self.nxs)
                self.bounding_surfaces.append(["xmin", j, k, l,(self.meshy[j]-self.meshy[j-1])*(self.meshz[k]-self.meshz[k-1])])
                self.bounding_surfaces.append(["xmax", j, k, l-self.nxs,(self.meshy[j]-self.meshy[j-1])*(self.meshz[k]-self.meshz[k-1])])
        
        print(f"$$- CAR2D: Surface identifiers at xmin bound: {self.surface_ids_xmin}")
        print(f"$$- CAR2D: Surface identifiers at xmax bound: {self.surface_ids_xmax}")
        
        # at the ymin/max bounds :
        self.nys = self.lx*self.lz # number of surfaces in y direction
        self.surface_ids_ymin = []
        self.surface_ids_ymax = []

        for i in range(1,self.lx+1):
            for k in range(1,self.lz+1):
                l = -(k+self.lz*(i-1))-2*self.nxs
                self.surface_ids_ymin.append(l)
                self.surface_ids_ymax.append(l-self.nys)
                self.bounding_surfaces.append([i, "ymin", k, l,(self.meshx[i]-self.meshx[i-1])*(self.meshz[k]-self.meshz[k-1])])
                self.bounding_surfaces.append([i, "ymax", k, l-self.nys,(self.meshx[i]-self.meshx[i-1])*(self.meshz[k]-self.meshz[k-1])])

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
        self.number_of_x_oriented_surfaces = (self.lx+1)*(self.ly)
        self.number_of_y_oriented_surfaces = (self.ly+1)*(self.lx)
        print(f"$$- CAR2D: {self.lx},{self.ly}: Number of x oriented surfaces: {self.number_of_x_oriented_surfaces}")
        print(f"$$- CAR2D: {self.lx},{self.ly}: Number of y oriented surfaces: {self.number_of_y_oriented_surfaces}")

        # x bounds
        print(f"$$- CAR2D: nxs = {self.nxs}, nys = {self.nys}")
        for i in range(self.lx+1):
            for j in range(1,self.ly+1):
                for k in range(1,self.lz+1):
                    l = (j+self.ly*(k-1))+i*self.nxs
                    self.all_surfaces.append([i, j, k, l])
                    self.all_surface_ids.append(l)
        # y bounds
        # Identified an issue with this as this way of numering inner surfaces produces a non-unique numbering, 
        # x and y oriented surfaces touching the symmetry axis are numbered twice. Find a way to avoid this ?
        for j in range(self.ly+1):
            for i in range(1,self.lx+1):
                for k in range(1,self.lz+1):
                    l = (k+self.lz*(i-1))+j*self.nys + self.number_of_x_oriented_surfaces
                    self.all_surfaces.append([i, j, k, l]) 
                    self.all_surface_ids.append(l)

        self.all_surfaces = np.array(self.all_surfaces)
        print(f"$$- CAR2D: Number of all surfaces: {len(self.all_surfaces)}")


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
            print(f"$$- CAR2D: Adding cell {cell.name} to region {host_region}")
            cell.setHostRegion(host_region)
            position = self.getPosFromRegion(host_region)
            print(f"$$- CAR2D: Position of the region: {position}")
            bounds = self.getBoundsFromPos(position[0], position[1], position[2])
            for correspondance_reg_to_surf in self.region_to_surfaces:
                if correspondance_reg_to_surf[0] == host_region:
                    surfaces = correspondance_reg_to_surf[1]
                    print(f"$$- CAR2D: Surfaces bounding region {host_region}: {surfaces}")
                    cell.setBoundingSurfaces(surfaces)
                    break
            self.sub_geometries.append(cell)
            print(f"$$- CAR2D: Cell {cell.name} added to region {host_region} with bounds {bounds}")


        else:
            print(f"$$- CAR2D: Error : host region {host_region} not found in the geometry")
        return
    
    def getPosFromRegion(self, region):
        """
        Get the position of the region in the geometry
        :param region: region number
        :return: position of the region in the geometry
        """
        for i in range(len(self.region_numbering)):
            if self.region_numbering[i][3] == region:
                return self.region_numbering[i][0], self.region_numbering[i][1], self.region_numbering[i][2]
        return
    def getBoundsFromPos(self, i, j, k):
        """
        Get the bounds of the region from its position
        :param i: x position
        :param j: y position
        :param k: z position
        :return: bounds of the region
        """
        return self.meshx[i-1], self.meshx[i], self.meshy[j-1], self.meshy[j], self.meshz[k-1], self.meshz[k]
    
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
            bounds = self.getBoundsFromPos(i, j, k)
            print(f"$$- CAR2D: Region {region} bounds: {bounds}")
            
            temp = [] # temporary list to store the surfaces bounding the region
            # Get the surfaces bounding the region
            for surface in self.all_surfaces:
                print
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


    def describeGeo(self):
        print(f"$$- CAR2D: Created geometry (GEO: CAR2D) object of name {self.name} and level {self.level}, with {self.lx} cells in x direction and {self.ly} cells in y direction")
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


class CARCEL(GEO):
    def __init__(self, name, level, lr, radii, meshx, meshy):
        """
        Definition of a 2D cartesian geometry according to the Dragon5 definition
        :param name: name of the geometry
        :param lr: number of cells in r direction
        :param radii: radii of the cells
        :param meshx: mesh in x direction : array of nx+1 values, bounds to the cartesian cell
        :param meshy: mesh in y direction : array of ny+1 values, bounds to the cartesian cell
        """
        super().__init__(name, level)
        self.type = "CARCEL"
        self.lr = lr
        self.radii = radii
        self.meshx = meshx
        self.meshy = meshy

        self.pitch_x = self.meshx[-1] - self.meshx[0]
        self.pitch_y = self.meshy[-1] - self.meshy[0]

        self.xmin = self.meshx[0]
        self.xmax = self.meshx[-1]
        self.ymin = self.meshy[0]
        self.ymax = self.meshy[-1]

        # Trying to catch exceptions when creating the CARCEL geometry object

        if len(self.meshx) != 2 or len(self.meshy) != 2:
            print("$$-CARCEL : Error in mesh definition, need to include 2 values for the mesh")
            print("$$-CARCEL : SPLITX and SPLITY not compatible with this prototype (aiming to build a simple self-shielding geometry)")

        if self.radii[0] != 0.0:
            print("$$-CARCEL : Error in radii definition, need to include 0.0 as the first value")
        if 2*self.radii[-1] > self.pitch_x or 2*self.radii[-1] > self.pitch_y:
            print("$$-CARCEL : Error in radii definition, outermost radius too large for the bounding cartesian box")

        if self.lr != len(self.radii)-1:
            print("$$-CARCEL : Error in radii definition")


        self.number_regions()
        self.number_outer_surfaces()

        self.createInnerSurfaces()

    def setHostRegion(self, host_region):
        if self.level == 1:
            print("$$-CARCEL : Error : host region not defined for level 1")
            return
        else:
            self.host_region = host_region
        return


    def number_regions(self):
        """
        In a CARCEL geometry, the regions are numbered from the center to the outermost region, from 1 to lr+1, SPLITX and SPLITY not compatible (aiming to build a simple self-shielding geometry)
        """
        self.regions = []
        for i in range(self.lr+1):
            self.regions.append(i+1)
        self.regions = np.array(self.regions)
        return
    
    def number_outer_surfaces(self):
        """
        In a CARCEL geometry without SPLITX and SPLITY, the outer surfaces are numbered from -1 to -4
        First surfaces is the xmin bound, second is the xmax bound, third is the ymin bound, fourth is the ymax bound
        """
        self.local_bounding_surfaces = np.array([self.xmin, self.xmax, self.ymin, self.ymax])
        # Surface connectivity array, first column is index of the surface in the surfaces array first column is the surface number, 
        self.local_surface_numbering = np.array([-1, -2, -3, -4])
        self.number_of_bounding_surfaces = len(self.local_bounding_surfaces)
        
        return

    def createInnerSurfaces(self):
        """
        The CARCEL object allows for a radial meshing, with the inner surfaces being defined as the interfaces between the regions
        These are the n-1 first regions in the numbering given by self.number_regions()
        The radii provided as inout can be used to create the inner cylindrical surfaces.
        Inner surfaces are numbered from 1 to lr-1, they are stored in a list of n-uples with :
            - The first element being the region number, 
            - The second being the radius, 
            - The third being the center's x coordinate,
            - The fourth being the center's y coordinate.

        This definition assumes that the geometry is not defined using the OFFCENTER or CLUSTER Dragon5 keywords.
        Will need to adapt this for defining the B4C cylinders of the AT10 control cross. 
        """
        self.inner_surfaces = []
        for i in range(self.lr):
            self.inner_surfaces.append([i+1, self.radii[i+1], (self.xmax-self.xmin)/2, (self.ymax-self.ymin)/2])
        self.inner_surfaces = np.array(self.inner_surfaces)
        return

    

    def describeGeo(self):
        print(f"$$-CARCEL : Created geometry (GEO: CARCEL) object of name {self.name} and level {self.level}, with {self.lr} mesh points in r direction and radii {self.radii}")
        print(f"$$-CARCEL : Regions numbering from center to outermost region: {self.regions}, {len(self.regions)} regions in total")
        print(f"$$-CARCEL : Outer surfaces (x/y): {self.local_bounding_surfaces}")
        print(f"$$-CARCEL : Bounding surface labelling : {self.local_surface_numbering}")
        for i in range(len(self.local_bounding_surfaces)):
            if i <= self.number_of_bounding_surfaces/2:
                print(f"$$-CARCEL : Surface x={self.local_bounding_surfaces[i]}, with number {self.local_surface_numbering[i]}")
            else:
                print(f"$$-CARCEL : Surface y={self.local_bounding_surfaces[i]}, with number {self.local_surface_numbering[i]}")
        if self.createInnerSurfaces:
            print(f"$$-CARCEL : Inner surfaces: {self.inner_surfaces}")

    def plotGeo(self):
        """
        Plot the geometry
        """
        fig, ax = plt.subplots()
        ax.set_xlim(self.xmin-0.5, self.xmax+0.5)
        ax.set_ylim(self.ymin-0.5, self.ymax+0.5)
        ax.set_aspect('equal')
        ax.set_title(f"Geometry {self.name}")
        # Plot the bounding surfaces
        print(f"$$- GEO: Plotting bounding surfaces: {self.bounding_surfaces}")

        for i in range(len(self.bounding_surfaces)):
            if i >=2:
                ax.plot([self.bounding_surfaces[0], self.bounding_surfaces[1]], [self.bounding_surfaces[i], self.bounding_surfaces[i]], 'r-')
            else:
                ax.plot([self.bounding_surfaces[i], self.bounding_surfaces[i]], [self.bounding_surfaces[2], self.bounding_surfaces[3]], 'r-')
        
        # Plot the inner surfaces
        print(f"$$- GEO: Plotting inner surfaces: {self.inner_surfaces}")
        for i in range(len(self.inner_surfaces)):
            circle = plt.Circle((self.inner_surfaces[i][2], self.inner_surfaces[i][3]), self.inner_surfaces[i][1], fill=False)
            ax.add_artist(circle)
        plt.savefig(f"GEO_{self.name}.png")
        plt.show()
        
        return
