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
    
    def setHostBoundingSurfaces(self, surfaces):
        self.host_bounding_surfaces = surfaces
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
                    tmp.append([self.level, region, cell.host_bounding_surfaces, cell.local_surface_numbering])
                    if cell.sub_geometries:
                        for sub_cell in cell.sub_geometries:
                            tmp.append([cell.level, cell.host_region, sub_cell.host_bounding_surfaces, sub_cell.local_surface_numbering])

            self.connectivity_map.append(tmp)

        print(f"$$- GEO: Connectivity map between regions and cells: {self.connectivity_map}")
        return
    
    def analyzeLevels(self):
        """
        Analyze the levels of the geometry and its sub-geometries
        """
        print(f"$$- GEO: Analyzing levels of the geometry {self.name}")
        print(f"$$- GEO: Level of the geometry: {self.level}")
        self.max_level = [1]
        if self.sub_geometries:
            for cell in self.sub_geometries:
                print(f"$$- GEO: Level of cell {cell.name}: {cell.level}")
                if cell.sub_geometries:
                    for sub_cell in cell.sub_geometries:
                        print(f"$$- GEO: Level of sub-cell {sub_cell.name}: {sub_cell.level}")
                        self.max_level.append(max(self.level, cell.level, sub_cell.level))

        self.max_level = max(self.max_level)                        
        return
    
    def plotGeo(self):
        """
        Recover information about the geometry then plot the geometry
        """


        fig, ax = plt.subplots()
        ax.set_xlim(self.xmin-0.2, self.xmax+0.2)
        ax.set_ylim(self.ymin-0.2, self.ymax+0.2)
        ax.set_aspect('equal')
        ax.set_title(f"Geometry {self.name}")
        # Plot the bounding surfaces
        print(f"$$- GEO: Plotting bounding surfaces: {self.bounding_surfaces}")

        for i in range(len(self.bounding_surfaces)):
            if self.bounding_surfaces[i][1]==self.ymin or self.bounding_surfaces[i][1]==self.ymax:
                ax.plot([self.xmin, self.xmax], [self.bounding_surfaces[i][1], self.bounding_surfaces[i][1]], 'r-')
            else:
                ax.plot([self.bounding_surfaces[i][0], self.bounding_surfaces[i][0]], [self.ymin, self.ymax], 'r-')
        plt.savefig(f"GEO_{self.name}_lvl1.png")

        # Plot the inner surfaces
        if self.sub_geometries:
            for i in range(len(self.sub_geometries)):
                print(f"GEO_plot: i = {i}")
                #self.sub_geometries[i].plotGeo() 
                print(f"$$$- GEO: Sub-Geometry {self.sub_geometries[i].name} with level {self.sub_geometries[i].level}")
                print(f"$$- GEO: Plotting inner bounding surfaces: {self.sub_geometries[i].bounding_surfaces}")
                print(f"$$- GEO: Plotting inner surfaces: {self.sub_geometries[i].all_surfaces}")
                for surf_id in range(len(self.sub_geometries[i].all_surfaces)):
                    print(f"$$$- GEO: Surface {self.sub_geometries[i].all_surfaces[surf_id]}")
                    is_x_oriented, a, b, c = self.sub_geometries[i].getSurfaceGeometricalData(self.sub_geometries[i].all_surfaces[surf_id][3])
                    if is_x_oriented:
                        x, y1, y2 = a, b, c
                        ax.plot([x, x], [y1, y2], 'b-')
                    else:
                        y, x1, x2 = a, b, c
                        ax.plot([x1, x2], [y, y], 'r-')
                    plt.savefig(f"GEO_{self.name}_lvl2.png")
                
                if self.sub_geometries[i].sub_geometries:
                    for j in range(len(self.sub_geometries[i].sub_geometries)):
                        print(f"$$$- GEO: Sub-Geometry {self.sub_geometries[i].sub_geometries[j].name} with level {self.sub_geometries[i].sub_geometries[j].level}")
                        print(f"$$- GEO: Plotting inner bounding surfaces: {self.sub_geometries[i].sub_geometries[j].host_bounding_surfaces}")
                        # Need to find proper translation to bring the circle's centers at the right position
                        x_displacements = []
                        y_displacements = []
                        for surf_index in range(len(self.sub_geometries[i].sub_geometries[j].host_bounding_surfaces)):
                            print("surf_index is ", surf_index)
                            print(f"surf_index in host surfaces = {self.sub_geometries[i].sub_geometries[j].host_bounding_surfaces[surf_index]}")
                            print(f"At n-1 level : host region = {self.sub_geometries[i].sub_geometries[j].host_region}, with local surfaces {self.sub_geometries[i].all_surfaces}")
                            print(f"surf_index in local surfaces = {self.sub_geometries[i].sub_geometries[j].local_surface_numbering[surf_index]}")
                            # Recover data for the host bounding surfaces :
                            is_x_oriented, a, b, c = self.sub_geometries[i].getSurfaceGeometricalData(self.sub_geometries[i].sub_geometries[j].host_bounding_surfaces[surf_index])
                            if is_x_oriented:
                                x, y1, y2 = a, b, c
                                x_displacements.append(x)
                            else:
                                y, x1, x2 = a, b, c
                                y_displacements.append(y)
                        print(f"$$- GEO: for geometry {self.sub_geometries[i].sub_geometries[j].name}, x_displacements = {x_displacements}, y_displacement = {y_displacements}")
                        print(f"$$- GEO: Plotting radial surfaces: {self.sub_geometries[i].sub_geometries[j].radial_surfaces}")
                        for k in range(len(self.sub_geometries[i].sub_geometries[j].radial_surfaces)):
                            if self.sub_geometries[i].sub_geometries[j].name == "C1":
                                c = "darkorange"
                            elif self.sub_geometries[i].sub_geometries[j].name == "C2":
                                c = "darkgreen"
                            elif self.sub_geometries[i].sub_geometries[j].name == "C3":
                                c = "indigo"
                            elif self.sub_geometries[i].sub_geometries[j].name == "C4":
                                c = "red"
                            elif self.sub_geometries[i].sub_geometries[j].name == "C6":
                                c = "violet"
                            elif self.sub_geometries[i].sub_geometries[j].name == "C7":
                                c = "gold"
                            
                            circle = plt.Circle((self.sub_geometries[i].sub_geometries[j].radial_surfaces[k][2]+min(x_displacements), self.sub_geometries[i].sub_geometries[j].radial_surfaces[k][3]+min(y_displacements)), self.sub_geometries[i].sub_geometries[j].radial_surfaces[k][1], fill=False, color=c)
                            ax.add_artist(circle)
                
            plt.savefig(f"GEO_{self.name}_lvl3.png")
            plt.show()

        return
    
    def getPositionFromRegion(self, region):
        """
        Get the position of the region in the geometry
        :param region: region number
        :return: position of the region in the geometry
        """
        for i in range(len(self.region_numbering)):
            if self.region_numbering[i][3] == region:
                return self.region_numbering[i][0], self.region_numbering[i][1], self.region_numbering[i][2]
        return
    
    def getBoundsFromPosition(self, i, j, k):
        """
        Get the bounds of the region from its position
        :param i: x position
        :param j: y position
        :param k: z position
        :return: bounds of the region
        """
        return self.meshx[i-1], self.meshx[i], self.meshy[j-1], self.meshy[j], self.meshz[k-1], self.meshz[k]

    def getRegionCenter(self, region):
        """
        Recover the region's geometrical center from its number/label,
        getPositionFromRegion() returns the position of the region in the geometry,
        getBoundsFromPosition() returns the bounds of the region from its position
        :param region: region number

        The geometrical center can then be defined as the middle point between the bounds of the region
        """
        i, j, k = self.getPositionFromRegion(region)
        bounds = self.getBoundsFromPosition(i, j, k)
        x_center = (bounds[0]+bounds[1])/2
        y_center = (bounds[2]+bounds[3])/2
        #z_center = (bounds[4]+bounds[5])/2
        return x_center, y_center
    
    def translateToCoordinateSystem(self, x, y, region):
        """
        Translate the coordinates x and y to the coordinate system of the region
        :param x: x coordinate
        :param y: y coordinate
        :param region: region number
        :return: translated coordinates
        """
        i, j, k = self.getPositionFromRegion(region)
        return x - self.meshx[i-1], y - self.meshy[j-1]

