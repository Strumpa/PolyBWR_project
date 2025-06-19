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

    def setHostGeometry(self, host_geometry):
        self.host_geometry = host_geometry
        return
    
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
                print(f"$$$- GEO: Sub-Geometry {self.sub_geometries[i].name} with level {self.sub_geometries[i].level}, of type {self.sub_geometries[i].type}")
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
                            print(f"Geo of type : {self.sub_geometries[i].sub_geometries[j].type}")
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
                            print(f"k= {k}")
                            print(f"$$$- GEO: Radial surface {self.sub_geometries[i].sub_geometries[j].radial_surfaces[k]}")
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
                            elif self.sub_geometries[i].sub_geometries[j].name == "C8":
                                c = "blue"
                            elif "C_UOX" in self.sub_geometries[i].sub_geometries[j].name:
                                c = "blue"
                            elif "C_Gd"in self.sub_geometries[i].sub_geometries[j].name:
                                c = "red"
                            print(self.sub_geometries[i].sub_geometries[j].name)
                            circle = plt.Circle((self.sub_geometries[i].sub_geometries[j].radial_surfaces[k][2]+min(x_displacements), self.sub_geometries[i].sub_geometries[j].radial_surfaces[k][3]+min(y_displacements)), self.sub_geometries[i].sub_geometries[j].radial_surfaces[k][1], fill=False, color=c)
                            ax.add_artist(circle)
                
            plt.savefig(f"GEO_{self.name}_lvl3.png")
            #plt.show()

        return
    

    def plotter(self, color_map = None):
        """
        Previous iteration of plotting function needs to be generalized
        Idea : test for what each sub-geometry is, and plot it accordingly
        for now can be either a CAR2D or a CARCEL

        """

        print(f"$$- GEO: Plotting geometry {self.name}")
        
        fig, axs = plt.subplots(1, 3, figsize=(15, 4))
        # ax2[0] = 1st level geometry
        # ax2[1] = 2nd level geometry
        # ax2[2] = 3rd level geometry
        axs[0].set_xlim(self.xmin-0.2, self.xmax+0.2)
        axs[1].set_xlim(self.xmin-0.2, self.xmax+0.2)
        axs[2].set_xlim(self.xmin-0.2, self.xmax+0.2)
        axs[0].set_ylim(self.ymin-0.2, self.ymax+0.2)
        axs[1].set_ylim(self.ymin-0.2, self.ymax+0.2)
        axs[2].set_ylim(self.ymin-0.2, self.ymax+0.2)
        axs[0].set_aspect('equal')
        axs[1].set_aspect('equal')
        axs[2].set_aspect('equal')
        axs[0].set_title(f"Geometry {self.name} level 1")
        axs[1].set_title(f"Geometry {self.name} level 2")
        axs[2].set_title(f"Geometry {self.name} level 3")
        fig_all, ax_all = plt.subplots(figsize=(10, 6)) # all levels in the same plot
        ax_all.set_xlim(self.xmin-0.2, self.xmax+0.2)
        ax_all.set_ylim(self.ymin-0.2, self.ymax+0.2)
        ax_all.set_aspect('equal')
        ax_all.set_title(f"Geometry {self.name} all levels")
        # Plot the bounding surfaces
        print(f"$$- GEO: Plotting bounding surfaces: {self.bounding_surfaces}")

        for i in range(len(self.bounding_surfaces)):
            if self.bounding_surfaces[i][1]==self.ymin or self.bounding_surfaces[i][1]==self.ymax:
                axs[0].plot([self.xmin, self.xmax], [self.bounding_surfaces[i][1], self.bounding_surfaces[i][1]], 'r-')
                ax_all.plot([self.xmin, self.xmax], [self.bounding_surfaces[i][1], self.bounding_surfaces[i][1]], 'r-')
            else:
                axs[0].plot([self.bounding_surfaces[i][0], self.bounding_surfaces[i][0]], [self.ymin, self.ymax], 'r-')
                ax_all.plot([self.bounding_surfaces[i][0], self.bounding_surfaces[i][0]], [self.ymin, self.ymax], 'r-')
        for i in range(len(self.all_surfaces)):
            print(f"$$$- GEO: Surface {self.all_surfaces[i]}")
            is_x_oriented, a, b, c = self.getSurfaceGeometricalData(self.all_surfaces[i][3])
            if is_x_oriented:
                x, y1, y2 = a, b, c
                axs[0].plot([x, x], [y1, y2], 'b-')
                ax_all.plot([x, x], [y1, y2], 'b-')
            else:
                y, x1, x2 = a, b, c
                axs[0].plot([x1, x2], [y, y], 'r-')
                ax_all.plot([x1, x2], [y, y], 'r-')

        # Plot the inner surfaces in a new plot
        if self.sub_geometries:
            print(f"$$- GEO: for self.name = {self.name},  number of sub-geometries is: {len(self.sub_geometries)}")
            # Plot the bounding surfaces
            for geom in self.sub_geometries: # level 2 geometries analysis
                print(f"$$- GEO: Plotting sub-geometry {geom.name} with level {geom.level}")
                region_id = geom.host_region
                geom_coords = self.getBoundsFromPosition(i=self.getPositionFromRegion(region_id)[0], 
                                                         j=self.getPositionFromRegion(region_id)[1],
                                                         k=self.getPositionFromRegion(region_id)[2])
                
                x_offset = geom_coords[0]
                y_offset = geom_coords[2]
                if geom.type == "CAR2D":
                    # Recover data for the sub-geometry's surfaces in the CAR2D case
                    for surf_id in range(len(geom.all_surfaces)):
                        print(f"$$$- GEO: Surface {geom.all_surfaces[surf_id]}")
                        is_x_oriented, a, b, c = geom.getSurfaceGeometricalData(geom.all_surfaces[surf_id][3])
                        if is_x_oriented:
                            x, y1, y2 = a + x_offset, b + y_offset, c + y_offset
                            axs[1].plot([x, x], [y1, y2], 'b-')
                            ax_all.plot([x, x], [y1, y2], 'b-')
                        else:
                            y, x1, x2 = a + y_offset, b + x_offset, c + x_offset
                            axs[1].plot([x1, x2], [y, y], 'r-')
                            ax_all.plot([x1, x2], [y, y], 'r-')
                elif geom.type == "CARCEL":
                    x_displacements = []
                    y_displacements = []
                    for surf_index in range(len(geom.host_bounding_surfaces)):
                        print("surf_index is ", surf_index)
                        print(f"surf_index in host surfaces = {geom.host_bounding_surfaces[surf_index]}")
                        print(f"surf_index in local surfaces = {geom.local_surface_numbering[surf_index]}")
                        # Recover data for the host bounding surfaces :
                        is_x_oriented, a, b, c = self.getSurfaceGeometricalData(geom.host_bounding_surfaces[surf_index])
                        if is_x_oriented:
                            x, y1, y2 = a + x_offset, b + y_offset, c + y_offset
                            x_displacements.append(x)
                            axs[1].plot([x, x], [y1, y2], 'b-')
                            ax_all.plot([x, x], [y1, y2], 'b-')
                        else:
                            y, x1, x2 = a + y_offset, b + x_offset, c + x_offset
                            y_displacements.append(y)
                            axs[1].plot([x1, x2], [y, y], 'r-')
                            ax_all.plot([x1, x2], [y, y], 'r-')

                    for radial_surf_index in range(len(geom.radial_surfaces)):
                            print(f"radial_surf_index= {radial_surf_index}")
                            print(f"$$$- GEO: Radial surface {geom.radial_surfaces[radial_surf_index]}")
                            if color_map:
                                print(f"geom.name = {geom.name}, geom.name[:-1] = {geom.name[:-1]}")
                                c = color_map[geom.name[:-1]]
                            else:
                                if "C_UOX" in geom.name:
                                    c = "blue"
                                elif "C_Gd"in geom.name:
                                    c = "red"
                                else: 
                                    c = "black"
                            circle = plt.Circle((geom.radial_surfaces[radial_surf_index][2]+min(x_displacements), geom.radial_surfaces[radial_surf_index][3]+min(y_displacements)), geom.radial_surfaces[radial_surf_index][1], fill=False, color=c)
                            circle2 = plt.Circle((geom.radial_surfaces[radial_surf_index][2]+min(x_displacements), geom.radial_surfaces[radial_surf_index][3]+min(y_displacements)), geom.radial_surfaces[radial_surf_index][1], fill=False, color=c)
                            axs[1].add_artist(circle)
                            ax_all.add_artist(circle2)
                else:
                    print(f"$$- GEO: Unknown geometry type {geom.type}")
            
                if geom.sub_geometries: # analysis of level 3 geometries
                    print(f"$$- GEO: for sub-geometry = {geom.name},  number of sub-sub-geometries is: {len(geom.sub_geometries)}")
                    print(f"$$- GEO: Plotting sub-geometries of {geom.name}")
                    # Plot the bounding surfaces
                    for geom_n3 in geom.sub_geometries:
                        print(f"$$- GEO: Plotting sub-geometry {geom_n3.name} with level {geom_n3.level}")
                        print(f"$$- GEO: Plotting sub level geometry {geom_n3.name} with level {geom_n3.level} and regionid {geom_n3.host_region}")
                        region_id = geom_n3.host_region
                        geom_coords = geom.getBoundsFromPosition(i=geom.getPositionFromRegion(region_id)[0], 
                                                            j=geom.getPositionFromRegion(region_id)[1],
                                                            k=geom.getPositionFromRegion(region_id)[2])
                        local_x_offset = geom_coords[0]
                        local_y_offset = geom_coords[2]
                        print(f"$$- GEO: in macro geometry {geom.name}, sub-geometry {geom_n3.name} has host region {region_id}")
                        print(f"$$- GEO: in macro geometry {geom.name}, x_offset = {x_offset}, y_offset = {y_offset}")
                        print(f"$$- GEO: local_x_offset = {local_x_offset}, local_y_offset = {local_y_offset}")
                        if geom_n3.type == "CAR2D":
                            # Recover data for the sub-geometry's surfaces in the CAR2D case
                            for surf_id in range(len(geom_n3.all_surfaces)):
                                print(f"$$$- GEO: Surface {geom_n3.all_surfaces[surf_id]}")
                                is_x_oriented, a, b, c = geom_n3.getSurfaceGeometricalData(geom_n3.all_surfaces[surf_id][3])
                                if is_x_oriented:
                                    x, y1, y2 = a + x_offset + local_x_offset, b + y_offset + local_y_offset, c + y_offset + local_y_offset
                                    axs[2].plot([x, x], [y1, y2], 'b-')
                                    ax_all.plot([x, x], [y1, y2], 'b-')
                                else:
                                    y, x1, x2 = a + y_offset + local_y_offset, b + x_offset + local_x_offset, c + x_offset + local_x_offset
                                    axs[2].plot([x1, x2], [y, y], 'r-')
                                    ax_all.plot([x1, x2], [y, y], 'r-')
                        elif geom_n3.type == "CARCEL":
                            print(f"$$- GEO: Plotting sub-geometry {geom_n3.name} with level {geom_n3.level} and regionid {geom_n3.host_region}, hosted in geometry {geom.name}")
                            print(f"geom_n3.host_bounding_surfaces = {geom_n3.host_bounding_surfaces}")
                            for surf_index in range(len(geom_n3.host_bounding_surfaces)):
                                # Prevent double counting of the same surface :
                                print("in plotter surf_index is ", surf_index)
                                # Recover data for the host bounding surfaces :
                                is_x_oriented, a, b, c = geom.getSurfaceGeometricalData(geom_n3.host_bounding_surfaces[surf_index])
                                print(f"$$$- GEO: Surface {geom_n3.host_bounding_surfaces[surf_index]}")
                                # Accounting for translation in 2 level geom :
                                if is_x_oriented:
                                    x, y1, y2 = a + x_offset , b + y_offset , c + y_offset
                                    if max(y1,y2) <= geom.ymax+y_offset and min(y1,y2) >= geom.ymin+y_offset:
                                        axs[2].plot([x, x], [y1, y2], 'b-')
                                        ax_all.plot([x, x], [y1, y2], 'b-')
                                else:
                                    y, x1, x2 = a + y_offset , b + x_offset , c + x_offset
                                    if max(x1,x2) <= geom.xmax+x_offset and min(x1,x2) >= geom.xmin+x_offset:
                                        axs[2].plot([x1, x2], [y, y], 'r-')
                                        ax_all.plot([x1, x2], [y, y], 'r-')

                            for radial_surf_index in range(len(geom_n3.radial_surfaces)):
                                if color_map:
                                    print(f"geom_n3.name = {geom_n3.name}, geom_n3.name[:-1] = {geom_n3.name[:-1]}")
                                    c = color_map[geom_n3.name[:-1]]
                                else:
                                    if "C_UOX" in geom_n3.name:
                                        c = "blue"
                                    elif "C_Gd"in geom_n3.name:
                                        c = "red"
                                    else: 
                                        c = "black"
                                pin_x_offset = x_offset + local_x_offset
                                pin_y_offset = y_offset + local_y_offset
                                circle = plt.Circle((geom_n3.radial_surfaces[radial_surf_index][2]+pin_x_offset, geom_n3.radial_surfaces[radial_surf_index][3]+pin_y_offset), geom_n3.radial_surfaces[radial_surf_index][1], fill=False, color=c)
                                circle2 = plt.Circle((geom_n3.radial_surfaces[radial_surf_index][2]+pin_x_offset, geom_n3.radial_surfaces[radial_surf_index][3]+pin_y_offset), geom_n3.radial_surfaces[radial_surf_index][1], fill=False, color=c)
                                axs[2].add_artist(circle)
                                ax_all.add_artist(circle2)
                        else:
                            print(f"$$- GEO: Unknown geometry type {geom_n3.type}")

            #plt.savefig(f"GEO_{self.name}_lvl2.png")
        fig.savefig(f"GEO_{self.name}_multi_level.png")
        fig_all.savefig(f"GEO_{self.name}_all_levels.png")


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

