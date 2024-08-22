# Python3 class to represent a DRAGON5 CARCEL geometry object
# Author: R. Guasch
# Date : 24 July 2024
# Dependencies : GEOM.py


from GEOM import GEO
from matplotlib import pyplot as plt
import numpy as np

class CARCEL(GEO):
    def __init__(self, name, level, nr, radii, meshx, meshy):
        """
        Definition of a 2D cartesian geometry according to the Dragon5 definition
        :param name: name of the geometry
        :param nr: number of cells in r direction
        :param radii: radii of the cells
        :param meshx: mesh in x direction : array of nx+1 values, bounds to the cartesian cell
        :param meshy: mesh in y direction : array of ny+1 values, bounds to the cartesian cell
        """
        super().__init__(name, level)
        self.type = "CARCEL"
        self.nr = nr
        self.radii = radii
        self.meshx = meshx
        self.meshy = meshy

        self.pitch_x = self.meshx[-1] - self.meshx[0]
        self.pitch_y = self.meshy[-1] - self.meshy[0]

        self.xmin = self.meshx[0]
        self.xmax = self.meshx[-1]
        self.ymin = self.meshy[0]
        self.ymax = self.meshy[-1]

        self.origin = [self.xmin, self.ymin]
        self.center = [(self.xmax-self.xmin)/2, (self.ymax-self.ymin)/2]

        # Trying to catch exceptions when creating the CARCEL geometry object

        if len(self.meshx) != 2 or len(self.meshy) != 2:
            print("$$-CARCEL : Error in mesh definition, need to include 2 values for the mesh")
            print("$$-CARCEL : SPLITX and SPLITY not compatible with this prototype (aiming to build a simple self-shielding geometry)")

        if self.radii[0] != 0.0:
            print("$$-CARCEL : Error in radii definition, need to include 0.0 as the first value")
        if 2*self.radii[-1] > self.pitch_x or 2*self.radii[-1] > self.pitch_y:
            print("$$-CARCEL : Error in radii definition, outermost radius too large for the bounding cartesian box")

        if self.nr != len(self.radii)-1:
            print("$$-CARCEL : Error in radii definition")


        self.number_regions()
        self.number_outer_surfaces()

        self.createRadialSurfaces()
        self.bounding_surfaces = self.local_bounding_surfaces

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
        for i in range(self.nr+1):
            self.regions.append(i+1)
        self.regions = np.array(self.regions)
        self.nb_regions = len(self.regions)
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

    def createRadialSurfaces(self):
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
        self.radial_surfaces = []
        for i in range(self.nr):
            self.radial_surfaces.append([i+1, self.radii[i+1], (self.xmax-self.xmin)/2, (self.ymax-self.ymin)/2])
        self.radial_surfaces = np.array(self.radial_surfaces)
        return

    

    def describeGeo(self):
        print(f"$$-CARCEL : Created geometry (GEO: CARCEL) object of name {self.name} and level {self.level}, with {self.nr} mesh points in r direction and radii {self.radii}")
        print(f"$$-CARCEL : Regions numbering from center to outermost region: {self.regions}, {len(self.regions)} regions in total")
        print(f"$$-CARCEL : Outer surfaces (x/y): {self.local_bounding_surfaces}")
        print(f"$$-CARCEL : Bounding surface labelling : {self.local_surface_numbering}")
        for i in range(len(self.local_bounding_surfaces)):
            if i <= self.number_of_bounding_surfaces/2:
                print(f"$$-CARCEL : Surface x={self.local_bounding_surfaces[i]}, with number {self.local_surface_numbering[i]}")
            else:
                print(f"$$-CARCEL : Surface y={self.local_bounding_surfaces[i]}, with number {self.local_surface_numbering[i]}")
        if self.radial_surfaces:
            print(f"$$-CARCEL : Inner surfaces: {self.radial_surfaces}")

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
        print(f"$$- GEO: Plotting radial surfaces: {self.radial_surfaces}")
        for i in range(len(self.radial_surfaces)):
            circle = plt.Circle((self.radial_surfaces[i][2], self.radial_surfaces[i][3]), self.radial_surfaces[i][1], fill=False)
            ax.add_artist(circle)
        plt.savefig(f"GEO_{self.name}.png")
        plt.show()
        
        return
