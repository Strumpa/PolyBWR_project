# Python3 script to prototype the NXT_CLUSTER class,
# Purpose : define a class capable of handling an input multi-level geometry,
# If needed : define the third geometry level as a cluster of pins.
# Author: R. Guasch
# Date : 24 July 2024


class NXTGAC:
    """
    NXTGAC class : NXT Geometry Analysis and Clustering
    """
    def __init__(self, geommetry):
        """
        Initialize the NXT_CLUSTER class, only input needed for our application is the geometry object,
        In real NXT applications, the tracking parameters would also be needed.
        """
        print("$$$- NXTGAC Initializing the NXT_CLUSTER class")
        self.geometry = geommetry
        geommetry.analyzeLevels()
        self.max_level = self.geometry.max_level
        self.analyse_geometry()
        

    def analyse_geometry(self):
        """
        Analyse the geometry and create the necessary data structures
        """
        print("$$- NXTGAC Analyzing the geometry")
        print(f"$$- NXTGAC Geometry has {self.max_level} levels, of type {self.geometry.type}")
        print(f"$$- NXTGAC Geometry is bounded by surfaces : {self.geometry.bounding_surfaces}")
        print(f"$$- NXTGAC Geometry has {self.geometry.nb_regions} regions : {self.geometry.regions}")
        print(f"$$- NXTGAC Geometry has {self.geometry.number_of_bounding_surfaces} outer surfaces")

        if self.max_level == 1:
            print("$$- NXT : Geometry is a simple 1 geometry")
            print("This geometry doesn't need further analysis or modifications : ")
            print("Can proceed with analysis and tracking using the regular NXT: procedure.")

        elif self.max_level == 2:
            print("$$- NXT : Geometry is a 2 level geometry")
            print("This geometry doesn't need further analysis or modifications : ")
            print("Can proceed with analysis and tracking using the regular NXT: procedure.")
        
        elif self.max_level == 3:
            print("$$- NXT : Geometry is a 3 level geometry")
            print("This geometry requires further analysis, the 3rd level is to be defined as a cluster of pins.")
            self.analyse_third_level()

        return

    def analyse_third_level(self):
        """
        
        """


        return



        