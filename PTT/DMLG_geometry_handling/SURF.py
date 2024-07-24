# Python3 class to define a surface object
# Author: R. Guasch
# Date : 23 July 2024
# In 2D applications, a surface is defined by a line segment.
# 

class surface:
    def __init__(self, label, region, level, points):
        self.label = label
        self.host_region = region
        self.geom_level = level
        self.points = points # 2 points are sufficient to define the surface

    def describeSurface(self):
        print("Surface id : ", self.label)
        print("Level : ", self.geom_level)
        print("Host region : ", self.host_region)
        print("Points defining the surface : ", self.points)

