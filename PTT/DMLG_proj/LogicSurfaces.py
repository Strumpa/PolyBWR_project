# Abstract relation class for surfaces.
# Author : R. Guasch
# Aim is to be classify surfaces according to their logical relations.


class LogicSurfaces:
    def __init__(self):
        self.surfaces = []
    
    def add_surface(self, surface):
        self.surfaces.append(surface)
    
    def get_surfaces(self):
        return self.surfaces
    
    def get_surface_idx(self, index):
        return self.surfaces[index]
    
    def set_relation(self, relation, list_surfaces):
        if relation == "Complement":
            print("Complement")
        elif relation == "Intersection":
            print("Intersection")
        elif relation == "Union":
            print("Union")
        elif relation == "ConnectedEdges":
            print("Connectivity")
    
    def get_relation(self):
        return self.relation
    
    

