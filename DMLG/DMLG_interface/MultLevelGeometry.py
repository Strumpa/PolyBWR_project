# Recursive definition of Multi-level geometry
# Author: R. Guasch
# Aim is to be able to define a multi-level geometry and create connectivity and numbering tables from it.

class MultiLevelGeom:
    def __init__(self, nb_levels, geom_type):
        self.nb_levels = nb_levels
        self.type = geom_type
        self.meshx = []
        self.meshy = []

        self.connectivity = []
        self.numbering = []

        self.cells = []

    def add_cell(self, cell):
        """
        Adding cell to geometry, cell is a geometry object, can itself be an assembly of cells.
        """

        self.cells.append(cell)


class Cell:
    def __init__(self, cell_type, cell_id, nb_levels):
        self.cell_type = cell_type
        self.cell_id = cell_id
        self.nb_levels = nb_levels
        self.meshx = []
        self.meshy = []

        self.connectivity = []
        self.numbering = []

        self.cells = []

    def add_cell(self, cell):
        """
        Adding cell to geometry, cell is a geometry object
        """

        self.cells.append(cell)



