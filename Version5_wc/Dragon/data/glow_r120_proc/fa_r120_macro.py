"""
Fuel Assembly with BP pins, representative of TL40.
"""
import math
import time
from glow.main import *
from glow.geometry_layouts.fillable_layouts import Fillable
from glow.geometry_layouts.lattices import Lattice
from glow.geometry_layouts.layouts import *
from glow.geometry_layouts.cells import HexCell
from glow.geometry_layouts.lattices import HexLattice
# from glow.main import analyse_and_generate_tdt
from glow.interface.geom_interface import *
from glow.geometry_layouts.geometries import Circle, Hexagon, build_parallelogram
from glow.support.utility import *
from glow.support.types import *


# Get the initial time
start_time = time.time()

# Declare the values of the hexagonal cells geometrical characteristics
pins_lattice_pitch = 1.32  #cm #twice the apothem
edge_length_f = pins_lattice_pitch/math.sqrt(3) #l=2a/sqrt(3)
radii_f = [0.45, 0.46, 0.53] #fuel cell
wrapper_i_2a = 19.8
wrapper_o_2a = 20.1
pitch_lattice_core = 20.5

edge_wrapper_i = wrapper_i_2a / math.sqrt(3)
edge_wrapper_o = wrapper_o_2a / math.sqrt(3)
edge_bypass = pitch_lattice_core / math.sqrt(3)

# Build hexagonal cells
fuel1_cell= HexCell(
    side=edge_length_f,
    name="Fuel1 Cell",
    base_props={PropertyType.MATERIAL: "LEAD"}
)
fuel1_cell.rotate(90)
properties_fuel_cell={
    PropertyType.MATERIAL: [
        "CLAD", "HELIUM", "FUEL1" 
    ]
}
# Add the circles representing the different zones for pellet and cladding
for i, r in enumerate(radii_f[::-1]):
    fuel1_cell.add(
        Region(
            Circle(radius=r),
            properties={
                PropertyType.MATERIAL: properties_fuel_cell[PropertyType.MATERIAL][i]}
        )
    )

lead_cell = HexCell(
    side=edge_length_f,
    name="Lead Cell",
    base_props={PropertyType.MATERIAL : "LEAD"}
)
lead_cell.rotate(90)

bp_cell = HexCell(side=edge_length_f,
                  name='BP Cell',
                  base_props={PropertyType.MATERIAL: "LEAD"})
bp_cell.rotate(90)
properties_bp_cell = {
    PropertyType.MATERIAL: [
        "CLAD", "HELIUM", "BP"
    ]
}
# Add the circles representing the different zones for pellet and cladding
for i, r in enumerate(radii_f[::-1]):
    bp_cell.add(
        Region(
            Circle(radius=r),
            properties={
                PropertyType.MATERIAL: properties_bp_cell[PropertyType.MATERIAL][i]}
        )
    )


# Update the viewer showing a color map for the MATERIAL property type
fuel1_cell.show(PropertyType.MATERIAL)
lead_cell.show(PropertyType.MATERIAL)
bp_cell.show(PropertyType.MATERIAL)

# Build the lattice
lattice_fa1 = HexLattice(cells=[lead_cell])
lattice_fa1.add_ring_of_cells(fuel1_cell, 1, 0)
lattice_fa1.add_ring_of_cells(fuel1_cell, 2, 0)
lattice_fa1.add_ring_of_cells(fuel1_cell, 3, 0)
lattice_fa1.add_ring_of_cells(fuel1_cell, 4, 0)
lattice_fa1.add_ring_of_cells(fuel1_cell, 5, 0)
lattice_fa1.add_ring_of_cells(fuel1_cell, 6, 0)
lattice_fa1.add_ring_of_cells(fuel1_cell, 7, 0)
lattice_fa1.add_ring_of_cells(fuel1_cell, 8, 0)

# Add burnable poison pins to the lattice

Rin = pins_lattice_pitch*3 / math.sqrt(3)
Rex = Rin*3

dx = lambda i: pins_lattice_pitch*4 * math.cos(math.radians(i*60))
dy = lambda i: pins_lattice_pitch*4 * math.sin(math.radians(i*60))

for i in range(3):
    x = Rin * math.cos(math.radians(30 + i*120))
    y = Rin * math.sin(math.radians(30 + i*120))
    lattice_fa1.add(bp_cell, (x, y, 0), 1)

    x = Rex * math.cos(math.radians(90 + i*120))
    y = Rex * math.sin(math.radians(90 + i*120))
    lattice_fa1.add(bp_cell, (x, y, 0), 1)

    lattice_fa1.add(bp_cell, (x+dx(i), y-dy(i), 0), 1)
    lattice_fa1.add(bp_cell, (x-dx(i), y+dy(i), 0), 1)


# Build the fuel assembly cell
fa_1 = HexCell(
    side=edge_bypass,
    name = 'Fuel Assembly 1',
    base_props={PropertyType.MATERIAL: 'LEAD'}
)
# Add the lattice first, then the cladding layer
fa_1.add(lattice_fa1)
fa_1.add(
    Region(
        Hexagon(edge_length=edge_wrapper_o)
        - Hexagon(edge_length=edge_wrapper_i),
        "Wrapper",
        {PropertyType.MATERIAL: 'CLAD'}
    )
)
fa_1.show(PropertyType.MATERIAL)

# Assign the 'MACRO' property so that the box layer and the background has
# the same value, whereas each cell of the lattice has a different value

fa_1.set_region_properties(
    {PropertyType.MACRO: "MAC_001"}, fa_1.layers[2][0] # wrapper
)

# Loop through the regions of the background of the assembly
for region in fa_1.layers[0]:
    fa_1.set_region_properties(
        {PropertyType.MACRO: "MAC_001"}, region
    )

macro_index = 1
# Loop through the layers of the lattice, and through those of each cell
for top_layer in fa_1.layers[1][0].layers: 
        for layer in top_layer: # cell
            macro_index +=1
            for cell_layer in layer.layers: 
                for cell_region in cell_layer:
                    cell_region.properties.update(
                        {PropertyType.MACRO: "MAC_00" + str(macro_index)}
                    )


# Apply a symmetry
fa_1.apply_symmetry(SymmetryType.THIRD)


# Display the assembly in the SALOME 3D viewer
print(f"--- Geometry generated in {time.time() - start_time} seconds ---")

fa_1.show(PropertyType.MACRO)

# Perform the geometry analysis and export the TDT file of the surface geometry
export_layout_to_tdt(fa_1, "fa_r120_tiso_macro", tdt_setup=TdtSetup(
    type_geo=LayoutGeometryType.ROTATION, symmetry_type=SymmetryType.THIRD, property_types= [PropertyType.MATERIAL, PropertyType.MACRO]
))

# Print the execution time
print(f"--- Code executed in {time.time() - start_time} seconds ---")
