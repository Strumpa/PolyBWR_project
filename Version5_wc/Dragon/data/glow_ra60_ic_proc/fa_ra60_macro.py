"""
Hexagonal assembly made of a central pin fuel cell sorrounded by 6 rings of
pin fuel cells. Symmetry is exploited by considering a one-sixth of the
complete hexagonal layout. Surfacic multicell approximation.
"""
import math
import time

# Import the GLOW geometry classes
from glow.geometry_layouts.layouts import *
from glow.geometry_layouts.cells import HexCell
from glow.geometry_layouts.lattices import HexLattice
from glow.main import *
from glow.interface.geom_interface import *
from glow.geometry_layouts.geometries import Circle, Hexagon
from glow.support.utility import *
from glow.support.types import *

# Get the initial time
start_time = time.time()

# Declare the values of the hexagonal cells geometrical characteristics
pins_lattice_pitch = 1.34  #cm #twice the apothem
edge_length_f = pins_lattice_pitch/math.sqrt(3) #l=2a/sqrt(3)
radii_f = [0.15, 0.35, 0.4, 0.45] #fuel cell
central_hex_o = 3.2     # twice the apothem of the outer hexagon
central_hex_i = 2.5     # twice the apothem of the inner hexagon
wrapper_i_2a = 20.1
wrapper_o_2a = 20.4
pitch_lattice_core = 21.0203115

edge_length_o = central_hex_o / math.sqrt(3)
edge_length_i = central_hex_i / math.sqrt(3)
edge_wrapper_i = wrapper_i_2a / math.sqrt(3)
edge_wrapper_o = wrapper_o_2a / math.sqrt(3)
edge_bypass = pitch_lattice_core / math.sqrt(3)

# Build hexagonal cells
fuel2_cell= HexCell(
    side=edge_length_f,
    name="Fuel1 Cell",
    base_props={PropertyType.MATERIAL: "LEAD"}
)
fuel2_cell.rotate(90)
properties={
    PropertyType.MATERIAL: [
        "CLADDING", "HELIUM", "FUEL", "HELIUM"
    ]
}
# Add the circles representing the different zones for pellet and cladding
for i, r in enumerate(radii_f[::-1]):
    fuel2_cell.add(
        Region(
            Circle(radius=r),
            properties={
                PropertyType.MATERIAL: properties[PropertyType.MATERIAL][i]}
        )
    )

lead_cell = HexCell(
    side=edge_length_f,
    name="Lead Cell",
    base_props={PropertyType.MATERIAL : "LEAD"}
)
lead_cell.rotate(90)

# Update the viewer showing a color map for the MATERIAL property type
fuel2_cell.show(PropertyType.MATERIAL)
lead_cell.show(PropertyType.MATERIAL)
fuel2_cell.sectorize([6, 6, 6, 6, 6], [0]*5)
lead_cell.sectorize([6], [0])

# Build the central channel cell
central_channel_o = HexCell(
    side=edge_length_o,
    name="Central Cell",
    base_props={PropertyType.MATERIAL: 'CLADDING'}
)
central_channel_o.add(
    Region(
        Circle(radius=edge_length_i),
        "Central Cell Inner",
        {PropertyType.MATERIAL: 'LEAD'}
    )
)
central_channel_o.show(PropertyType.MATERIAL)

# Build the lattice
lattice_fa2 = HexLattice(cells=[lead_cell])
lattice_fa2.add_ring_of_cells(lead_cell, 1, 0)
lattice_fa2.add_ring_of_cells(fuel2_cell, 2, 0)
lattice_fa2.add_ring_of_cells(fuel2_cell, 3, 0)
lattice_fa2.add_ring_of_cells(fuel2_cell, 4, 0)
lattice_fa2.add_ring_of_cells(fuel2_cell, 5, 0)
lattice_fa2.add_ring_of_cells(fuel2_cell, 6, 0)
lattice_fa2.add_ring_of_cells(fuel2_cell, 7, 0)
lattice_fa2.add_ring_of_cells(fuel2_cell, 8, 0)
# Add the cell with greater size
lattice_fa2.add(central_channel_o)

# Build the fuel assembly cell
fa_2 = HexCell(
    side=edge_bypass,
    name = 'Fuel Assembly 2',
    base_props={PropertyType.MATERIAL: 'LEAD'}
)
# Add the lattice first, then the cladding layer
fa_2.add(lattice_fa2)
fa_2.add(
    Region(
        Hexagon(edge_length=edge_wrapper_o)
        - Hexagon(edge_length=edge_wrapper_i),
        "Wrapper",
        {PropertyType.MATERIAL: 'CLADDING'}
    )
)
fa_2.show(PropertyType.MATERIAL)

# Assign the 'MACRO' property so that the box layer and the background has
# the same value, whereas each cell of the lattice has a different value

fa_2.set_region_properties(
    {PropertyType.MACRO: "MAC_001"}, fa_2.layers[2][0] # wrapper
)

# Loop through the regions of the background of the assembly
for region in fa_2.layers[0]:
    fa_2.set_region_properties(
        {PropertyType.MACRO: "MAC_001"}, region
    )

macro_index = 2
# Loop through the layers of the lattice, and through those of each cell
for idx, top_layer in enumerate(fa_2.layers[1][0].layers):

    if idx == 0:
        # first 6 (lead around central channel)
        for layer in top_layer[:6]:
            for cell_layer in layer.layers:
                for cell_region in cell_layer:
                    cell_region.properties.update(
                        {PropertyType.MACRO: "MAC_002"}
                    )

        # all other cells of the lattice
        for layer in top_layer[6:]:
            macro_index += 1
            for cell_layer in layer.layers:
                for cell_region in cell_layer:
                    cell_region.properties.update(
                        {PropertyType.MACRO: "MAC_00" + str(macro_index)}
                    )

    else:
        for layer in top_layer: # in this case, only central channel
            macro_index += 1
            for cell_layer in layer.layers:
                for cell_region in cell_layer:
                    cell_region.properties.update(
                        {PropertyType.MACRO: "MAC_002"}
                    )

fa_2.apply_symmetry(SymmetryType.SIXTH)
# Show the assembly regions according to the 'MACRO' colour map
fa_2.show(PropertyType.MACRO)

# # Display the assembly in the SALOME 3D viewer
print(f"--- Geometry generated in {time.time() - start_time} seconds ---")

# Perform the geometry analysis and export the TDT file of the surface
# geometry. Assume RA60 symmetry.
export_layout_to_tdt(
    fa_2,
    "fa_ra60_macro",
    tdt_setup=TdtSetup(
        type_geo=LayoutGeometryType.ROTATION,
        symmetry_type=SymmetryType.SIXTH,
        property_types=[PropertyType.MACRO, PropertyType.MATERIAL]
    )
)

# Print the execution time
print(f"--- Code executed in {time.time() - start_time} seconds ---")
