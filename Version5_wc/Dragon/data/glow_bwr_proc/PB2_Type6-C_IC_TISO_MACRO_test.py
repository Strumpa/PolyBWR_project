## Non-regression test for glow geometry generation capabilities :
# Goal : cover BWR controlled assemblies, channel box with rounded corners : TISO + MACRO for IC version
# Assembly represented : simplified Peach Bottom 2 Type 6 from 
# Ref: "BWR Progression Problems", C. Lawing, S. Palmtag, and M. Asgari, (2021),
# Oak Ridge National Laboratory (ORNL), Oak Ridge, TN (United States), 
# https://doi.org/10.2172/1838995

# Author : R. Guasch
# Date : 31-05-2026

## Glow imports
from glow import *
from glow.geometry_layouts.layouts import associate_colors_to_regions, build_compound_regions
from glow.support.types import *
from glow.geometry_layouts.cells import CartesianCell
from glow.geometry_layouts.layouts import Region
from glow.geometry_layouts.lattices import CartesianLattice
from glow.geometry_layouts.geometries import Rectangle, Circle
from glow.main import TdtSetup, export_layout_to_tdt
from glow.interface.geom_interface import *
from glow.support.types import GeometryType, LayoutGeometryType, PropertyType, \
    SymmetryType
from glow.interface.geom_entities import wrap_shape
import numpy as np
import os

SYM_TO_LAYOUT_AND_BOUNDARY = {
    SymmetryType.FULL: {"TISO": {"layout": LayoutGeometryType.ISOTROPIC, "boundary": None},
                        "TSPC": {"layout": LayoutGeometryType.RECTANGLE_SYM, "boundary": BoundaryType.AXIAL_SYMMETRY}},
    SymmetryType.DIAG: {"TISO": {"layout": LayoutGeometryType.SYMMETRIES_TWO, "boundary": BoundaryType.AXIAL_SYMMETRY},
                        "TSPC": {"layout": LayoutGeometryType.RECTANGLE_EIGHT, "boundary": BoundaryType.AXIAL_SYMMETRY}},
}


tdt_output_path = ".." # modify according to where you wish to save output tdt files
# check of output path exists and create if not
cwd = os.getcwd()
if not os.path.isabs(tdt_output_path):
    tdt_output_path = os.path.join(cwd, tdt_output_path)
if not os.path.exists(tdt_output_path):
    os.makedirs(tdt_output_path)


apply_diagonal_symmetry = True
### Geometric specifications : 
# fuel rod scale (ROD) : 
fuel_radius = 0.52070
gap_radius = 0.53213
clad_radius = 0.61341
pin_pitch = 1.6256

# water rod (WROD) 
water_rod_inner_radius = 0.67437
water_rod_outer_radius = 0.75057

# Lattice description (simplified to only 2 enrichments): 
# RODs are ordered in x-increasing, y-increasing order.
lattice_description = [
    ["ROD4", "ROD3", "ROD2", "ROD2", "ROD2", "ROD2", "ROD2", "ROD3"],
    ["ROD3", "ROD2", "ROD1", "ROD5", "ROD1", "ROD1", "ROD1", "ROD2"],
    ["ROD2", "ROD1", "ROD1", "ROD1", "ROD1", "ROD1", "ROD5", "ROD1"],
    ["ROD2", "ROD5", "ROD1", "ROD1", "WROD", "ROD1", "ROD1", "ROD1"],
    ["ROD2", "ROD1", "ROD1", "WROD", "ROD1", "ROD1", "ROD1", "ROD1"],
    ["ROD2", "ROD1", "ROD1", "ROD1", "ROD1", "ROD1", "ROD1", "ROD1"],
    ["ROD2", "ROD1", "ROD5", "ROD1", "ROD1", "ROD1", "ROD5", "ROD1"],
    ["ROD3", "ROD2", "ROD1", "ROD1", "ROD1", "ROD1", "ROD1", "ROD2"],
    ]

n_rows = 8
n_cols = 8

macro_groupings_2x2 = [
    ["MACRO01", "MACRO02", "MACRO02", "MACRO03", "MACRO03", "MACRO04", "MACRO04", "MACRO05"],
    ["MACRO06", "MACRO07", "MACRO07", "MACRO08", "MACRO08", "MACRO09", "MACRO09", "MACRO10"],
    ["MACRO06", "MACRO07", "MACRO07", "MACRO08", "MACRO08", "MACRO09", "MACRO09", "MACRO10"],
    ["MACRO11", "MACRO12", "MACRO12", "MACRO13", "MACRO13", "MACRO14", "MACRO14", "MACRO15"],
    ["MACRO11", "MACRO12", "MACRO12", "MACRO13", "MACRO13", "MACRO14", "MACRO14", "MACRO15"],
    ["MACRO16", "MACRO17", "MACRO17", "MACRO18", "MACRO18", "MACRO19", "MACRO19", "MACRO20"],
    ["MACRO16", "MACRO17", "MACRO17", "MACRO18", "MACRO18", "MACRO19", "MACRO19", "MACRO20"],
    ["MACRO21", "MACRO22", "MACRO22", "MACRO23", "MACRO23", "MACRO24", "MACRO24", "MACRO25"],
    ]

# Assembly scale : 
assembly_pitch = 15.24
# Some BWR fuel channels are not centered : the wide-wide corner leaves more room for control cross insertion 
gap_wide = 0.9017 
gap_narrow = 0.42418
channel_box_thickness = 0.254
inner_corner_radius_of_curvature = 0.9652
outer_corner_radius_of_curvature = inner_corner_radius_of_curvature + channel_box_thickness
Gd_rod_ids = ["ROD5"]
non_fuel_rod_ids = ["WROD"]
number_of_water_rods= 2

# Control cross position and dimensions : 
# To make use of the diagonal symmetry, the control cross is positionned at the wide-wide corner, selected as the "south-west" corner.
# The geometry presented in "BWR Progression Problems" is effectively rotated by pi/2 in the trigonometric direction.

cross_center = "south-west"
# control cross dimensions
number_tubes_per_wing = 21
blade_half_span = 12.3825 # length from cross center to tip, 
blade_thickness = 0.79248 # total thickness of the blade, needs to be divided by 2 to split along axis of symmetry
tip_radius = 0.39624 # radius of the rounded tip of the blade, optional, if not specified, the blade will have a sharp tip
central_structure_half_span = 1.98501 # length from center to begining of absorber tubes array,
sheath_thickness = 0.14224 # thickness of the sheath around the absorber tubes, distance from inner blade moderator to outer assembly moderator.
# absorber tube dimensions
absorber_tube_outer_radius = 0.23876
absorber_tube_inner_radius = 0.17526
absorber_material = "B4C"
sheath_material = "SS304"

## Compute useful dimensions to build assembly model :
channel_box_outer_side = assembly_pitch - gap_wide - gap_narrow
channel_box_inner_side = assembly_pitch - 2 * channel_box_thickness - gap_wide - gap_narrow 
intra_assembly_coolant_width = (channel_box_inner_side - 8 * pin_pitch) / 2.0
## Compute x / y offsets to switch from lattice to assembly coordinates :
offset_x = gap_wide + channel_box_thickness + intra_assembly_coolant_width
offset_y = offset_x

center=(assembly_pitch / 2.0, assembly_pitch / 2.0, 0.0)

##  Derived dimensions for control cross
inner_sheath_width = blade_thickness - 2 * sheath_thickness
wing_length = blade_half_span - central_structure_half_span

# Tube spacing: compute automatically if not provided.
# The tubes are distributed evenly inside the inner sheath
# region of the wing (from central structure edge + sheath
# thickness to wing tip - sheath thickness).
inner_wing_length = wing_length - sheath_thickness
tube_spacing = inner_wing_length / float(number_tubes_per_wing)

extra_moderator_gap = (tube_spacing - 2.0*absorber_tube_outer_radius) / 2.0
first_tube_offset = central_structure_half_span + absorber_tube_outer_radius + extra_moderator_gap


def computeVolumeBasedRadii(fuel_radius, gap_radius, clad_radius, fuel_volume_fractions):
    """
    Helper to define fuel region radii for fuel pins
    fuel_radius:
        radius of the fuel pellet, in cm
    gap_radius:
        radius of the expansion gap, in cm
    clad_radius:
        radius of the cladding, in cm
    fuel_volume_fractions:
        volume fractions for the fuel subdivisions, ordered radially
    """
    pin_radii = []
    if gap_radius is not None and gap_radius < fuel_radius:
        # in case expansion gap is inner most region, add it as the first radius and then compute the fuel radii based on the remaining fuel volume after subtracting the gap volume
        pin_radii.append(gap_radius)
    for volume_fract in fuel_volume_fractions:
        pin_radii.append(
            (volume_fract**0.5) * fuel_radius,
        )
    if gap_radius is not None and gap_radius > fuel_radius:
        # if gap radius is inner most region, add the clad radius as the last radius after the fuel radii
        pin_radii.append(gap_radius)    
    if clad_radius is not None:
        pin_radii.append(clad_radius)
    if gap_radius is not None and gap_radius > fuel_radius and clad_radius is None:
        ## print warning as this would be (fuel zones + gap zone) without clad. 
        print("Warning : gap radius is larger than fuel radius but no clad radius provided, this would lead to a pin with fuel zones and outer gap zone but no clad. Please check the input radii.")
    return pin_radii


def is_degenerate(geometry):
    """
    Check if a geometry is degenerate (zero area, empty, or None).

    Parameters
    ----------
    geometry : any
        Geometry object (Rectangle, Circle, or compound geometry)

    Returns
    -------
    bool
        True if geometry is None, has zero area, or is empty
    """
    if geometry is None:
        return True

    # Check for zero-area geometries
    # Rectangles have dimensions property
    if hasattr(geometry, 'dimensions'):
        dims = geometry.dimensions
        if dims[0] < 1e-10 or dims[1] < 1e-10:
            return True

    # Check for area attribute
    if hasattr(geometry, 'area'):
        if geometry.area < 1e-10:
            return True

    return False


def generate_macro_subregions(macro_rect, macro_name, coolant_rect, channel_box_rect):
    """
    Decompose a MACRO rectangle into subregions using set operations.

    For each MACRO rectangle, this function computes the intersection with subregion
    boundaries to produce subregion-specific geometries. The subregions are:
    - COOLANT: Inside coolant_rect
    - CHANNEL_BOX: Between coolant_rect and channel_box_rect
    - MODERATOR: Outside channel_box_rect

    Parameters
    ----------
    macro_rect : Rectangle
        Rectangle geometry for the MACRO region
    macro_name : str
        Identifier for the MACRO (e.g., "BOT_1", "CORNER_BL")
    coolant_rect : Rectangle
        Reference boundary for inner coolant region
    channel_box_rect : Rectangle
        Reference boundary for outer channel box region

    Returns
    -------
    dict
        Maps subregion type to geometry:
        {
            'COOLANT': geometry or None,
            'CHANNEL_BOX': geometry or None,
            'MODERATOR': geometry or None
        }

        Each geometry may be:
        - Rectangle: if subregion aligns with MACRO bounds
        - Wrapped compound geometry: if MACRO spans subregion boundary
        - None: if subregion doesn't exist for this MACRO
    """
    result = {}

    # Subregion 1: COOLANT = macro_rect ∩ coolant_rect (intersection)
    try:
        coolant_geom = macro_rect * coolant_rect  # intersection operator
        # Wrap geometry to handle COMPOUND types
        coolant_geom = wrap_shape(coolant_geom)
        if coolant_geom is not None and not is_degenerate(coolant_geom):
            result['COOLANT'] = coolant_geom
            print(f"  {macro_name}: COOLANT geometry computed")
    except Exception as e:
        print(f"  Warning: Failed to compute COOLANT subregion for {macro_name}: {e}")

    # Subregion 2: CHANNEL_BOX = (macro_rect ∩ channel_box_rect) - coolant_rect (intersection minus)
    try:
        temp_geom = macro_rect * channel_box_rect  # intersection
        channel_box_geom = temp_geom - coolant_rect  # difference
        # Wrap geometry to handle COMPOUND types
        channel_box_geom = wrap_shape(channel_box_geom)
        if channel_box_geom is not None and not is_degenerate(channel_box_geom):
            result['CHANNEL_BOX'] = channel_box_geom
            print(f"  {macro_name}: CHANNEL_BOX geometry computed")
    except Exception as e:
        print(f"  Warning: Failed to compute CHANNEL_BOX subregion for {macro_name}: {e}")

    # Subregion 3: MODERATOR = macro_rect - channel_box_rect (difference)
    try:
        moderator_geom = macro_rect - channel_box_rect  # difference operator
        # Wrap geometry to handle COMPOUND types
        moderator_geom = wrap_shape(moderator_geom)
        if moderator_geom is not None and not is_degenerate(moderator_geom):
            result['MODERATOR'] = moderator_geom
            print(f"  {macro_name}: MODERATOR geometry computed")
    except Exception as e:
        print(f"  Warning: Failed to compute MODERATOR subregion for {macro_name}: {e}")

    return result


def generate_corner_channel_box_regions(x0, y0, x1, y1, ap, n_rows, n_cols, coolant_rect, channel_box_rect, center):
    """
    Generate channel box material regions at assembly corners that intersect with lattice footprint.

    These are the 4 corner regions where:
    - The channel box extends into the lattice footprint
    - They are OUTSIDE coolant_rect, INSIDE channel_box_rect, INSIDE lattice footprint
    - They should be assigned to the corner fuel pin macros (MACRO00, MACRO0{n_cols-1}, etc.)

    Parameters
    ----------
    x0, y0, x1, y1 : float
        Lattice footprint bounds
    ap : float
        Assembly pitch
    n_rows, n_cols : int
        Number of rows and columns in lattice
    coolant_rect : Rectangle
        Reference boundary for inner coolant region
    channel_box_rect : Rectangle
        Reference boundary for outer channel box region
    center : tuple
        Assembly center (x_center, y_center, z_center)

    Returns
    -------
    dict
        Maps corner identifier to dict with geometry and macro name:
        {
            'BL': {'geometry': Rectangle or None, 'macro': 'MACRO00'},
            'BR': {'geometry': Rectangle or None, 'macro': 'MACRO0{n_cols-1}'},
            'TL': {'geometry': Rectangle or None, 'macro': 'MACRO{n_rows-1}0'},
            'TR': {'geometry': Rectangle or None, 'macro': 'MACRO{n_rows-1}{n_cols-1}'}
        }
    """
    z_center = center[2]
    result = {}
    pin_pitch_x = (x1 - x0) / n_cols  # Assuming uniform pin pitch in x
    pin_pitch_y = (y1 - y0) / n_rows  # Assuming uniform pin pitch in y
    # Define the 4 corner regions in assembly space
    corners = {
        'BL': {'x_min': x0, 'x_max': x0 + pin_pitch_x, 'y_min': y0, 'y_max': y0 + pin_pitch_y, 'macro_row': 0, 'macro_col': 0},
        'BR': {'x_min': x1 - pin_pitch_x, 'x_max': x1, 'y_min': y0, 'y_max': y0 + pin_pitch_y, 'macro_row': 0, 'macro_col': n_cols - 1},
        'TL': {'x_min': x0, 'x_max': x0 + pin_pitch_x, 'y_min': y1 - pin_pitch_y, 'y_max': y1, 'macro_row': n_rows - 1, 'macro_col': 0},
        'TR': {'x_min': x1 - pin_pitch_x, 'x_max': x1, 'y_min': y1 - pin_pitch_y, 'y_max': y1, 'macro_row': n_rows - 1, 'macro_col': n_cols - 1}
    }

    for corner_id, corner_def in corners.items():
        x_min, x_max = corner_def['x_min'], corner_def['x_max']
        y_min, y_max = corner_def['y_min'], corner_def['y_max']
        macro_row = corner_def['macro_row']
        macro_col = corner_def['macro_col']

        try:
            # Create corner rectangle in assembly coordinates
            width = x_max - x_min
            height = y_max - y_min
            corner_center_x = (x_min + x_max) / 2.0
            corner_center_y = (y_min + y_max) / 2.0

            corner_rect = Rectangle(
                name=f"CORNER_CB_{corner_id}",
                height=height,
                width=width,
                center=(corner_center_x, corner_center_y, z_center)
            )

            # Compute channel box material portion: corner_rect ∩ channel_box_rect - coolant_rect
            temp_geom = corner_rect * channel_box_rect  # intersection
            corner_cb_geom = temp_geom - coolant_rect   # difference (remove coolant)

            # Wrap and validate
            corner_cb_geom = wrap_shape(corner_cb_geom)

            if corner_cb_geom is not None and not is_degenerate(corner_cb_geom):
                macro_name = macro_groupings_2x2[macro_row][macro_col]
                
                result[corner_id] = {
                    'geometry': corner_cb_geom,
                    'macro': macro_name
                }
            else:
                result[corner_id] = {
                    'geometry': None,
                    'macro': macro_name
                }

        except Exception as e:
            print(f"  Warning: Failed to compute channel box region for corner {corner_id}: {e}")
            macro_name = macro_groupings_2x2[macro_row][macro_col]
            result[corner_id] = {
                'geometry': None,
                'macro': macro_name
            }

    return result


def _corner_transform(corner, x, y, ap):
    """
    Map canonical north-west coordinates ``(x, y)`` to the actual
    assembly corner.

    In the canonical (north-west) system the cross centre sits at the
    top-left corner ``(0, ap)``.  This helper mirrors the coordinates
    for the other three corners.

    Parameters
    ----------
    corner : str
        ``"north-west"``, ``"north-east"``, ``"south-west"``, or
        ``"south-east"``.
    x, y : float
        Coordinates in the canonical north-west system.
    ap : float
        Assembly pitch.

    Returns
    -------
    (float, float)
        Transformed ``(x, y)``.
    """
    if corner == "north-west":
        return (x, y)
    elif corner == "north-east":
        return (ap - x, y)
    elif corner == "south-west":
        return (x, ap - y)
    elif corner == "south-east":
        return (ap - x, ap - y)
    else:
        raise ValueError(f"Unknown corner '{corner}'.")


def _remap_rounded_corner_indices(corner_indices, cross_corner):
    """
    Remap glow rounded-corner indices for a wing-tip rectangle after
    applying the corner transform.

    In the canonical north-west system the wing tips have a rounded
    corner at glow index 1 (= bottom-right of the rectangle in glow
    convention).  When the cross is placed at a different assembly
    corner the rectangle is mirrored and the rounded-corner index
    must change accordingly.

    Parameters
    ----------
    corner_indices : list of (int, float)
        List of ``(glow_corner_index, radius)`` pairs in the canonical
        system (north-west).
    cross_corner : str
        The actual assembly corner.

    Returns
    -------
    list of (int, float)
        Remapped corner/radius pairs.
    """
    # Mapping: NW is identity.  Mirror in x flips left <-> right (0 <-> 1, 3 <-> 2).
    # Mirror in y flips top <-> bottom (0 <-> 3, 1 <-> 2).
    _mirror_x = {0: 1, 1: 0, 2: 3, 3: 2}
    _mirror_y = {0: 3, 1: 2, 2: 1, 3: 0}

    result = list(corner_indices)
    if cross_corner in ("north-east", "south-east"):
        result = [(_mirror_x[idx], r) for idx, r in result]
    if cross_corner in ("south-west", "south-east"):
        result = [(_mirror_y[idx], r) for idx, r in result]
    return result


def build_control_cross_elements(ap, assembly_center=(None, None, None)):
    """
    Build the elements describing the control cross
    
    Parameters
    ----------
    ap : float
        Assembly pitch.
    assembly_center : tuple, optional
        (x, y, z) center of the assembly_universe. Defaults to (ap/2, ap/2, 0).
        Region centers are converted from absolute assembly coords to local
        (center-relative) coords.


    Returns 
    -------
    control_cross_universe : CartesianCell
        A CartesianCell containing all the Regions that define the control cross geometry, centered at 
        either the NW or SW corner depending on the identified lattice symmetry.

    """

    # retrieve dimensions and properties of the control cross model

    if assembly_center[0] is None:
        assembly_center = (ap / 2.0, ap / 2.0, 0.0)

    absorber_tube_inner_radius
    extra_moderator_gap = (tube_spacing - 2.0*absorber_tube_outer_radius) / 2.0
    half_square_side = absorber_tube_outer_radius + extra_moderator_gap
    distance_from_first_center_to_last_center = (number_tubes_per_wing - 1) * 2.0 * half_square_side
    width_controlled_section = blade_half_span - sheath_thickness - central_structure_half_span

    first_tube_offset = central_structure_half_span + absorber_tube_outer_radius + extra_moderator_gap # In assembly coordinates.

    assembly_bounding_rect = Rectangle(
        height=ap,
        width=ap,
        center=(ap/2,ap/2,0.0)
    )

    elements = {}
    elements["CTRL_H_SHEATH"] = {}
    elements["CTRL_V_SHEATH"] = {}
    # Create a universe cell for the control cross and add all regions
    cx, cy = _corner_transform(cross_center, blade_half_span/2.0, ap, ap)
    control_cross_horizontal_wing = Rectangle(
        name="CTRL_CROSS_H",
        width=blade_half_span,
        height=blade_thickness,
        center=(cx, cy, ap, 0.0)
        #rounded_corners=_remap_rounded_corner_indices([(1, tip_radius), (2, tip_radius)], cross_center) if tip_radius > 0.0 else None
    )
    control_cross_horizontal_wing_rounded = Rectangle(
        name="CTRL_CROSS_H_rounded",
        width=blade_half_span,
        height=blade_thickness,
        center=(cx, cy, ap, 0.0),
        rounded_corners=_remap_rounded_corner_indices([(1, tip_radius), (2, tip_radius)], cross_center) if tip_radius > 0.0 else None
    )

    # compute delta between centers of horizontal wing with rounded corners and horizontal wing without rounded corners
    delta_centers_x_horizontal_wing = get_point_coordinates(control_cross_horizontal_wing_rounded.o)[0] - get_point_coordinates(control_cross_horizontal_wing.o)[0]
    delta_centers_y_horizontal_wing = get_point_coordinates(control_cross_horizontal_wing_rounded.o)[1] - get_point_coordinates(control_cross_horizontal_wing.o)[1]
    #wing_elements[(blade_half_span/2.0 + delta_centers_x_horizontal_wing, ap + delta_centers_y_horizontal_wing, 0.0)] = control_cross_horizontal_wing_rounded
    
    h_wing_in_assembly_footprint = control_cross_horizontal_wing_rounded * assembly_bounding_rect

    elements["CTRL_H_SHEATH"]["geometry"] = h_wing_in_assembly_footprint
    elements["CTRL_H_SHEATH"]["offset"] = (delta_centers_x_horizontal_wing, delta_centers_y_horizontal_wing, 0.0)
    elements["CTRL_H_SHEATH"]["macro"] = "CTRL_H"
    elements["CTRL_H_SHEATH"]["material"] = sheath_material
    print(f"offset is : {elements['CTRL_H_SHEATH']['offset']}")

    cx, cy = _corner_transform(cross_center, 0.0, ap - blade_half_span/2.0, ap)
    control_cross_vertical_wing = Rectangle(
        name="CTRL_CROSS_V",
        width=blade_thickness,
        height=blade_half_span,
        center=(cx, cy, 0.0)
    )
    control_cross_vertical_wing_rounded = Rectangle(
        name="CTRL_CROSS_V_rounded",
        width=blade_thickness,
        height=blade_half_span,
        center=(cx, cy, 0.0),
        rounded_corners=_remap_rounded_corner_indices([(0, tip_radius), (1, tip_radius)], cross_center) if tip_radius > 0.0 else None
    )

    delta_centers_x_vertical_wing = get_point_coordinates(control_cross_vertical_wing_rounded.o)[0] - get_point_coordinates(control_cross_vertical_wing.o)[0]
    delta_centers_y_vertical_wing = get_point_coordinates(control_cross_vertical_wing_rounded.o)[1] - get_point_coordinates(control_cross_vertical_wing.o)[1]
    
    v_wing_in_assembly_footprint = control_cross_vertical_wing_rounded * assembly_bounding_rect
    elements["CTRL_V_SHEATH"]["geometry"] = v_wing_in_assembly_footprint
    elements["CTRL_V_SHEATH"]["offset"] = (delta_centers_x_vertical_wing, delta_centers_y_vertical_wing, 0.0)
    elements["CTRL_V_SHEATH"]["macro"] = "CTRL_V"
    elements["CTRL_V_SHEATH"]["material"] = sheath_material

    ### Corner at intersection between H and V wings
    
    elements["CTRL_CORNER"] = {}
    corner = control_cross_vertical_wing * control_cross_horizontal_wing
    corner_in_assembly_footprint = corner * assembly_bounding_rect
    elements["CTRL_CORNER"]["geometry"] = corner_in_assembly_footprint
    elements["CTRL_CORNER"]["macro"] = "CTRL_CORNER"
    elements["CTRL_CORNER"]["material"] = sheath_material

    # build hollow sheath region and absorber pins 
    base_material = "MODERATOR"
    elements["CTRL_H_HOLLOW"] = {}
    elements["CTRL_V_HOLLOW"] = {}
    # build inner region for hollow sheath (moderator-filled cavity inside the blade)
    cx, cy = _corner_transform(cross_center, central_structure_half_span + (blade_half_span - sheath_thickness - central_structure_half_span)/2.0, ap, ap)
    inner_sheath_horizontal_wing = Rectangle(
        name="CTRL_CROSS_H_INNER",
        width=(blade_half_span - sheath_thickness - central_structure_half_span),
        height=(blade_thickness - 2.0*sheath_thickness),
        center=(cx, cy, 0.0)
    )
    # get center for the horizontal wing inner sheath without rounded corners :
    inner_sheath_horizontal_wing_rounded = Rectangle(
        name="CTRL_CROSS_H_INNER_HOLLOW",
        width=(blade_half_span - sheath_thickness - central_structure_half_span),
        height=(blade_thickness - 2.0*sheath_thickness),
        center=(cx, cy, 0.0),
        rounded_corners=_remap_rounded_corner_indices([(1, tip_radius - sheath_thickness), (2, tip_radius - sheath_thickness)], cross_center) if tip_radius > 0.0 else None
    )
    delta_centers_x_inner_sheath = get_point_coordinates(inner_sheath_horizontal_wing_rounded.o)[0] - get_point_coordinates(inner_sheath_horizontal_wing.o)[0]
    delta_centers_y_inner_sheath = get_point_coordinates(inner_sheath_horizontal_wing_rounded.o)[1] - get_point_coordinates(inner_sheath_horizontal_wing.o)[1]
    
    intersection_with_assembly_footprint = inner_sheath_horizontal_wing_rounded * assembly_bounding_rect
    elements["CTRL_H_HOLLOW"]["geometry"] = intersection_with_assembly_footprint
    elements["CTRL_H_HOLLOW"]["offset"] = (delta_centers_x_inner_sheath, delta_centers_y_inner_sheath, 0.0)
    elements["CTRL_H_HOLLOW"]["macro"] = "CTRL_H"
    elements["CTRL_H_HOLLOW"]["material"] = base_material


    # inner region for hollow sheath in vertical wing :
    cx, cy = _corner_transform(cross_center, 0.0, ap - (central_structure_half_span + (blade_half_span - sheath_thickness - central_structure_half_span)/2.0), ap)
    inner_sheath_vertical_wing = Rectangle(
        name="CTRL_CROSS_V_INNER",
        width=(blade_thickness - 2.0*sheath_thickness),
        height=(blade_half_span - sheath_thickness - central_structure_half_span),
        center=(cx, cy, 0.0)
    )
    inner_sheath_vertical_wing_rounded = Rectangle(
        name="CTRL_CROSS_V_INNER_HOLLOW",
        width=(blade_thickness - 2.0*sheath_thickness),
        height=(blade_half_span - sheath_thickness - central_structure_half_span),
        center=(cx, cy, 0.0),
        rounded_corners=_remap_rounded_corner_indices([(0, tip_radius - sheath_thickness), (1, tip_radius - sheath_thickness)], cross_center) if tip_radius > 0.0 else None
    )

    delta_centers_x_inner_sheath_vertical = get_point_coordinates(inner_sheath_vertical_wing_rounded.o)[0] - get_point_coordinates(inner_sheath_vertical_wing.o)[0]
    delta_centers_y_inner_sheath_vertical = get_point_coordinates(inner_sheath_vertical_wing_rounded.o)[1] - get_point_coordinates(inner_sheath_vertical_wing.o)[1]
    
    intersection_with_assembly_footprint = inner_sheath_vertical_wing_rounded * assembly_bounding_rect
    elements["CTRL_V_HOLLOW"]["geometry"] = intersection_with_assembly_footprint
    elements["CTRL_V_HOLLOW"]["offset"] = (delta_centers_x_inner_sheath_vertical, delta_centers_y_inner_sheath_vertical, 0.0)
    elements["CTRL_V_HOLLOW"]["macro"] = "CTRL_V"
    elements["CTRL_V_HOLLOW"]["material"] = base_material
    
    for i in range(number_tubes_per_wing):
        elements[f"ABS_TUBE_H_{i}"] = {}
        elements[f"ABS_TUBE_V_{i}"] = {}
        offset = first_tube_offset + float(i * tube_spacing)
        # Horizontal wing tube (along x)
        tx_h, ty_h = _corner_transform(cross_center, offset, ap, ap)
        tx_v, ty_v = _corner_transform(cross_center, 0.0, ap - offset, ap)
        # Hollow tube (GE-14 style): inner absorber + outer cladding
        abs_tube_h = Circle(radius=absorber_tube_inner_radius, center=(tx_h, ty_h, 0.0))
        h_tube_intersection_with_af = abs_tube_h * assembly_bounding_rect
        elements[f"ABS_TUBE_H_{i}"]["geometry"] = h_tube_intersection_with_af
        elements[f"ABS_TUBE_H_{i}"]["macro"] = "CTRL_H"
        elements[f"ABS_TUBE_H_{i}"]["material"] = absorber_material

        abs_tube_v = Circle(radius=absorber_tube_inner_radius, center=(tx_v, ty_v, 0.0))
        v_tube_intersection_with_af = abs_tube_v * assembly_bounding_rect
        elements[f"ABS_TUBE_V_{i}"]["geometry"] = v_tube_intersection_with_af
        elements[f"ABS_TUBE_V_{i}"]["macro"] = "CTRL_V"
        elements[f"ABS_TUBE_V_{i}"]["material"] = absorber_material

        elements[f"SHEATH_TUBE_H_{i}"] = {}
        elements[f"SHEATH_TUBE_V_{i}"] = {}
        sheath_tube_h = Circle(radius=absorber_tube_outer_radius, center=(tx_h, ty_h, 0.0))
        h_tube_intersection_with_af = sheath_tube_h * assembly_bounding_rect - abs_tube_h
        elements[f"SHEATH_TUBE_H_{i}"]["geometry"] = h_tube_intersection_with_af
        elements[f"SHEATH_TUBE_H_{i}"]["macro"] = "CTRL_H"
        elements[f"SHEATH_TUBE_H_{i}"]["material"] = sheath_material

        sheath_tube_v = Circle(radius=absorber_tube_outer_radius, center=(tx_v, ty_v, 0.0))
        v_tube_intersection_with_af = sheath_tube_v * assembly_bounding_rect - abs_tube_v
        elements[f"SHEATH_TUBE_V_{i}"]["geometry"] = v_tube_intersection_with_af
        elements[f"SHEATH_TUBE_V_{i}"]["macro"] = "CTRL_V"
        elements[f"SHEATH_TUBE_V_{i}"]["material"] = sheath_material

    cx, cy = _corner_transform(cross_center, blade_half_span - sheath_thickness/2.0, ap - blade_thickness/4.0, ap)
    tip_rectangle_h = Rectangle(
        name="CTRL_H_TIP",
        width=sheath_thickness,
        height=blade_thickness/2,
        center=(cx, cy, 0.0)
    )
    rounded_modertator_tip_h = control_cross_horizontal_wing*assembly_bounding_rect - control_cross_horizontal_wing_rounded*assembly_bounding_rect
    rounded_sheath_tip = tip_rectangle_h*h_wing_in_assembly_footprint
    moder_around_tip = tip_rectangle_h - h_wing_in_assembly_footprint
    extra_moder_outside_of_tip_rectangle = rounded_modertator_tip_h - tip_rectangle_h
    # build regions information and store in elements
    macro = "CTRL_H"
    elements["CTRL_H_TIP_SHEATH"] = {}
    elements["CTRL_H_TIP_SHEATH"]["macro"] = macro
    elements["CTRL_H_TIP_SHEATH"]["material"] = sheath_material
    elements["CTRL_H_TIP_SHEATH"]["geometry"] = rounded_sheath_tip

    elements["CTRL_H_TIP_MODER"] = {}
    elements["CTRL_H_TIP_MODER"]["macro"] = macro 
    elements["CTRL_H_TIP_MODER"]["material"] = "MODERATOR"
    elements["CTRL_H_TIP_MODER"]["geometry"] = moder_around_tip

    elements["CTRL_H_TIP_EXTRA_MODER"] = {}
    elements["CTRL_H_TIP_EXTRA_MODER"]["macro"] = macro 
    elements["CTRL_H_TIP_EXTRA_MODER"]["material"] = "MODERATOR"
    elements["CTRL_H_TIP_EXTRA_MODER"]["geometry"] = extra_moder_outside_of_tip_rectangle

    cx, cy = _corner_transform(cross_center, blade_thickness/4.0, ap - (blade_half_span - sheath_thickness/2.0), ap)
    tip_rectangle_v = Rectangle(
        name="CTRL_V_TIP",
        height=sheath_thickness,
        width=blade_thickness/2,
        center=(cx, cy, 0.0)
    )
    rounded_modertator_tip_v = control_cross_vertical_wing*assembly_bounding_rect - control_cross_vertical_wing_rounded*assembly_bounding_rect
    rounded_sheath_tip = tip_rectangle_v*v_wing_in_assembly_footprint
    moder_around_tip = tip_rectangle_v - v_wing_in_assembly_footprint
    extra_moder_outside_of_tip_rectangle = rounded_modertator_tip_v - tip_rectangle_v

    # build regions information and store in elements
    macro = "CTRL_V"
    elements["CTRL_V_TIP_SHEATH"] = {}
    elements["CTRL_V_TIP_SHEATH"]["macro"] = macro
    elements["CTRL_V_TIP_SHEATH"]["material"] = sheath_material
    elements["CTRL_V_TIP_SHEATH"]["geometry"] = rounded_sheath_tip

    elements["CTRL_V_TIP_MODER"] = {}
    elements["CTRL_V_TIP_MODER"]["macro"] = macro 
    elements["CTRL_V_TIP_MODER"]["material"] = "MODERATOR"
    elements["CTRL_V_TIP_MODER"]["geometry"] = moder_around_tip

    elements["CTRL_V_TIP_EXTRA_MODER"] = {}
    elements["CTRL_V_TIP_EXTRA_MODER"]["macro"] = macro 
    elements["CTRL_V_TIP_EXTRA_MODER"]["material"] = "MODERATOR"
    elements["CTRL_V_TIP_EXTRA_MODER"]["geometry"] = extra_moder_outside_of_tip_rectangle



    ### Introduce rectangular macros in the prolongation of cross arms
    # horizontal :
    cx, cy = _corner_transform(cross_center, blade_half_span + ((n_cols-1)*pin_pitch+offset_x - blade_half_span)/2, ap - blade_thickness/4.0, ap)
    prolong_rect_h = Rectangle(
        name = "prolong_H",
        height=blade_thickness/2,
        width=((n_cols-1)*pin_pitch+offset_x - blade_half_span),
        center=(cx, cy, 0.0)
    )
    elements["CTRL_PROLONG_H"] = {}
    elements["CTRL_PROLONG_H"]["macro"] = "CTRL_H"
    elements["CTRL_PROLONG_H"]["material"] = "MODERATOR"
    elements["CTRL_PROLONG_H"]["geometry"] = prolong_rect_h

    # vertical :
    cx, cy = _corner_transform(cross_center, blade_thickness/4.0, ap - (blade_half_span + ((n_cols-1)*pin_pitch+offset_y - blade_half_span)/2), ap)
    prolong_rect_v = Rectangle(
        name = "prolong_V",
        width=blade_thickness/2,
        height=((n_cols-1)*pin_pitch+offset_y - blade_half_span),
        center=(cx, cy, 0.0)
    )
    elements["CTRL_PROLONG_V"] = {}
    elements["CTRL_PROLONG_V"]["macro"] = "CTRL_V"
    elements["CTRL_PROLONG_V"]["material"] = "MODERATOR"
    elements["CTRL_PROLONG_V"]["geometry"] = prolong_rect_v

    return elements

# Create an empty lattice
lattice = CartesianLattice(
    name=f"PB2_Type6-C_lattice",
    centre=center,
    cells=[],
)

# Iterate over lattice description to create all fuel cells, water rod cells and add to a Lattice
water_rod_counter = 0
row_idx = -1
for row in lattice_description:
    row_idx += 1
    cell_idx = -1
    for cell in row:
        cell_idx += 1
        macro_name = macro_groupings_2x2[row_idx][cell_idx]
        if cell in non_fuel_rod_ids:
            ## Build water rod geometry
            water_rod_counter += 1
            tmp_cell = CartesianCell(
                    name=f"WROD_{water_rod_counter}",
                    width_height=(
                        pin_pitch,
                        pin_pitch,
                ),
                center=(0.0, 0.0, 0.0),
                base_props={PropertyType.MATERIAL: "COOLANT",
                            PropertyType.MACRO: macro_name},
            )
            tmp_cell.add(Region(Circle(radius=water_rod_outer_radius), 
                                properties={PropertyType.MATERIAL: "CLAD",
                                            PropertyType.MACRO: macro_name}))
            tmp_cell.add(Region(Circle(radius=water_rod_inner_radius), 
                                properties={PropertyType.MATERIAL: "MODERATOR",
                                            PropertyType.MACRO: macro_name}))
        else:
            if cell in Gd_rod_ids:
                fuel_material_names = [f"Gd_{cell}", f"Gd_{cell}", f"Gd_{cell}",f"Gd_{cell}", f"Gd_{cell}", f"Gd_{cell}"]
                fuel_volume_fractions = [0.2, 0.4, 0.6, 0.8, 0.95, 1.0]
            else:
                fuel_material_names = [f"UOX_{cell}", f"UOX_{cell}", f"UOX_{cell}", f"UOX_{cell}"]
                fuel_volume_fractions = [0.5, 0.8, 0.95, 1.0]

            radii = computeVolumeBasedRadii(fuel_radius, gap_radius, clad_radius, fuel_volume_fractions)
            list_of_material_names = [fuel_material_name for fuel_material_name in fuel_material_names] + ["GAP", "CLAD"]
            tmp_cell = CartesianCell(
                    name=fuel_material_names[0],
                    width_height=(pin_pitch, pin_pitch),
                    center=(0.0, 0.0, 0.0),
                    base_props={PropertyType.MATERIAL:"COOLANT",
                                PropertyType.MACRO: macro_name}
                )
            for radius, mat in zip(radii[::-1],list_of_material_names[::-1]):
                tmp_cell.add(
                    Region(Circle(radius=radius), properties={PropertyType.MATERIAL:mat, 
                                                                PropertyType.MACRO:macro_name})
            )
        lattice.add(
                tmp_cell, position=((cell_idx + 0.5) * pin_pitch + offset_x,
                                (row_idx + 0.5) * pin_pitch + offset_y,
                                0.0)
                )
        

## Create a base CartesianCell to hold the assembly model :
assembly = CartesianCell(
    name="PB2_Type6-C_assembly",
    width_height=(assembly_pitch, assembly_pitch),
    center=center,
    base_props={PropertyType.MATERIAL: "COOLANT",
                PropertyType.MACRO:"BASE_CELL"},
)

assembly.add(lattice)

### Compute lattice footprint (square bounding box)
x0 = offset_x
y0 = offset_y
x1 = x0 + n_cols*pin_pitch
y1 = y0 + n_rows*pin_pitch


### 

### Generate rectangles for all MACRO regions surrounding the pin lattice.
macro_rects = {}
z_center = center[2]

# ======================================================================
# BOTTOM STRIP: [x0, x1] × [0, y0]
# Subdivided into n_cols regions (one per pin column)
# ======================================================================
strip_height = y0  # Distance from bottom edge (y=0) to lattice bottom (y=y0)

for col_idx in range(n_cols):
    region_name_key = f"BOT_{col_idx}"
    # Rectangle spans one pin column width
    col_x_min = x0 + col_idx * pin_pitch
    col_x_max = col_x_min + pin_pitch

    # Center of this bottom strip rectangle
    rect_center_x = (col_x_min + col_x_max) / 2.0
    rect_center_y = strip_height / 2.0  # Centered in the bottom strip

    macro_name = macro_groupings_2x2[0][col_idx]
    if col_idx == 0:
        #macro_name = "TOTO"
        print("In region of bottom strip that should be macro00")
    rect_params = (
        macro_name,                                  # name
        strip_height,                                # height (y-direction)
        pin_pitch,                                   # width (x-direction)
        (rect_center_x, rect_center_y, z_center),   # center (x, y, z)
        None                                         # rounded_corners (None = no rounding)
    )
    macro_rects[region_name_key] = rect_params

# ======================================================================
# TOP STRIP: [x0, x1] × [y1, assembly_pitch]
# Subdivided into n_cols regions (one per pin column)
# ======================================================================
strip_height = assembly_pitch - y1  # Distance from lattice top (y=y1) to top edge (y=assembly_pitch)

for col_idx in range(n_cols):
    region_name_key = f"TOP_{col_idx}"
    col_x_min = x0 + col_idx * pin_pitch
    col_x_max = col_x_min + pin_pitch

    rect_center_x = (col_x_min + col_x_max) / 2.0
    rect_center_y = y1 + (assembly_pitch - y1) / 2.0  # Centered in the top strip

    macro_name = macro_groupings_2x2[-1][col_idx]

    rect_params = (
        macro_name,
        strip_height,
        pin_pitch,
        (rect_center_x, rect_center_y, z_center),
        None                                         # rounded_corners (None = no rounding)
    )
    macro_rects[region_name_key] = rect_params

# ======================================================================
# LEFT STRIP: [0, x0] × [y0, y1]
# Subdivided into n_rows regions (one per pin row)
# ======================================================================
strip_width = x0  # Distance from left edge (x=0) to lattice left (x=x0)

for row_idx in range(n_rows):
    region_name_key = f"LEFT_{row_idx}"
    row_y_min = y0 + row_idx * pin_pitch
    row_y_max = row_y_min + pin_pitch

    rect_center_x = strip_width / 2.0  # Centered in the left strip
    rect_center_y = (row_y_min + row_y_max) / 2.0
    macro_name = macro_groupings_2x2[row_idx][0]

    rect_params = (
        macro_name,
        pin_pitch,                                   # height (y-direction)
        strip_width,                                 # width (x-direction)
        (rect_center_x, rect_center_y, z_center),
        None                                         # rounded_corners (None = no rounding)
    )
    macro_rects[region_name_key] = rect_params

# ======================================================================
# RIGHT STRIP: [x1, assembly_pitch] × [y0, y1]
# Subdivided into n_rows regions (one per pin row)
# ======================================================================
strip_width = assembly_pitch - x1  # Distance from lattice right (x=x1) to right edge (x=assembly_pitch)

for row_idx in range(n_rows):
    region_name_key = f"RIGHT_{row_idx}"
    row_y_min = y0 + row_idx * pin_pitch
    row_y_max = row_y_min + pin_pitch

    rect_center_x = x1 + (assembly_pitch - x1) / 2.0  # Centered in the right strip
    rect_center_y = (row_y_min + row_y_max) / 2.0
    macro_name = macro_groupings_2x2[row_idx][-1]

    rect_params = (
        macro_name,
        pin_pitch,
        strip_width,
        (rect_center_x, rect_center_y, z_center),
        None                                         # rounded_corners (None = no rounding)
    )
    macro_rects[region_name_key] = rect_params

# ======================================================================
# CORNERS
# ======================================================================
corner_width_bl = x0
corner_height_bl = y0
macro_name = macro_groupings_2x2[0][0]

# Bottom-Left
macro_rects["CORNER_BL"] = (
    macro_name,
    corner_height_bl,
    corner_width_bl,
    (corner_width_bl / 2.0, corner_height_bl / 2.0, z_center),
    None                                             # rounded_corners (None = no rounding)
)

# Bottom-Right
corner_width_br = assembly_pitch - x1
corner_height_br = y0
macro_name = macro_groupings_2x2[0][-1]

macro_rects["CORNER_BR"] = (
    macro_name,
    corner_height_br,
    corner_width_br,
    (x1 + corner_width_br / 2.0, corner_height_br / 2.0, z_center),
    None                                             # rounded_corners (None = no rounding)
)

# Top-Left
corner_width_tl = x0
corner_height_tl = assembly_pitch - y1
macro_name = macro_groupings_2x2[-1][0]

macro_rects["CORNER_TL"] = (
    macro_name,
    corner_height_tl,
    corner_width_tl,
    (corner_width_tl / 2.0, y1 + corner_height_tl / 2.0, z_center),
    None
)

# Top-Right
corner_width_tr = assembly_pitch - x1
corner_height_tr = assembly_pitch - y1
macro_name = macro_groupings_2x2[-1][-1]

macro_rects["CORNER_TR"] = (
    macro_name,
    corner_height_tr,
    corner_width_tr,
    (x1 + corner_width_tr / 2.0, y1 + corner_height_tr / 2.0, z_center),
    None
)

### Create lattice, coolant and channel box bounding rectangles, with rounded corners

offset = (gap_wide - gap_narrow) / 2.0

rect_center = (center[0] + offset, center[1] + offset, center[2])


# Build rectangles with potentially asymmetric dimensions
coolant_rect = Rectangle(
    name="intra_assembly_coolant",
    height=channel_box_inner_side,
    width=channel_box_inner_side,
    center=rect_center,
    rounded_corners=[
        (0,inner_corner_radius_of_curvature),
        (1,inner_corner_radius_of_curvature),
        (2,inner_corner_radius_of_curvature),
        (3,inner_corner_radius_of_curvature),
    ],
)

channel_box_rect = Rectangle(
    name="channel_box",
    height=channel_box_outer_side,
    width=channel_box_outer_side,
    center=rect_center,
    rounded_corners=[
        (0,outer_corner_radius_of_curvature),
        (1,outer_corner_radius_of_curvature),
        (2,outer_corner_radius_of_curvature),
        (3,outer_corner_radius_of_curvature),
    ],
)

lattice_rect = Rectangle(
    name="lattice_bounding_rectangular_box",
    height=n_rows*pin_pitch,
    width=n_cols*pin_pitch,
    center=rect_center
)



for name, rect_params in macro_rects.items():
    macro_name, height, width, rect_center, rounded = rect_params
    print("Creating regions for macro_rect with :")
    print(f"name : {name}")
    print(f"dimensions [h x w] = {height} x {width}")
    print(f"center : {rect_center}")
    # Create Rectangle geometry for this MACRO
    macro_rect = Rectangle(
        name=name,
        height=height,
        width=width,
        center=rect_center,  # Keep center for boolean operations to work correctly
        rounded_corners=rounded
    )

    # Generate sub-regions using set operations
    subregion_geometries = generate_macro_subregions(
        macro_rect,
        macro_name,
        coolant_rect,
        channel_box_rect
    )

    # Create Region objects for each subregion portion
    for material_type, geometry in subregion_geometries.items():
        if geometry is None:  # Skip if subregion doesn't exist for this MACRO
            print(f"WARNING : None geometry for material_type : {material_type}")
            continue

        try:
            # Create Region
            # All portions of a MACRO (e.g., BOT_1_COOLANT and BOT_1_CHANNEL_BOX)
            # share the same MACRO property identifier (BOT_1)
            region = Region(
                geometry,
                properties={
                    PropertyType.MATERIAL: material_type,
                    PropertyType.MACRO: macro_name
                }
            )
            subregion_center = get_point_coordinates(region.o)
            assembly.add(region, position=subregion_center)

        except Exception as e:
            print(f"  Warning: Failed to create Region for {macro_name}/{material_type}: {e}")


### Generate corner channel box regions
# These regions bridge between lattice footprint and channel box at corners,
# assigning them to the corner fuel pin macros (MACRO00, MACRO0{n_cols-1}, etc.)
corner_cb_regions = generate_corner_channel_box_regions(
    x0, y0, x1, y1, assembly_pitch,
    n_rows, n_cols,
    coolant_rect, channel_box_rect,
    center
)

for corner_id, corner_data in corner_cb_regions.items():
    geometry = corner_data['geometry']
    macro_name = corner_data['macro']

    if geometry is None:
        continue

    try:
        region = Region(
            geometry,
            properties={
                PropertyType.MATERIAL: "CHANNEL_BOX",
                PropertyType.MACRO: macro_name
            }
        )
        subregion_center = get_point_coordinates(region.o)
        assembly.add(region, position=subregion_center)
    except Exception as e:
        print(f"  Warning: Failed to create corner CB region {corner_id}: {e}")


elements = build_control_cross_elements(
    ap=assembly_pitch,
    assembly_center=center
)
for element_name in elements.keys():
    geom_obj = elements[element_name]["geometry"]
    material = elements[element_name]["material"]
    macro = elements[element_name]["macro"]
    
    if "offset" in elements[element_name].keys():
        offset = elements[element_name]["offset"] 
    else:
        offset = (0.0, 0.0, 0.0)

    try: 
        region = Region(
            geom_obj=geom_obj,
            properties={
                PropertyType.MATERIAL: material,
                PropertyType.MACRO: macro
            }
        )
        region_center = get_point_coordinates(region.o)
        offset_center = (region_center[0]+offset[0],region_center[1]+offset[1], region_center[2]+offset[2])
        assembly.add(region, position=offset_center)
    except Exception as e:
        print(f"  Warning failed to create Region for control cross element : {element_name}")

assembly.update_hierarchical_structure(True)

if apply_diagonal_symmetry:
    assembly.apply_symmetry(SymmetryType.DIAG)
    symmetry_id = "DIAG"
    symetry_type = SymmetryType.DIAG

else:
    assembly.apply_symmetry(SymmetryType.FULL)
    symmetry_id = "FULL"
    symetry_type = SymmetryType.FULL

assembly.show(GeometryType.TECHNOLOGICAL, PropertyType.MACRO)

## Export the lower diagonal applying TISO conditions and exporting both MATERIAL and MACRO properties

export_layout_to_tdt(
    assembly, f"{tdt_output_path}/PB2_Type6-C_{symmetry_id}_IC_MACRO_TISO_glow_test", 
            TdtSetup(GeometryType.TECHNOLOGICAL,
                    property_types=[PropertyType.MACRO, PropertyType.MATERIAL],
                    type_geo=SYM_TO_LAYOUT_AND_BOUNDARY[symetry_type]["TISO"]["layout"],
                    symmetry_type=symetry_type))    