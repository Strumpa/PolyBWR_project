## PyGan script used to initialize tracking operators in DRAGON5, and recover the tracking files according to the methods selected.
# The available tracking operators are modules SALT:, NXT: and SYBILT:
# Hoever the last two are not compatible (yet) with 3-level BWR assemblies, so only SALT: is used here for now.
# Author : R. Guasch
# Date : 2025-06-02

import lcm
import cle2000
import lifo


def trackFluxGeomSALT(pyGEOM, num_angles, line_density, refl_type, anisotropy_level, solution_door, moc_angular_quandrature, nmu, batch, ps_file):
    """
    Function to track a main flux geometry object using the SALT: module in DRAGON5.
    Parameters : 
    ----------
    pyGEOM : (LCM object)
        Native geometry object to be tracked, generated by GEO: module, [flux geometry]. 
    num_angles : (int)
        Number of angles used for tracking.
    line_density : (float)
        Line density used for tracking.
    refl_type : (str)
        Type of reflection used at the tracking step, can be "TISO" for isotropic relfections or "TSPC" for specular reflections.
        Note : "TSPC" is recommended to properly model boundaries in irregular lattices.
    anisotropy_level : (int)
        Level of anisotropy for the tracking, can be 1 (isotropic), 2 (linearly anisotropic), 3 (anisotropy order P_2), or 4 (anisotropy order P_3).
    solution_door : (string)
        Flag to indicate whether the tracking should be modified for a MOC solution.
    moc_angular_quandrature : (str)
        Angular quadrature scheme used for MOC tracking, can be "LCMD", "OPP1", "OGAU", "GAUS", "DGAU", "CACA", or "CACB".
    num_polar_angles : (int)
        Number of polar angles used for MOC tracking.
    batch : (int)
        Number of batches to be used in the tracking process.
        To be optimized for parallel processing. --> depends on the number of tracks ie angles and line density, as well as number of regions.
    ps_file : (str)
        Name of the postscript file to be generated with the tracking results.

    Returns :
    ----------
    track_lcm : (lcm object)
        LCM object to store the tracking results in LCM format.
    track_binary : (binary)
        Sequential binary tracking file used to store the tracks lengths.
    figure : (ASCII)
        postscript file containing the tracking results in ASCII format, used for data visualization.
    
    """
    print("Tracking geometry with SALT: module...")
    print("Number of angles: ", num_angles)
    print("Line density: ", line_density)
    print("Reflection type: ", refl_type)
    print("Batch size: ", batch)
    print("Postscript file: ", ps_file)

    # Create a new LIFO object to manage the tracking process
    myLifo = lifo.new()

    # Push the geometry object onto the LIFO stack
    myLifo.pushEmpty('TRACK', "LCM")
    myLifo.pushEmpty('TF_EXC', "BINARY")
    myLifo.push(pyGEOM)
    myLifo.pushEmpty(ps_file, "ASCII")
    # Push the parameters for tracking onto the LIFO stack
    myLifo.push(num_angles)
    myLifo.push(line_density)
    myLifo.push(refl_type)
    myLifo.push(anisotropy_level)
    myLifo.push(solution_door)
    myLifo.push(moc_angular_quandrature)
    myLifo.push(nmu)
    myLifo.push(batch)

    # Create a new tracking operator using the SALT: module
    track_proc = cle2000.new('TRK_A_SALT', myLifo, 1)
    # Execute the tracking operation
    track_proc.exec()
    
    # Recover the results from the LIFO stack
    myLifo.lib()
    track_lcm = myLifo.node('TRACK')
    track_binary = myLifo.node('TF_EXC')
    figure = myLifo.node(ps_file)

    # empty the Lifo stack for TRK_A_SALT
    while myLifo.getMax() > 0:
        myLifo.pop();
    

    return track_lcm, track_binary, figure

def trackSSHGeomSALT(pyGEOMSSH, num_angles_ssh, line_density_ssh, refl_type, batch, ps_file_ssh):
    """
    Function to track a self-shielding geometry object using the SALT: module in DRAGON5.
    Parameters : 
    ----------
    pyGEOMSSH : (LCM object)
        Native geometry object to be tracked, generated by GEO: module, [self-shielding geometry]. 
    num_angles_ssh : (int)
        Number of angles used for self-shielding tracking.
    line_density_ssh : (float)
        Line density used for self-shielding tracking.
    refl_type : (str)
        Type of reflection used at the tracking step, can be "TISO" for isotropic relfections or "TSPC" for specular reflections.
        Note : "TSPC" is recommended to properly model boundaries in irregular lattices.
    batch : (int) 
        Number of batches to be used in the tracking process.
        To be optimized for parallel processing. --> depends on the number of tracks ie angles and line density, as well as number of regions.
    ps_file_ssh : (str)
        Name of the postscript file to be generated with the tracking results.


    Returns :
    ----------
    track_lcm : (lcm object)
        LCM object to store the tracking results in LCM format.
    track_binary : (binary)
        Sequential binary tracking file used to store the tracks lengths.
    figure : (ASCII)
        postscript file containing the tracking results in ASCII format, used for data visualization.
    
    """
    # Create a new LIFO object to manage the tracking process
    myLifo = lifo.new()

    # Push the geometry object onto the LIFO stack
    myLifo.pushEmpty('TRACK_SSH', "LCM")
    myLifo.pushEmpty('TF_EXC_SSH', "BINARY")
    myLifo.push(pyGEOMSSH)
    myLifo.pushEmpty(ps_file_ssh, "ASCII")
    # Push the parameters for tracking onto the LIFO stack
    myLifo.push(num_angles_ssh)
    myLifo.push(line_density_ssh)
    myLifo.push(refl_type)
    myLifo.push(batch)

    # Create a new tracking operator using the SALT: module
    track_proc = cle2000.new('TRK_SSH_SALT', myLifo, 1)

    # Execute the tracking operation
    track_proc.exec()
    myLifo.lib()
    # Recover the results from the LIFO stack
    ssh_track_lcm = myLifo.node('TRACK_SSH')
    ssh_track_binary = myLifo.node('TF_EXC_SSH')
    ssh_figure = myLifo.node(ps_file_ssh)

    # empty the Lifo stack for TRK_A_SALT
    while myLifo.getMax() > 0:
        myLifo.pop();
    

    return ssh_track_lcm, ssh_track_binary, ssh_figure



