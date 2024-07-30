# Python3 script used to check analytical vs serpent2 vs dragon5 volumes for the ATRIUM-10 assembly.
# Purpose : vaidate geometry definition and volume integration.
# Author : R. Guasch
# Date : 2024-06-19

import DMLGInterface as DMLG
import GeometricTools as GeoT
import numpy as np

# Exection parameters : 

check_AT10_assbly = True
check_AT10_assbly_ctrl = False


def check_SerpentvsDragon_vols(Serpent2_case, Dragon5_case, material_assocoation_dict={}):
    """
    Serpent2_case : DMLG object corresponding to Serpent2 case
    Dragon5_case : DMLG object corresponding to equivalent Dragon5 case
    material_association_dict : dictionary associating Serpent2 materials to region numbers from Dragon5. 
    """
    if Serpent2_case.mode == "output":
        """
        In this case, serpent2 output file has been read, there might be issues with the volume of water since it is not properly bounded.
        """
        D5_regions_volumes = Dragon5_case.Dragon5_geom.getVolumesAndRegions()
        S2_regions_volumes = Serpent2_case.S2_mat_properties.getMaterialsandVolumes()
        print(D5_regions_volumes)
        print(S2_regions_volumes)
        errors=[]
        for mat in material_assocoation_dict.keys():
            region=material_assocoation_dict[mat]
            if mat == "H2O":
                errors.append(D5_regions_volumes[f"region {region}"]-8.47125E-01)
            else:
                errors.append(D5_regions_volumes[f"region {region}"]-S2_regions_volumes[mat])

        sum=0
        for err in errors:
            sum+=err**2
        rms=np.sqrt(sum/len(errors))
        print(f"errors on volumes are {errors}, with an RMS error of {rms}")
    elif Serpent2_case.mode == "check_volumes":
        """
        In this case, serpent2 volume evaluation has been performed using the -checkvolumes command. 
        All volumes are finite, however this does not allow for the creation of the "Materials properties" structure.
        Checking volumes with Dragon5 case : 
        """
        D5_regions_volumes = Dragon5_case.Dragon5_geom.getOrderedVolumesandRegions()
        print(f"D5 vols are: {D5_regions_volumes}")
        S2_material_volumes = Serpent2_case.S2_geom.getOrderedMaterialVols()
        print(f"S2 vols are: {S2_material_volumes}")
        # Check consistency between lenghts 
        if len(D5_regions_volumes) == len(S2_material_volumes):
            print(f"Consistent material volumes dict lenghts : S2 and D5 case both have {len(S2_material_volumes)} material volumes/Merged Mixes")
            errors=[]
            rel_errors =[]
            S2_vols = list(S2_material_volumes.values())
            D5_vols = list(D5_regions_volumes.values())
            for i in range(len(S2_vols)):
                error = D5_vols[i]-S2_vols[i]
                relative_error = (D5_vols[i]*2-S2_vols[i])*100/S2_vols[i]
                errors.append(error)
                rel_errors.append(relative_error)
                if abs(relative_error)>=0.01:
                    print(f"Warning in region {list(S2_material_volumes.keys())[list(S2_material_volumes.values()).index(S2_vols[i])]}, with S2 volume {S2_vols[i]/2} and D5 volume {D5_vols[i]*2}")
            sum=0
            for err in errors:
                sum+=err**2
            rms=np.sqrt(sum/len(errors))
            sum=0
            for rel_err in rel_errors:
                sum+=rel_err**2
            rel_rms = np.sqrt(sum/len(rel_errors))
            print(f"Absolute errors on volumes are {errors}, \n with an RMS error of {rms}")
            print(f"Relative errors on volumes are {rel_errors} (%), \n with an RMS error of {rel_rms} (%)")
        else:
            print("Unconsistent definition of S2 material volumes vs D5 regions, cannot compare volumes")





    return


if check_AT10_assbly:
    # Check un-controlled assembly volumes :

    assbly_serp_vols = "../Serpent2/AT10_ASSBLY_mc.mvol"
    assbly_drag_vols = "../Dragon5/SALT_ASSBLY_vols.txt"

    print("$$$ ---- Checking Dragon5 vs Serpent2 volumes for controlled assembly ---- $$$")
    AT10_ASSBLY_Serp = DMLG.DMLG_Interface(assbly_serp_vols, code="Serpent2", mode="check_volumes")
    AT10_ASSBLY_Serp.createS2_geom("ATRIUM-10 NL-24 UOX bundle", height=2)


    AT10_ASSBLY_drag = DMLG.DMLG_Interface(assbly_drag_vols, code="Dragon", mode="output")
    AT10_ASSBLY_drag.Dragon5_geom.ComputeOrderedVolumesandRegions()
    check_SerpentvsDragon_vols(AT10_ASSBLY_Serp, AT10_ASSBLY_drag)


    print("$$$ ---- Checking Analytical vs Serpent2 volumes for un-controlled assembly ---- $$$")
    # Checking with geometric data : Analytical evaluation of volumes by geom_ASSBLY class :
    Channel_box_out = 2.3975+1.1025
    Channel_box_in = 2.3175+1.0225
    pitch_A=7.62*2
    Box_o = 6.87*2
    Box_i = 6.7*2
    Channel_box_xL_out = -1.1025
    Channel_box_XR_out = 2.3975

    pins_names=["24UOx", "32UOx", "42UOx", "45UOx", "48UOx", "50UOx", "45GADO", "42GADO"]
    pin_radii =[0.4435,0.4520,0.5140] # Fuel, gap, clad radii
    AT10_volume_check = GeoT.geom_ASSBLY(pitch_A,pins_names, Box_o, Channel_box_out, Box_i, Channel_box_in)
    AT10_volume_check.setPins(pin_radii[0], pin_radii[1], pin_radii[2])

    pins_number_dict={"24UOx": 4, "32UOx" : 8, "42UOx": 10, "45UOx": 21, "48UOx": 6, "50UOx":26, "45GADO": 14, "42GADO": 2}
    AT10_volume_check.setNumberofPins(pins_number_dict)
    AT10_volume_check.computeVolumes()
    Analytical_volumes = AT10_volume_check.Volumes

    D5_regions_volumes = AT10_ASSBLY_drag.Dragon5_geom.getVolumesAndRegions()
    S2_regions_volumes = AT10_ASSBLY_Serp.S2_geom.getOrderedMaterialVols()
    print(f"Analytical volumes are {Analytical_volumes}")
    print(f"Dragon5 volumes are {D5_regions_volumes}")
    print(f"Serpent2 volumes are {S2_regions_volumes}")
    for region_name in Analytical_volumes.keys():
        if region_name in S2_regions_volumes.keys():
            S2_error = (S2_regions_volumes[region_name] - Analytical_volumes[region_name])*100/Analytical_volumes[region_name]
            print(f"Statistical error on volume of region {region_name} from Serpent2 is {S2_error:.03f} %")
        

if check_AT10_assbly_ctrl:
    # Check controlled assembly volumes :
    print("$$$ ---- Checking Analytical vs Serpent2 volumes for controlled assembly ---- $$$")
    assbly_ctrl_serp_vols = "../Serpent2/AT10_ASSBLY_CTRL_mc_pitchissue.mvol"

    Channel_box_out = 2.3975+1.1025
    Channel_box_in = 2.3175+1.0225
    pitch_A=7.62*2
    Box_o = 6.87*2
    Box_i = 6.7*2
    Channel_box_xL_out = -1.1025
    Channel_box_XR_out = 2.3975


    cross_half_length = 4.79654 + 7.62
    cross_half_width = 7.62 - 7.26784
    nb_ctrl_rods = 18*2
    ctrl_rod_radius = 0.21082
    ctrl_wings_symmetry = True

    # Analytical volumes for controlled assembly :
    AT10_CTRL_volume_check = GeoT.geom_ASSBLY(pitch_A,pins_names, Box_o, Channel_box_out, Box_i, Channel_box_in, controlled=True)
    AT10_CTRL_volume_check.setPins(pin_radii[0], pin_radii[1], pin_radii[2])
    AT10_CTRL_volume_check.setNumberofPins(pins_number_dict)
    AT10_CTRL_volume_check.setCtrlCross(cross_half_length, cross_half_width, nb_ctrl_rods, ctrl_rod_radius, ctrl_wings_symmetry)
    AT10_CTRL_volume_check.computeVolumes()
    Analytical_volumes_CTRL = AT10_CTRL_volume_check.Volumes


    AT10_ASSBLY_CTRL_Serp = DMLG.DMLG_Interface(assbly_ctrl_serp_vols, code="Serpent2", mode="check_volumes")
    AT10_ASSBLY_CTRL_Serp.createS2_geom("ATRIUM-10 NL-24 Controlled UOX bundle", height=1)
    S2_regions_volumes_CTRL = AT10_ASSBLY_CTRL_Serp.S2_geom.getOrderedMaterialVols()

    for region_name in Analytical_volumes_CTRL.keys():
        if region_name in S2_regions_volumes_CTRL.keys():
            print(f"$ -- Seperent2 volume for region {region_name} is {S2_regions_volumes_CTRL[region_name]}")
            print(f"$ -- Analytical volume for region {region_name} is {Analytical_volumes_CTRL[region_name]}")
            S2_error = (S2_regions_volumes_CTRL[region_name] - Analytical_volumes_CTRL[region_name])*100/Analytical_volumes_CTRL[region_name]
            print(f"Statistical error on volume of region {region_name} from Serpent2 is {S2_error:.03f} %")
