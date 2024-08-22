# Python3 script used to check analytical vs serpent2 vs dragon5 volumes for 2x2 UOX sub-assembly ATRIUM-10 geometry.
# Purpose : vaidate geometry definition and volume integration.
# Author : R. Guasch
# Date : 2024-06-19

import DMLGInterface as DMLG
import GeometricTools as GeoT
import numpy as np


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
                relative_error = (D5_vols[i]-S2_vols[i])*100/S2_vols[i]
                errors.append(error)
                rel_errors.append(relative_error)
                if abs(relative_error)>=0.01:
                    print(f"Warning in region {list(S2_material_volumes.keys())[list(S2_material_volumes.values()).index(S2_vols[i])]}, with S2 volume {S2_vols[i]} and D5 volume {D5_vols[i]}")
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


assbly_serp_vols = "../Serpent2/AT10_2x2_UOX_mc.mvol"
assbly_drag_vols = "../Dragon5/SALT_2x2_UOX_vols.txt"

AT10_ASSBLY_Serp = DMLG.DMLG_Interface(assbly_serp_vols, code="Serpent2", mode="check_volumes")
AT10_ASSBLY_Serp.createS2_geom("ATRIUM-10 2x2 UOX bundle", height=1)


AT10_ASSBLY_drag = DMLG.DMLG_Interface(assbly_drag_vols, code="Dragon", mode="output")
AT10_ASSBLY_drag.Dragon5_geom.ComputeOrderedVolumesandRegions()
check_SerpentvsDragon_vols(AT10_ASSBLY_Serp, AT10_ASSBLY_drag)

