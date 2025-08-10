import DMLGInterface as DMLG
import GeometricTools as GeoT
import numpy as np
import os


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
        # Check consistency between lengths 
        print(f"Length of D5 volumes = {len(D5_regions_volumes)}")
        print(f"Length of S2 volumes = {len(S2_material_volumes)}")
        if len(D5_regions_volumes) == len(S2_material_volumes):
            print(f"Consistent material volumes dict lenghts : S2 and D5 case both have {len(S2_material_volumes)} material volumes/Merged Mixes")
            errors=[]
            rel_errors =[]
            S2_vols = list(S2_material_volumes.values())
            D5_vols = list(D5_regions_volumes.values())
            for i in range(len(S2_vols)):
                error = D5_vols[i]*8-S2_vols[i]/2
                relative_error = (D5_vols[i]*8-S2_vols[i]/2)*100/S2_vols[i]/2
                errors.append(error)
                rel_errors.append(relative_error)
                if abs(relative_error)>=0.01:
                    print(f"Warning in region {list(S2_material_volumes.keys())[list(S2_material_volumes.values()).index(S2_vols[i])]}, with S2 volume {S2_vols[i]/2} and D5 volume {D5_vols[i]*8}")
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


assbly_serp_vols = f"{os.environ['SERPENT_RESULTS']}/../data/OECD_NEA_PHASE_IIIB_inp.mvol"
assbly_drag_vols_SSH = "case_data/OECD_NEA_IIIB_SALT_volumes_SSH.txt"
assbly_drag_vols_FLX = "case_data/OECD_NEA_IIIB_SALT_volumes_FLX.txt"

ASSBLY_Serp = DMLG.DMLG_Interface(assbly_serp_vols, type="Serpent2",mode="check_volumes")
ASSBLY_Serp.createS2_geom("OECD NEA Phase IIIB benchmark", 1)
#print(AT10_ASSBLY_Serp.S2_geom.getOrderedMaterialVols()["box"])


ASSBLY_SSH_drag = DMLG.DMLG_Interface(assbly_drag_vols_SSH, type="Dragon",mode="output")
ASSBLY_SSH_drag.Dragon5_geom.ComputeOrderedVolumesandRegions()
check_SerpentvsDragon_vols(ASSBLY_Serp, ASSBLY_SSH_drag)

assbly_drag_vols_FLX = DMLG.DMLG_Interface(assbly_drag_vols_FLX, type="Dragon",mode="output")
assbly_drag_vols_FLX.Dragon5_geom.ComputeOrderedVolumesandRegions()
check_SerpentvsDragon_vols(ASSBLY_Serp, assbly_drag_vols_FLX)


# Checking with geometric data : Analytical evaluation of volumes by geom_ASSBLY class :
moder_box_out = 1.6
moder_box_in = 1.5
pitch_A=7.62*2
Box_o = 6.774*2
Box_i = 6.52*2
pins_names=["23UOX_1", "30UOX_2", "36UOX_3", "36UOX_4", "36UOX_5", "45Gd_6", "49UOX_7", "49UOX_8", "49UOX_9"] # 

pin_radii =[0.529,0.529,0.615] # Fuel, (gap radius = fuel radius : no gap in OECD benchmark), clad radii
Phase_IIIB_volume_check = GeoT.geom_ASSBLY(pitch_A,pins_names, Box_o, moder_box_out, Box_i, moder_box_in, moderating_box_type="cylinder")
Phase_IIIB_volume_check.setPins(pin_radii[0], pin_radii[1], pin_radii[2])

pins_number_dict={"23UOX_1": 4, "30UOX_2": 8, "36UOX_3": 8, "36UOX_4": 8, "36UOX_5": 4, "45Gd_6": 8, "49UOX_7": 8, "49UOX_8": 4, "49UOX_9": 8}
Phase_IIIB_volume_check.setNumberofPins(pins_number_dict)
Phase_IIIB_volume_check.computeVolumes()
Analytical_volumes = Phase_IIIB_volume_check.Volumes
print(Analytical_volumes)

D5_regions_ssh_volumes = ASSBLY_SSH_drag.Dragon5_geom.getVolumesAndRegions()
S2_regions_volumes = ASSBLY_Serp.S2_geom.getOrderedMaterialVols()

for region_name in Analytical_volumes.keys():
    print(f"Analyzing region volume : {region_name}")
    if "UOX" in region_name or "Gd" in region_name:
        #print(f"In region UOX or Gd, region_name = {region_name}")
        #print(f"analytical volume region name = {region_name}")
        #print(f"region_name.split() = {region_name.split('_')}")
        S2_region_name = f"{region_name.split('_')[0]}_{region_name.split('_')[2]}_{region_name.split('_')[1]}"
    elif region_name == "gap":
        region_name = "gap_1"
        S2_region_name = "gap_1"
        continue
    else: 
        print(f"region_name = {region_name}")
        S2_region_name = region_name
        print(f"Monte Carlo volume = {S2_regions_volumes[S2_region_name]/2}")
        print(f"Analytical volume = {Analytical_volumes[region_name]}")
    S2_error = (S2_regions_volumes[S2_region_name]/2 - Analytical_volumes[region_name])*100/Analytical_volumes[region_name]
    print(f"Statistical error on volume of region {region_name} from Serpent2 is {S2_error:.03f} %")
    