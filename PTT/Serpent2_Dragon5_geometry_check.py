import DMLGInterface as DMLG
import numpy as np
"""
input_f = "MCNP_AT10_sanitized.inp"
AT10_CRTL_MCNP = DMLG_Interface(input_f, type="MCNP",mode="input")
#print(AT10_CRTL_MCNP.Cell_Cards)
AT10_CRTL_MCNP.getMCNP_card_data(print_cells=False, print_surfaces=False, print_materials=True)
#print(AT10_CRTL_MCNP.Cell_Cards)
#print(AT10_CRTL_MCNP.Cell_Cards['water box centered at ( 8.267, 6.973)'].material_densities[2])

cell_serpent_output = "Serpent2/data_AT10_24UOX_try/AT10_24_test_mc.out"
AT10_24UOX_cell_serp = DMLG.DMLG_Interface(cell_serpent_output, type="Serpent2",mode="output")
cell_dragon_output = "../Version5_ev3232/Dragon/Linux_x86_64/SALT_TSPC.result"
AT10_24UOX_cell_drag = DMLG.DMLG_Interface(cell_dragon_output, type="Dragon",mode="output")


material_assocoation_dict = {"UOx_A":1, "UOx_B":2, "UOx_C":3, "UOx_D":4, "gap":5, "clad":6, "H2O":7}
#print(AT10_24UOX_cell_serp.Serpent2_cards)
AT10_24UOX_cell_serp.createS2_mat_properties("AT10_24UOX_cell",height=3.8E2)
print(f"Total fuel mass is : {AT10_24UOX_cell_serp.S2_mat_properties.getTotalFuelMass():.3f} g")
print(f"Total fuel mass dens is : {AT10_24UOX_cell_serp.S2_mat_properties.getTotalFuelMassDens():.3f} g/cm3")
"""

def check_SerpentvsDragon_vols(Serpent2_case, Dragon5_case, material_assocoation_dict):
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
                errors.append(D5_vols[i]-S2_vols[i])
                rel_errors.append((D5_vols[i]-S2_vols[i])*100/S2_vols[i])
            sum=0
            for err in errors:
                sum+=err**2
            rms=np.sqrt(sum/len(errors))
            sum=0
            for rel_err in rel_errors:
                sum+=rel_err**2
            rel_rms = np.sqrt(sum/len(rel_errors))
            print(f"Absolute errors on volumes are {errors}, with an RMS error of {rms}")
            print(f"Relative errors on volumes are {rel_errors} (%), with an RMS error of {rel_rms} (%)")





    return

#check_SerpentvsDragon_vols(AT10_24UOX_cell_serp, AT10_24UOX_cell_drag, material_assocoation_dict)

assbly_serp_vols = "/home/loutre/RESEARCH/PolyBWR_project/PTT/Serpent2/AT10_ASSBLY_mc.mvol"
AT10_ASSBLY_Serp = DMLG.DMLG_Interface(assbly_serp_vols, type="Serpent2",mode="check_volumes")
AT10_ASSBLY_Serp.createS2_geom("ATRIUM-10 bundle", 1)

assbly_drag_vols = "/home/loutre/RESEARCH/PolyBWR_project/PTT/Dragon5/SALT_volumes.txt"
AT10_ASSBLY_drag = DMLG.DMLG_Interface(assbly_drag_vols, type="Dragon",mode="output")
AT10_ASSBLY_drag.Dragon5_geom.ComputeOrderedVolumesandRegions()
check_SerpentvsDragon_vols(AT10_ASSBLY_Serp, AT10_ASSBLY_drag, material_assocoation_dict={})


