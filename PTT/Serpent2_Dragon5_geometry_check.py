import DMLGInterface as DMLG
"""
input_f = "MCNP_AT10_sanitized.inp"
AT10_CRTL_MCNP = DMLG_Interface(input_f, type="MCNP",mode="input")
#print(AT10_CRTL_MCNP.Cell_Cards)
AT10_CRTL_MCNP.getMCNP_card_data(print_cells=False, print_surfaces=False, print_materials=True)
#print(AT10_CRTL_MCNP.Cell_Cards)
#print(AT10_CRTL_MCNP.Cell_Cards['water box centered at ( 8.267, 6.973)'].material_densities[2])
"""
cell_serpent_output = "Serpent2/data_AT10_24UOX_try/AT10_24_test_mc.out"
AT10_24UOX_cell_serp = DMLG.DMLG_Interface(cell_serpent_output, type="Serpent2",mode="output")
cell_dragon_output = "../Version5_ev3232/Dragon/Linux_x86_64/SALT_TSPC.result"
AT10_24UOX_cell_drag = DMLG.DMLG_Interface(cell_dragon_output, type="Dragon",mode="output")

association_dict = {"UOx_A":1, "UOx_B":2, "UOx_C":3, "UOx_D":4, "gap":5, "clad":6, }
print(AT10_24UOX_cell_serp.Serpent2_cards)
