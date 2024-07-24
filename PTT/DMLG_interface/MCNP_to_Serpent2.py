# Python 3 script to allow for the conversion of an MCNP input file to a Serpent2 input file
# Makes use of the MCNP input file parser from the DMLG interface tool.
# Author : R. Guasch
# Date : 2024-05-02

import numpy as np
import DMLGInterface as DMLG


materials_dict = {"COOL": 2, "MODE":3, "CHAN":5, "BOX":6, "CONTROL_CROSS":8, "CONTROL_CROSS_CORNER":9, "VOID":22, 
                  "UOX1":23, "UOX2":24, "UOX3":25, "UOX4":26, "UOX5":27, "UOX6":28, "GADO1":29, "GADO2":30}
read_MCNP_case = DMLG.DMLG_Interface("MCNP_AT10_sanitized.inp", code="MCNP",mode="input", reactor_type= "BWR")
#read_MCNP_case.getMCNP_card_data(print_cells=True, print_surfaces=True, print_materials=True)
read_MCNP_case.choose_output("Serpent2", materials_dict, printlevel = 1)