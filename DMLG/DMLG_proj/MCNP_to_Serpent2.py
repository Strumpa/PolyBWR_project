# Python 3 script to allow for the conversion of an MCNP input file to a Serpent2 input file
# Makes use of the MCNP input file parser from the DMLG interface tool.
# Author : R. Guasch
# Date : 2024-05-02

import numpy as np
import DMLGInterface as DMLG

read_MCNP_case = DMLG.DMLG_Interface("MCNP_AT10_sanitized.inp", type="MCNP",mode="input")
#read_MCNP_case.getMCNP_card_data(print_cells=True, print_surfaces=True, print_materials=True)
read_MCNP_case.choose_output("Serpent2")