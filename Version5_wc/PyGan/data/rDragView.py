#!/usr/bin/env python

# Purpose:
#  Launcher file.
# 
# Copyright:
#  Copyright (C) 2025 Polytechnique Montreal
#  This library is free software; you can redistribute it and/or modify it under the terms 
#  of the GNU Lesser General Public License as published by the Free Software Foundation; 
#  either version 2.1 of the License, or (at your option) any later version
# 
# Author(s): Atyab A. Calloo

import os
import sys

from utils import parserfunc
from get_lcm import read_lcm_obj, run_pydragon
from vtu_tools import create_vtu

verbose=1
test_name="3DHEX"

Geom, Mcro, Trck, Flux = run_pydragon(test_name,verbose)

### CREATE MATERIAL AND FLUX VTU FILES
create_vtu(Geom, Mcro, Trck, Flux, test_name, verbose)

print("test rDragView completed")
