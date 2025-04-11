#!/usr/bin/env python

# Purpose:
#  Utility functions
# 
# Copyright:
#  Copyright (C) 2025 Polytechnique Montreal
#  This library is free software; you can redistribute it and/or modify it under the terms 
#  of the GNU Lesser General Public License as published by the Free Software Foundation; 
#  either version 2.1 of the License, or (at your option) any later version
# 
# Author(s): Atyab A. Calloo

import argparse 


def parserfunc():
   """
   Parse the input on the command line to pick up relevant parameters

   Args:
      None
   Returns:
      dict: of arguments
   Raises:
      None
   Notes:
      None
   """
   parser = argparse.ArgumentParser(prog='rDragView', description='Create vtu objects for '
      'visualisation of DRAGON5 flux with Paraview. Currently works for SN 2D/D in Cartesian and '
      'hexagonal geometries. Default run assumes that the required LCM objects are in the '
      'directory from where this is run. ')
   parser.add_argument('-c', '--calc_type', required=False,
                       type=str, action="store", nargs=1,
                       default=["view"], choices={"view", "dragon_view"},
                       help='Specify type of calculation. "view": (default) use already available '
                       'LCM objects to compute vtu files. "dragon_view": run dragon followed by '
                       'vtu creation.')
   
   parser.add_argument('-n', '--test_name', required=False, metavar='CALC_NAME',
                       type=str, action="store", nargs=1,
                       help='Specify a test name if desired. Used for naming vtu files. ')
   
   parser.add_argument('-v', '--verbose', required=False, metavar='VERB',
                       type=int, action="store", nargs=1,
                       default=[1],
                       help='Verbose level, ie, level of print output.')

   args = parser.parse_args()

   ### ENSURE THAT TEST NAME GIVEN FOR dragon_view CALCULATION
   if args.calc_type[0]=='dragon_view' and (args.test_name is None):
      parser.error("dragon_view CALCULATION REQUIRES A TEST NAME.")

   return args
