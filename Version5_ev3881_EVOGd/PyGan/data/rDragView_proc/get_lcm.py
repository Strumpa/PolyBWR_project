#!/usr/bin/env python

# Purpose:
#  Read or create lcm objects.
# 
# Copyright:
#  Copyright (C) 2025 Polytechnique Montreal
#  This library is free software; you can redistribute it and/or modify it under the terms 
#  of the GNU Lesser General Public License as published by the Free Software Foundation; 
#  either version 2.1 of the License, or (at your option) any later version
# 
# Author(s): Atyab A. Calloo

from pathlib import Path

import os
import lifo
import lcm
import cle2000


def read_lcm_obj(verbose):
   """
   Reads Geom, Macr, Trck, Flux LCM objects from current working directory. 
   Expects file name to begin with an underscore and be capitalised. 

   Args:
      verbose         : controls write output level
   Returns:
      Geom            : geometry LCM object
      Macr            : macrolib LCM object
      Trck            : tracking LCM object
      Flux            : flux LCM object
   """
   ### FETCH VERBOSE LEVEL
   lcm_verbose = verbose-2 if verbose > 2 else 0

   ### PICK UP LCM OBJECTS
   Geom=lcm.new(pytype='LCM_INP',name='GEOM',iact=0,impx=lcm_verbose)
   Macr=lcm.new(pytype='LCM_INP',name='MACR',iact=0,impx=lcm_verbose)
   Trck=lcm.new(pytype='LCM_INP',name='TRCK',iact=0,impx=lcm_verbose)
   Flux=lcm.new(pytype='LCM_INP',name='FLUX',iact=0,impx=lcm_verbose)

   if verbose > 0:
      print('>>>>>>> SUCCESSFULLY LOADED LCM OBJECTS.')
      print('-'*74 + '\n')

   return Geom, Macr, Trck, Flux


def run_pydragon(test_name,verbose):
   """
   Executes c2m procedure using DRAGON/DONJON through the cle200 class in PyGan
   and returns Geom, Macr, Trck, Flux LCM objects.

   Args:
      test_name       : file name of the c2m procedure to be computed
      verbose         : controls write output level
   Returns:
      Geom            : geometry LCM object
      Macr            : macrolib LCM object
      Trck            : tracking LCM object
      Flux            : flux LCM object
   Raises:
      FileNotFound    : if provided file name (test_name) does not match any of the
                        files in current working directory.
   """

   ### CHECK IF FILE EXISTS
   extnsn = '.c2m'
   template_file = test_name+extnsn
   path_file = Path(template_file)
   if not path_file.exists():
        for entry in os.scandir('.'):
            if entry.is_file():
                print(entry.name)
        raise FileNotFoundError('FILE '+template_file+' NOT FOUND IN CURRENT DIRECTORY.')

   ### CONSTRUCT LIFO STACK
   ipLifo1=lifo.new()
   ipLifo1.pushEmpty("Geom", "LCM")
   ipLifo1.pushEmpty("Macr", "LCM")
   ipLifo1.pushEmpty("Trck", "LCM")
   ipLifo1.pushEmpty("Syst", "LCM")
   ipLifo1.pushEmpty("Flux", "LCM")

   ### RUN CLE-2000 PROCEDURE
   CLE2KFILE = cle2000.new(test_name,ipLifo1,1)
   CLE2KFILE.exec()

   ### RECOVER OUTPUT LCM OBJECTS
   Geom = ipLifo1.node("Geom")
   Macr = ipLifo1.node("Macr")
   Trck = ipLifo1.node("Trck")
   Flux = ipLifo1.node("Flux")

   ### EMPTY LIFO STACK
   while ipLifo1.getMax() > 0:
      ipLifo1.pop()

   if verbose > 0:
      print('>>>>>>> SUCCESSFULLY RAN ',template_file,' DRAGON/DONJON CALCULATION.')
      print('-'*74 + '\n')
      
   return Geom, Macr, Trck, Flux