#!/usr/bin/env python

# Purpose:
#  Functions related to creating vtu files.
# 
# Copyright:
#  Copyright (C) 2025 Polytechnique Montreal
#  This library is free software; you can redistribute it and/or modify it under the terms 
#  of the GNU Lesser General Public License as published by the Free Software Foundation; 
#  either version 2.1 of the License, or (at your option) any later version
# 
# Author(s): Atyab A. Calloo

import numpy as np
import vtk
from math import sqrt


def create_vtu(Geom, Mcro, Trck, Flux, test_name, verbose):
   """
   Uses Geom, Macr, Trck, Flux LCM objects to create the material and flux vtu
   objects that can then be viewed in Paraview. 

   Args:
      Geom            : geometry LCM object
      Macr            : macrolib LCM object
      Trck            : tracking LCM object
      Flux            : flux LCM object
      test_name       : file name of the c2m procedure to be computed
      verbose         : controls write output level
   """

   ### RECOVER STATEVECTORS FROM LCM OBJECTS
   SV_geom = Geom["STATE-VECTOR"]
   SV_mcro = Mcro["STATE-VECTOR"]
   SV_trck = Trck["STATE-VECTOR"]
   SV_flux = Flux["STATE-VECTOR"]
   GV_mcro = Mcro["GROUP"]

   ### RECOVER RELEVANT MATERIAL, TRACKING AND FLUX DATA
   n_mat = SV_mcro[1]
   n_reg = SV_trck[0]
   n_unk = SV_trck[1]
   n_dim = SV_trck[8]
   spltl = SV_trck[25]
   n_grp = SV_flux[0]
   geom_type = SV_trck[5] # 5 - Car2D; 7 - Car3D; 8 - Hex2D; 9 - Hex3D
   ### RECOVER OTHER DATA FROM LCM OBJECTS
   mat_array  = Trck["MATCOD"]

   ### INITIALIZE FLUX AND NUSIGF ARRAYS
   flx_array = np.zeros((n_grp,n_unk), dtype='f')
   nsf_mat_array = np.zeros((n_grp,n_mat), dtype='f')
   nsf_reg_array = np.zeros((n_grp,n_reg), dtype='f')
   ### PICK UP FLUX AND NUSIGF ARRAYS 
   for i_g in range(n_grp):
      flx_array[i_g,:] = Flux["FLUX"][i_g]
      nsf_mat_array[i_g,:] = GV_mcro[i_g]["NUSIGF"]
      ### CREATE NUSIGF ARRAY BY REGION
      for i_reg in range(n_reg):
         ind_mat = mat_array[i_reg]
         if ind_mat == 0 : 
            nsf_reg_array[i_g,i_reg] = 0
         elif ind_mat > 0 : 
            ind_mat=ind_mat-1
            nsf_reg_array[i_g,i_reg] = nsf_mat_array[i_g,ind_mat]

   ### REBUILD MESHES
   lx, meshx, ly, meshy, lz, meshz = rebuild_meshes(Geom)

   ### FIND HEXAGON CENTRES IF NEEDED
   if geom_type==8 or geom_type==9:
      side = Geom["SIDE"][0]
      n_hex=int(n_reg/(lz*3*spltl**2))
      n_ring = int((sqrt(  float((4*n_hex-1)/3)  )+1.)/2.)
      hex_centres = compute_hex_centres(n_hex,n_ring,side)

   ### CALL TO MAKE_VTU GEOM-DEPENDENT
   nscat = SV_trck[6]
   ielem = SV_trck[7]
   if (geom_type == 5) or (geom_type == 7):
      if verbose > 0: 
         print('\n>>>>>>> GENERATING VTUs FOR ',test_name,'.')
         print('-'*74 + '\n')
      make_vtu_cartesian_geom(lx,meshx,ly,meshy,lz,meshz,mat_array,"material",test_name )
      for i_g in range(n_grp):
         make_vtu_cartesian_geom(lx,meshx,ly,meshy,lz,meshz,flx_array[i_g,:],"flux",test_name ,nscat=nscat,ielem=ielem,i_g=i_g)
   elif (geom_type == 8) or (geom_type == 9):
      if verbose > 0: 
         print('\n>>>>>>> GENERATING VTUs FOR ',test_name,'.')
         print('-'*74 + '\n')
      make_vtu_hexagonal_geom(n_hex,side,hex_centres,mat_array,"material",test_name ,spltl=spltl,meshz=meshz,lz=lz)
      for i_g in range(n_grp):
         make_vtu_hexagonal_geom(n_hex,side,hex_centres,flx_array[i_g,:],"flux",test_name ,meshz=meshz,lz=lz,nscat=nscat,ielem=ielem,spltl=spltl,i_g=i_g)

   print("\n>>>>>>> VTU CREATION COMPLETED!")
   print('-'*74 + '\n')


########################################################################
def rebuild_meshes(Geom):
   """
   Uses Geom LCM objects to return number of mesh elements and mesh arrays along each 
   axis, taking split into account.

   Args:
      Geom            : geometry LCM object
   Returns:
      lx              : number of mesh elements along x
      meshx           : mesh array along x
      ly              : number of mesh elements along y
      meshy           : mesh array along y
      lz              : number of mesh elements along z
      meshz           : mesh array along z
   """

   SV_geom = Geom["STATE-VECTOR"]
   geom_type = SV_geom[0] # 5 - Car2D; 7 - Car3D; 8 - Hex2D; 9 - Hex3D

   mesh_dict = {}

   ### FOR CARTESIAN GEOMETRIES
   if geom_type==5 or geom_type==7:
      lx    = SV_geom[2]
      meshx = Geom["MESHX"]
      spltx = Geom["SPLITX"]
      lx, meshx = get_mesh_oneaxis(lx,meshx,spltx)
      ly    = SV_geom[3]
      meshy = Geom["MESHY"]
      splty = Geom["SPLITY"]
      ly, meshy = get_mesh_oneaxis(ly,meshy,splty)
   else:
      lx   =1
      meshx=1
      ly   =1
      meshy=1

   ### FOR 3D GEOMETRIES
   if geom_type==7 or geom_type==9:
      lz    = SV_geom[4]
      meshz = Geom["MESHZ"]
      spltz = Geom["SPLITZ"]
      lz, meshz = get_mesh_oneaxis(lz,meshz,spltz)
   else:
      lz   =1
      meshz=1
      spltz=1

   return lx, meshx, ly, meshy, lz, meshz


########################################################################
def get_mesh_oneaxis(lx,mesh,splt):
   """
   Compute number of mesh elements and mesh arrays taking split into account.

   Args:
      lx              : number of mesh elements along axis w/o splits
      meshx           : mesh array along axis w/o splits
      splt            : split value along axis
   Returns:
      lx              : number of mesh elements along axis w/ splits
      meshx           : mesh array along axis w/ splits
   """

   tmp_mesh=np.empty(1, dtype='f') 
   val = 0.
   tmp_mesh[0] = val
   for i in range(len(splt)):
      step=(mesh[i+1]-mesh[i])/splt[i]
      for j in range(splt[i]):
         val += step
         tmp_mesh = np.append(tmp_mesh,val)

   lx=len(tmp_mesh)-1

   return lx, tmp_mesh


########################################################################
def compute_hex_centres(n_hex,n_ring,side):
   """
   Compute centres of hexagons within hexagonal domain. Centre hexagon is zeroth
   ring.
   Args:
      n_hex           : number of hexagons in domain
      n_ring          : number of rings in domain; centre hexagon is zeroth ring
      side            : length of one side of hexagon
   Returns:
      hex_centres     : 2D array with x and y coordinates of hexagon centres
   """

   hex_centres = np.zeros((2,n_hex), dtype='f')
   ind_hex =0
   hex_centres[0, 0] = 0.0
   hex_centres[1, 0] = 0.0

   ### LOOP OVER EACH RING OF HEXAGONS
   for ind_ring in range(1, n_ring):
      numhex_ring = ind_ring * 6
      xcoord = ind_ring * 1.5 * side
      ycoord = ind_ring * (sqrt(3) / 2) * side

      hex_centres[0, ind_hex + 1] = xcoord
      hex_centres[1, ind_hex + 1] = ycoord
      ind_hex += 1

      ind_dir = 1
      ind_step = 0
      numhex_perdir = ind_ring

      ### LOOP OVER HEXAGONS IN CURRENT RING
      ### For each of the six major axis directions, move in the direction
      ### of that axis for numhex_perdir hexagons (which happens to be the
      ### same as the ring index).
      for i_hex in range(1, numhex_ring):
         ind_step += 1
         if ind_dir == 1:
            xcoord -= 1.5 * side
            ycoord += (sqrt(3) / 2) * side
         elif ind_dir == 2:
            xcoord -= 1.5 * side
            ycoord -= (sqrt(3) / 2) * side
         elif ind_dir == 3:
            ycoord -= sqrt(3) * side
         elif ind_dir == 4:
            xcoord += 1.5 * side
            ycoord -= (sqrt(3) / 2) * side
         elif ind_dir == 5:
            xcoord += 1.5 * side
            ycoord += (sqrt(3) / 2) * side
         elif ind_dir == 6:
            ycoord += sqrt(3) * side

         hex_centres[0, ind_hex + 1] = xcoord
         hex_centres[1, ind_hex + 1] = ycoord
         ind_hex += 1

         if ind_step == numhex_perdir:
            ind_dir += 1
            ind_step = 0

   return hex_centres


########################################################################
def writeXML(ugrid, file_name):
   """
   Save unstructured grid vtk object (ugrid) to file.
   Args:
      ugrid           : vtkUnstructuredGrid object
      file_name       : file name under which to save ugrid
   Returns:
      None
   """
   writer = vtk.vtkXMLUnstructuredGridWriter()
   writer.SetFileName(file_name)
   writer.SetInputData(ugrid)
   writer.Write()


########################################################################
def make_vtu_cartesian_geom(
   lx,meshx,ly,meshy,lz,meshz,
   data_np,
   dataset_name,test_name,
   **kwargs
   ):
   """
   Build vtk unstructured grid object for 2D or 3D Cartesian mesh and populate with
   given data, material number or flux.
   Args:
      lx              : number of mesh elements along x
      meshx           : mesh array along x
      ly              : number of mesh elements along y
      meshy           : mesh array along y
      lz              : number of mesh elements along z
      meshz           : mesh array along z
      data_np         : given data, material number or flux
      dataset_name    : data name, usually 'material' or 'flux'
      test_name       : file name of the computed c2m procedure in case of a prior
                        dragon run or given name by user or default name
   Returns:
      None
   """

   ### ASSIGN OPTIONAL ARGS
   nscat = kwargs.get("nscat", 1)
   ielem = kwargs.get("ielem", 1)
   i_g = kwargs.get("i_g", 0)
   if dataset_name=="rate": 
      nsf_reg=kwargs.get("nusigf", 1)

   ### 2D CONSIDERATIONS
   if lz==1: meshz = [0,0]
   ndim = 2 if lz==1 else 3

   ### INITIALISE VTK OBJECTS
   points = vtk.vtkPoints()
   element = vtk.vtkVoxel()
   cells = vtk.vtkCellArray()
   ugrid = vtk.vtkUnstructuredGrid()
   ### INITIALISE SCALAR CELL DATA
   DRAGONdataset = vtk.vtkFloatArray()
   DRAGONdataset.SetName(dataset_name)
   ugrid.GetCellData().SetScalars(DRAGONdataset)

   ### COMPUTE CONSTANTS
   flux_offset = 1 if dataset_name=="material" else nscat*(ielem)**ndim

   i_cell=-1
   ### LOOP OVER Z-AXIS MESH
   for kk in range(lz):

      ### LOOP OVER Y-AXIS MESH
      for jj in range(ly):
         ycoord = [meshy[jj], meshy[jj], meshy[jj+1], meshy[jj+1]]

         ### LOOP OVER X-AXIS MESH
         for ii in range(lx):
            xcoord = [meshx[ii], meshx[ii+1], meshx[ii+1], meshx[ii]]

            ### INSERT NEW POINTS AND CELLS (DO NOT MERGE THE TWO LOOPS BELOW!!!)
            for i_crd in range(0,4):
               points.InsertNextPoint( xcoord[i_crd],ycoord[i_crd], meshz[kk])
            for i_crd in range(0,4):
               points.InsertNextPoint( xcoord[i_crd],ycoord[i_crd], meshz[kk+1])
            for i in range(0, 8):
               i_cell += 1
               element.GetPointIds().SetId(i, i_cell)
            cells.InsertNextCell(element)

            ugrid.SetPoints(points)
            ugrid.SetCells(vtk.VTK_HEXAHEDRON, cells)

            ### SET SCALAR CELL DATA
            ind1 = (kk*lx*ly) + (jj*lx) + ii
            ind2 = (kk*lx*ly)*flux_offset + (jj*lx)*flux_offset + ii*flux_offset
            if dataset_name=="rate":
               DRAGONdataset.InsertTuple1(ind1, nsf_reg[ind1]*data_np[ind2])
            else:
               DRAGONdataset.InsertTuple1(ind1, data_np[ind2])

   print(test_name, dataset_name,' PROCESSED. SAVING ...')
   filename=test_name+'_'+dataset_name+str(i_g)+".vtu"
   ### SAVE VTU FILE
   writeXML(ugrid, filename)
   return ugrid


########################################################################
def make_vtu_hexagonal_geom(
   n_hex,side,hex_centres,
   data_np,
   dataset_name,test_name,
   **kwargs
   ):
   """
   Build vtk unstructured grid object for 2D or 3D hexagon mesh and populate with
   given data, material number or flux.
   Args:
      n_hex           : number of hexagons in domain
      side            : length of one side of hexagon
      hex_centres     : 2D array with x and y coordinates of hexagon centres
      data_np         : given data, material number or flux
      dataset_name    : data name, usually 'material' or 'flux'
      test_name       : file name of the computed c2m procedure in case of a prior
                        dragon run or given name by user or default name
   Returns:
      None
   """
   ### ASSIGN OPTIONAL ARGS
   spltl = kwargs.get("spltl", 1)
   lz = kwargs.get("lz", 1)
   meshz = kwargs.get("meshz", 1)
   nscat = kwargs.get("nscat", 1)
   ielem = kwargs.get("ielem", 1)
   i_g = kwargs.get("i_g", 0)

   ### 2D CONSIDERATIONS
   if lz==1: meshz = [0,0]
   ndim = 2 if lz==1 else 3

   ### INITIALISE VTK OBJECTS
   points = vtk.vtkPoints()
   lozenge = vtk.vtkVoxel()
   cells = vtk.vtkCellArray()
   ugrid = vtk.vtkUnstructuredGrid()
   ### INITIALISE SCALAR CELL DATA
   DRAGONdataset = vtk.vtkFloatArray()
   DRAGONdataset.SetName(dataset_name)
   ugrid.GetCellData().SetScalars(DRAGONdataset)

   ### COMPUTE CONSTANTS
   x_step=(1.0)*(side/spltl)
   y_step=(sqrt(3)/2)*(side/spltl)
   n_loz = 3*spltl**2
   flux_offset = 1 if dataset_name=="material" else nscat*(ielem)**ndim

   i_cell = -1
   ### LOOP OVER Z-AXIS LAYERS
   for i_z in range(0,lz):
      ### LOOP OVER HEXAGONS
      for i_hex in range(n_hex):
         ### FETCH CENTRE OF EACH HEX
         cx, cy = hex_centres[0, i_hex], hex_centres[1, i_hex]

         ### LOOP OVER EACH HEX SECTION (EACH OF THE THREE MAIN LOZENGES)
         for hex_section in range(3):
            ### SET INITIAL VALUES OF XCOORD AND YCOORD
            if hex_section == 0:
               x_coords = [0.0, x_step*0.5, x_step, x_step*0.5]
               y_coords = [0.0, -y_step, 0.0, y_step]
            elif hex_section == 1:
               x_coords = [-side, -side+x_step, -side+(x_step*1.5), -side+(x_step * 0.5)]
               y_coords = [0.0, 0.0, y_step, y_step]
            else:
               x_coords = [-side/2, (-side/2)+x_step, (-side/2)+(x_step*0.5), (-side/2)-(x_step*0.5)]
               y_coords = [-side*sqrt(3)/2, -side*sqrt(3)/2, (-side*sqrt(3)/2)+y_step, (-side*sqrt(3)/2)+y_step]

            ### LOOP OVER Y-MESH INSIDE LOZENGE, INCREMENT XCOORD, YCOORD
            for i_mshy in range(spltl):
               if i_mshy != 0:
                  if hex_section == 0:
                     x_coords = x_coords - ((x_step/2)*(spltl-2))
                     y_coords = y_coords + y_step*(spltl)
                  elif hex_section == 1:
                     x_coords = x_coords - ((x_step)*(spltl-1.5))
                     y_coords = y_coords + y_step
                  else:
                     x_coords = x_coords - ((x_step)*(spltl-0.5))
                     y_coords = y_coords + y_step

               ### LOOP OVER X-MESH INSIDE LOZENGE, INCREMENT XCOORD, YCOORD
               for i_mshx in range(spltl):
                  if i_mshx != 0:
                     x_coords = x_coords + (x_step/2) if hex_section == 0 else x_coords + (x_step)
                     y_coords = y_coords - (y_step) if hex_section == 0 else y_coords

                  ### INSERT NEW POINTS AND CELLS (DO NOT MERGE THE TWO LOOPS BELOW!!!)
                  for i_crd in range(0,4):
                     points.InsertNextPoint(x_coords[i_crd]+cx, y_coords[i_crd]+cy, meshz[i_z])
                  for i_crd in range(0,4):
                     points.InsertNextPoint(x_coords[i_crd]+cx, y_coords[i_crd]+cy, meshz[i_z+1])
                  for i in range(0, 8):
                     i_cell += 1
                     lozenge.GetPointIds().SetId(i, i_cell)
                  cells.InsertNextCell(lozenge)

         ugrid.SetPoints(points)
         ugrid.SetCells(vtk.VTK_HEXAHEDRON, cells)

         ### SET SCALAR CELL DATA
         for i_loz in range(0,n_loz):
            ind1 = (i_z*n_hex*n_loz) + (i_hex*n_loz) + i_loz
            ind2 = (i_z*n_hex*n_loz*flux_offset) + (i_hex*n_loz*flux_offset) + (i_loz*flux_offset)
            DRAGONdataset.InsertTuple1(ind1, data_np[ind2])

   print(test_name, dataset_name,' PROCESSED. SAVING ...')
   filename=test_name+'_'+dataset_name+str(i_g)+".vtu"
   ### SAVE VTU FILE
   writeXML(ugrid, filename)
   return ugrid
