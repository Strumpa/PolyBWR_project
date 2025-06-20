
/**********************************/
/* C API for hdf5 file support    */
/* (auxiliary functions)          */
/* author: A. Hebert (30/11/2021) */
/**********************************/

/*
   Copyright (C) 2021 Ecole Polytechnique de Montreal

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
 */

#if defined(HDF5_LIB)
#include <stdlib.h>
#include <string.h>
#include "hdf5_aux.h"
#define MAX_NAME 1024

hid_t h5tools_get_native_type(hid_t type);
herr_t print_info(hid_t ifile, const char *name, void *opdata);  /* Operator function */

herr_t iretcd;

static char AbortString[164];
static char name1024[1024];

void hdf5_open_file_c(const char *fname, hid_t *ifile, int_32 irdonly) {
/*
 *----------------------------------------------------------------------
 *
 * Open an HDF5 file
 *
 * input parameters:
 *   fname : hdf5 file namr
 *   ifile : hdf5 file identificator.
 *
 *----------------------------------------------------------------------
 */
  char *nomsub="hdf5_open_file_c";
  switch (irdonly) {
  case 0:
    *ifile = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
    if (*ifile < 0) {
      sprintf(AbortString,"%s: H5Fopen failure on read-write HDF5 file '%s'.",nomsub,fname);
      xabort_c(AbortString);
     }
    break;
  case 1:
    *ifile = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (*ifile < 0) {
      sprintf(AbortString,"%s: H5Fopen failure on read-only HDF5 file '%s'.",nomsub,fname);
      xabort_c(AbortString);
     }
    break;
  default:
    sprintf(AbortString,"%s: Invalid action on HDF5 file '%.72s'.",nomsub,fname);
    xabort_c(AbortString);
    break;
  }
}

void hdf5_close_file_c(hid_t *ifile) {
/*
 *----------------------------------------------------------------------
 *
 * Close an HDF5 file
 *
 * input parameters:
 *   ifile : hdf5 file identificator.
 *
 *----------------------------------------------------------------------
 */
  char *nomsub="hdf5_close_file_c";
  iretcd = H5Fclose(*ifile);
  if (iretcd != 0) {
    sprintf(AbortString,"%s: HDF5 close failure. iretcd=%d\n", nomsub, iretcd);
    xabort_c(AbortString);
  }
}

void hdf5_list_c(hid_t *ifile, const char *namp) {
/*
 *----------------------------------------------------------------------
 *
 * List the root table of contents
 *
 * input parameters:
 *   ifile : hdf5 file identificator.
 *
 *----------------------------------------------------------------------
 */
  char *nomsub="hdf5_list_c";
  H5G_stat_t statbuf;
  hid_t loc_id=*ifile;
  sprintf(name1024,"/%s",namp);
  printf("\nTable of contents --'%s'--\n", name1024);
  if(!hdf5_group_exists_c(ifile, namp)) {
    printf("%s: HDF5 missing group=%s\n", nomsub, namp);
    return;
  }
  H5Gget_objinfo(loc_id, name1024, 0, &statbuf);
  switch (statbuf.type) {
  case H5G_GROUP:
    printf("dataset/group name........................                                  ");
    printf("type............   bytes..    shape.................\n");
    H5Giterate(loc_id, name1024, NULL, print_info, NULL);
    break;
  default:
    printf("%s is not a H5G_GROUP.\n", namp);
    break;
  }
}

void hdf5_get_dimensions_c(hid_t *ifile, const char *namp, int_32 *rank) {
/*
 *----------------------------------------------------------------------
 *
 * Find dataset rank (number of dimensions)
 *
 * input parameters:
 *   ifile : hdf5 file identificator.
 *    namp : character name of the dataset.
 *
 * output parameter:
 *   rank  : number of dimensions.
 *
 *----------------------------------------------------------------------
 */
  char *nomsub="hdf5_get_dimensions_c";
  hid_t loc_id=*ifile;
  
  hid_t dataset = H5Dopen(loc_id,namp,H5P_DEFAULT);
  hid_t filespace = H5Dget_space (dataset);
  *rank = (int_32)H5Sget_simple_extent_ndims(filespace);
  if (*rank < 0) {
      sprintf(AbortString,"%s: H5Sget_simple_extent_ndims failure on object '%.72s'.",nomsub,namp);
      xabort_c(AbortString);
  }
  H5Dclose(dataset);
}

void hdf5_get_num_group_c(hid_t *ifile, const char *namp, int_32 *nbobj) {
/*
 *----------------------------------------------------------------------
 *
 * Find the number of objects (daughter datasets and daughter groups) in a group
 *
 * input parameters:
 *   ifile : hdf5 file identificator.
 *    namp : character name of the group.
 *
 * output parameter:
 *   nbobj : number of objects
 *
 *----------------------------------------------------------------------
 */
  char *nomsub="hdf5_get_num_group_c";
  hid_t loc_id=*ifile;

  hid_t group = H5Oopen(loc_id,namp,H5P_DEFAULT);
  if (group < 0) {
      sprintf(AbortString,"%s: H5Oopen failure on object '%.72s'.",nomsub,namp);
      xabort_c(AbortString);
  }
  hsize_t nobj;
  iretcd = H5Gget_num_objs(group, &nobj);
  if (iretcd < 0) {
      sprintf(AbortString,"%s: H5Gget_num_objs failure on object '%.72s'.",nomsub,namp);
      xabort_c(AbortString);
  }
  *nbobj = (int_32)nobj;
  H5Oclose(group);
}

void hdf5_list_datasets_c(hid_t *ifile, const char *namp, int_32 *ndsets, char *idata) {
/*
 *----------------------------------------------------------------------
 *
 * Recover character daughter dataset names in a group.
 *
 * input parameters:
 *   ifile : hdf5 file identificator.
 *    namp : character name of the group.
 *
 * output parameter:
 *  ndsets : number of daughter datasets.
 *   idata : dataset names.
 *
 *----------------------------------------------------------------------
 */
  char *nomsub="hdf5_list_datasets_c";
  hid_t loc_id=*ifile;
  int idx;

  hid_t group = H5Oopen(loc_id,namp,H5P_DEFAULT);
  if (group < 0) {
      sprintf(AbortString,"%s: H5Oopen failure on object '%.72s'.",nomsub,namp);
      xabort_c(AbortString);
  }
  hsize_t nobj;
  int otype;
  char memb_name[MAX_NAME];
  iretcd = H5Gget_num_objs(group, &nobj);
  if (iretcd < 0) {
      sprintf(AbortString,"%s: H5Gget_num_objs failure on object '%.72s'.",nomsub,namp);
      xabort_c(AbortString);
  }
  *ndsets = 0;
  for (idx = 0; idx < nobj; idx++) {
    H5Gget_objname_by_idx(group, (hsize_t)idx, memb_name, (size_t)MAX_NAME);
    otype =  H5Gget_objtype_by_idx(group, (size_t)idx);
    switch(otype) {
    case H5G_DATASET:
      strncpy (idata+MAX_NAME*(*ndsets), memb_name, MAX_NAME);
      (*ndsets)++;
      break;
    }
  }
  H5Oclose(group);
}

void hdf5_list_groups_c(hid_t *ifile, const char *namp, int_32 *ndsets, char *idata) {
/*
 *----------------------------------------------------------------------
 *
 * Recover character daughter group names in a group.
 *
 * input parameters:
 *   ifile : hdf5 file identificator.
 *    namp : character name of the group.
 *
 * output parameter:
 *  ndsets : number of daughter groups.
 *   idata : group names.
 *
 *----------------------------------------------------------------------
 */
  char *nomsub="hdf5_list_groups_c";
  hid_t loc_id=*ifile;
  int idx;

  hid_t group = H5Oopen(loc_id,namp,H5P_DEFAULT);
  if (group < 0) {
      sprintf(AbortString,"%s: H5Oopen failure on object '%.72s'.",nomsub,namp);
      xabort_c(AbortString);
  }
  hsize_t nobj;
  int otype;
  char memb_name[MAX_NAME];
  iretcd = H5Gget_num_objs(group, &nobj);
  if (iretcd < 0) {
      sprintf(AbortString,"%s: H5Gget_num_objs failure on object '%.72s'.",nomsub,namp);
      xabort_c(AbortString);
  }
  *ndsets = 0;
  for (idx = 0; idx < nobj; idx++) {
    H5Gget_objname_by_idx(group, (hsize_t)idx, memb_name, (size_t)MAX_NAME);
    otype =  H5Gget_objtype_by_idx(group, (size_t)idx);
    switch(otype) {
    case H5G_GROUP:
      strncpy (idata+MAX_NAME*(*ndsets), memb_name, MAX_NAME);
      (*ndsets)++;
      break;
    }
  }
  H5Oclose(group);
}

void hdf5_info_c(hid_t *ifile, const char *namp, int_32 *rank, int_32 *type, int_32 *nbyte, int_32 *dimsr) {
/*
 *----------------------------------------------------------------------
 *
 * Find dataset information.
 *
 * input parameters:
 *   ifile : hdf5 file identificator.
 *    namp : character name of the dataset.
 *
 * output parameter:
 *   rank  : number of dimensions.
 *   type  : type of dataset.
 *   nbyte : number of bytes per value.
 *   dimsr : shape.
 *
 *----------------------------------------------------------------------
 */
  char *nomsub="hdf5_info_c";
  hid_t loc_id=*ifile;
  hsize_t dimsr_t[10];
  int i;
  
  hid_t dataset = -1;
  H5E_BEGIN_TRY {
    dataset = H5Dopen(loc_id,namp,H5P_DEFAULT);
  } H5E_END_TRY;
  if (dataset < 0) {
      *rank=0;
      *type=99;
      *nbyte=0;
      return;
  }
  hid_t filespace = H5Dget_space (dataset);
  *rank = H5Sget_simple_extent_ndims (filespace);
  if (*rank < 0) {
      sprintf(AbortString,"%s: H5Sget_simple_extent_ndims failure on object '%.72s'.",nomsub,namp);
      xabort_c(AbortString);
  } else if (*rank > 10) {
      sprintf(AbortString,"%s: the object '%.72s' has rank= %d > 10.",nomsub,namp,*rank);
      xabort_c(AbortString);
  }
  hid_t htype = H5Dget_type(dataset);
  if (htype < 0) {
      sprintf(AbortString,"%s: H5Dget_type failure on object '%.72s'.",nomsub,namp);
      xabort_c(AbortString);
  }
  size_t precision = 0;
  switch (H5Tget_class(htype)) {
  case H5T_INTEGER: 
       *type=1;
       precision = H5Tget_size(htype);
       break;
  case H5T_FLOAT: 
       /*precision = H5Tget_precision(htype);*/
       precision = H5Tget_size(htype);
       if (precision == 4) {
         *type=2;
       } else if (precision == 8) {
         *type=4;
       } else {
         *type=7;
       }
       break;
  case H5T_STRING: 
       *type=3;
       precision = H5Tget_size(htype);
       break;
  default:
       *type=7;
  }
  *nbyte = precision;
  iretcd = H5Sget_simple_extent_dims(filespace, dimsr_t, NULL);
  if (iretcd == -1) {
      sprintf(AbortString,"%s: H5Sget_simple_extent_dims failure on object '%.72s'.",nomsub,namp);
      xabort_c(AbortString);
  }
  for(i=0; i<*rank; i++) dimsr[i] = (int_32)dimsr_t[(*rank)-i-1];
  
  H5Tclose(htype);
  H5Dclose(dataset);
}

int_32 hdf5_group_exists_c(hid_t *ifile, const char *namp)
/*
 *-----------------------------------------------------------------------
 *
 * Test for existence of a group.
 *
 * input parameters:
 *   ifile : hdf5 file identificator.
 *    namp : character name of the group.
 *
 * output parameter:
 *  hdf5_group_exists_c : =0: the group does exists; =1: does not exists.
 *
 *-----------------------------------------------------------------------
 */
{
  hid_t loc_id=*ifile;
  H5G_stat_t statbuf;
H5E_BEGIN_TRY
  if(strlen(namp)==0) {
    iretcd = H5Gget_objinfo(loc_id, "/", 0, &statbuf);
  } else {
    iretcd = H5Gget_objinfo(loc_id, namp, 0, &statbuf);
  }
H5E_END_TRY
  if (iretcd >= 0) {
    if (statbuf.type == H5G_GROUP) return 1;
  }
  return 0;
}

void hdf5_create_group_c(hid_t *ifile, const char *namp)
/*
 *-----------------------------------------------------------------------
 *
 * Create a group.
 *
 * input parameters:
 *   ifile : hdf5 file identificator.
 *    namp : character name of the group to create.
 *
 *-----------------------------------------------------------------------
 */
{
  char *nomsub="hdf5_create_group_c";
  hid_t loc_id=*ifile;
  H5G_stat_t statbuf;
H5E_BEGIN_TRY
  iretcd = H5Gget_objinfo(loc_id, namp, 0, &statbuf);
H5E_END_TRY
  if (iretcd >= 0) {
    if (statbuf.type == H5G_GROUP) return;
  }
  hid_t groupid = H5Gcreate(loc_id, namp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (groupid < 0)  {
    sprintf(AbortString,"%s: H5Gcreate failure on object '%.72s'.",nomsub,namp);
    xabort_c(AbortString);
  }
  H5Gclose(groupid);
}

void hdf5_delete_c(hid_t *ifile, const char *namp)
/*
 *-----------------------------------------------------------------------
 *
 * Delete a group or a dataset.
 *
 * input parameters:
 *   ifile : hdf5 file identificator.
 *    namp : character name of the group or dataset to delete.
 *
 *-----------------------------------------------------------------------
 */
{
  char *nomsub="hdf5_delete_c";
  hid_t loc_id=*ifile;
  iretcd = H5Ldelete(loc_id,namp,H5P_DEFAULT);
  if (iretcd < 0) {
    sprintf(AbortString,"%s: HDF5 delete failure. iretcd=%d\n", nomsub, iretcd);
    xabort_c(AbortString);
  }
}

void hdf5_copy_c(hid_t *ifile_s, const char *namp_s, hid_t *ifile_d, const char *namp_d)
/*
 *-----------------------------------------------------------------------
 *
 * Copy a group or a dataset from one location to another. The source and
 * destination need not be in the same file.
 *
 * input parameters:
 *   ifile_s : hdf5 source file identificator.
 *    namp_s : character name of the source group or dataset to copy.
 *
 * output parameters:
 *   ifile_d : hdf5 destination file identificator.
 *    namp_d : character name of the destination group or dataset.
 *
 *-----------------------------------------------------------------------
 */
{
  char *nomsub="hdf5_copy_c";
  if(hdf5_group_exists_c(ifile_s, namp_s)) {
    /* namp_s is a group */
    hid_t loc_id_s=*ifile_s;
    hid_t loc_id_d=*ifile_d;
    iretcd = H5Ocopy(loc_id_s,namp_s,loc_id_d,namp_d,H5P_DEFAULT,H5P_DEFAULT);
    if (iretcd < 0) {
      sprintf(AbortString,"%s: HDF5 copy failure. iretcd=%d\n", nomsub, iretcd);
      xabort_c(AbortString);
    }
  } else {
    /* namp_s is a dataset */
    int i;
    int_32 rank, type, nbyte, dimsr[10], dimsr2[10], length;
    hdf5_info_c(ifile_s, namp_s, &rank, &type, &nbyte, dimsr);
    length = 1;
    for(i=0; i<rank; i++) {
      dimsr2[i] = (int_32)dimsr[(rank)-i-1];
      length = length*dimsr2[i];
    }
    if (type == 1) {
      int_32 *idata = (int_32 *)malloc(sizeof(int_32)*length);
      hdf5_read_data_int_c(ifile_s, namp_s, idata);
      hdf5_write_data_int_c(ifile_d, namp_s, rank, dimsr2, idata);
      free(idata);
    } else if (type == 2) {
      float *rdata = (float *)malloc(sizeof(float)*length);
      hdf5_read_data_real4_c(ifile_s, namp_s, rdata);
      hdf5_write_data_real4_c(ifile_d, namp_s, rank, dimsr2, rdata);
      free(rdata);
    } else if (type == 3) {
      char *ihdata = (char *)malloc(nbyte*length);
      hdf5_read_data_string_c(ifile_s, namp_s, ihdata);
      hdf5_write_data_string_c(ifile_d, namp_s, rank, nbyte, dimsr2, ihdata);
      free(ihdata);
    } else if (type == 4) {
      double *rddata = (double *)malloc(sizeof(double)*length);
      hdf5_read_data_real8_c(ifile_s, namp_s, rddata);
      hdf5_write_data_real8_c(ifile_d, namp_s, rank, dimsr2, rddata);
      free(rddata);
    } else {
      sprintf(AbortString,"%s: dataset '%.60s' has the wrong type(2).", nomsub, namp_s);
      xabort_c(AbortString);
    }
  }
}

void hdf5_read_data_int_c(hid_t *ifile, const char *namp, int_32 *idata) {
/*
 *----------------------------------------------------------------------
 *
 * Copy an integer dataset from hdf5 into memory.
 *
 * input parameters:
 *   ifile : hdf5 file identificator.
 *    namp : character name of the dataset.
 *
 * output parameter:
 *   idata : information elements.
 *
 *----------------------------------------------------------------------
 */
  char *nomsub="hdf5_read_data_int_c";
  hid_t loc_id=*ifile;
  hid_t dataset = H5Dopen(loc_id,namp,H5P_DEFAULT);
  iretcd = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, idata);
  if (iretcd != 0) {
      sprintf(AbortString,"%s: the object '%.72s' cannot be readed.",nomsub,namp);
      xabort_c(AbortString);
  }
  H5Dclose(dataset);
}

void hdf5_read_data_real4_c(hid_t *ifile, const char *namp, float *rdata) {
/*
 *----------------------------------------------------------------------
 *
 * Copy a real(4) dataset from hdf5 into memory.
 *
 * input parameters:
 *   ifile : hdf5 file identificator.
 *    namp : character name of the dataset.
 *
 * output parameter:
 *   rdata : information elements.
 *
 *----------------------------------------------------------------------
 */
  char *nomsub="hdf5_read_data_real4_c";
  hid_t loc_id=*ifile;
  hid_t dataset = H5Dopen(loc_id,namp,H5P_DEFAULT);
  iretcd = H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata);
  if (iretcd != 0) {
      sprintf(AbortString,"%s: the object '%.72s' cannot be readed.",nomsub,namp);
      xabort_c(AbortString);
  }
  H5Dclose(dataset);
}

void hdf5_read_data_real8_c(hid_t *ifile, const char *namp, double *rdata) {
/*
 *----------------------------------------------------------------------
 *
 * Copy a real(8) dataset from hdf5 into memory.
 *
 * input parameters:
 *   ifile : hdf5 file identificator.
 *    namp : character name of the dataset.
 *
 * output parameter:
 *   rdata : information elements.
 *
 *----------------------------------------------------------------------
 */
  char *nomsub="hdf5_read_data_real8_c";
  hid_t loc_id=*ifile;
  hid_t dataset = H5Dopen(loc_id,namp,H5P_DEFAULT);
  iretcd = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata);
  if (iretcd != 0) {
      sprintf(AbortString,"%s: the object '%.72s' cannot be readed.",nomsub,namp);
      xabort_c(AbortString);
  }
  H5Dclose(dataset);
}

void hdf5_read_data_string_c(hid_t *ifile, const char *namp, char *idata) {
/*
 *----------------------------------------------------------------------
 *
 * Copy a string dataset from hdf5 into memory.
 *
 * input parameters:
 *   ifile : hdf5 file identificator.
 *    namp : character name of the dataset.
 *
 * output parameter:
 *   idata : information elements.
 *
 *----------------------------------------------------------------------
 */
  char *nomsub="hdf5_read_data_string_c";
  hid_t loc_id=*ifile;
  
  hid_t dataset = H5Dopen(loc_id,namp,H5P_DEFAULT);
  hid_t datatype = H5Dget_type(dataset);
  if (datatype < 0) {
      sprintf(AbortString,"%s: H5Dget_type failure on object '%.72s'.",nomsub,namp);
      xabort_c(AbortString);
  }
  iretcd = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, idata);
  if (iretcd != 0) {
      sprintf(AbortString,"%s: the object '%.72s' cannot be readed.",nomsub,namp);
      xabort_c(AbortString);
  }
  H5Tclose(datatype);
  H5Dclose(dataset);
}

void hdf5_write_data_int_c(hid_t *ifile, const char *namp, int_32 rank, int_32 *dimsr, int_32 *idata) {
/*
 *----------------------------------------------------------------------
 *
 * Copy an integer array from memory into a hdf5 dataset
 *
 * input parameters:
 *   ifile : hdf5 file identificator.
 *    namp : character name of the dataset.
 *    rank : number of dimensions.
 *   dimsr : number of information along each dimension.
 *   idata : information elements.
 *
 *----------------------------------------------------------------------
 */
  char *nomsub="hdf5_write_data_int_c";
  hid_t loc_id=*ifile;
  int i;
  hsize_t dimsr_t[10];
  if (rank > 10) {
      sprintf(AbortString,"%s: rank > 10 on object '%.72s'.",nomsub,namp);
      xabort_c(AbortString);
  }
  hid_t dataset = -1;
  H5E_BEGIN_TRY {
  dataset = H5Dopen(loc_id,namp,H5P_DEFAULT);
  } H5E_END_TRY;
  if (dataset  >= 0) {
    iretcd = H5Ldelete(loc_id,namp,H5P_DEFAULT);
    if (iretcd < 0) {
      sprintf(AbortString,"%s: HDF5 delete failure. iretcd=%d\n", nomsub, iretcd);
      xabort_c(AbortString);
    }
    H5Dclose(dataset);
  }
  for (i = 0; i < rank; ++i) dimsr_t[i] = dimsr[i];
  hid_t dataspace = H5Screate_simple(rank, dimsr_t, NULL); 
  hid_t datatype = H5Tcopy(H5T_NATIVE_INT);
  dataset = H5Dcreate(loc_id, namp, datatype, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset < 0) {
      sprintf(AbortString,"%s: H5Dcreate failure on object '%.72s'.",nomsub,namp);
      xabort_c(AbortString);
  }
  iretcd = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, idata);
  if (iretcd != 0) {
      sprintf(AbortString,"%s: the object '%.72s' cannot be written.",nomsub,namp);
      xabort_c(AbortString);
  }
  H5Sclose(dataspace);
  H5Tclose(datatype);
  H5Dclose(dataset);
}

void hdf5_write_data_real4_c(hid_t *ifile, const char *namp, int_32 rank, int_32 *dimsr, float *rdata) {
/*
 *----------------------------------------------------------------------
 *
 * Copy a real(4) array from memory into a hdf5 dataset
 *
 * input parameters:
 *   ifile : hdf5 file identificator.
 *    namp : character name of the dataset.
 *    rank : number of dimensions.
 *   dimsr : number of information along each dimension.
 *   rdata : information elements.
 *
 *----------------------------------------------------------------------
 */
  char *nomsub="hdf5_write_data_real4_c";
  hid_t loc_id=*ifile;
  int i;
  hsize_t dimsr_t[10];
  if (rank > 10) {
      sprintf(AbortString,"%s: rank > 10 on object '%.72s'.",nomsub,namp);
      xabort_c(AbortString);
  }
  hid_t dataset = -1;
  H5E_BEGIN_TRY {
  dataset = H5Dopen(loc_id,namp,H5P_DEFAULT);
  } H5E_END_TRY;
  if (dataset  >= 0) {
    iretcd = H5Ldelete(loc_id,namp,H5P_DEFAULT);
    if (iretcd < 0) {
      sprintf(AbortString,"%s: HDF5 delete failure. iretcd=%d\n", nomsub, iretcd);
      xabort_c(AbortString);
    }
    H5Dclose(dataset);
  }
  for (i = 0; i < rank; ++i) dimsr_t[i] = dimsr[i];
  hid_t dataspace = H5Screate_simple(rank, dimsr_t, NULL); 
  hid_t datatype = H5Tcopy(H5T_NATIVE_FLOAT);
  dataset = H5Dcreate(loc_id, namp, datatype, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset < 0) {
      sprintf(AbortString,"%s: H5Dcreate failure on object '%.72s'.",nomsub,namp);
      xabort_c(AbortString);
  }
  iretcd = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata);
  if (iretcd != 0) {
      sprintf(AbortString,"%s: the object '%.72s' cannot be written.",nomsub,namp);
      xabort_c(AbortString);
  }
  H5Sclose(dataspace);
  H5Tclose(datatype);
  H5Dclose(dataset);
}

void hdf5_write_data_real8_c(hid_t *ifile, const char *namp, int_32 rank, int_32 *dimsr, double *rdata) {
/*
 *----------------------------------------------------------------------
 *
 * Copy a real(8) array from memory into a hdf5 dataset
 *
 * input parameters:
 *   ifile : hdf5 file identificator.
 *    namp : character name of the dataset.
 *    rank : number of dimensions.
 *   dimsr : number of information along each dimension.
 *   rdata : information elements.
 *
 *----------------------------------------------------------------------
 */
  char *nomsub="hdf5_write_data_real8_c";
  hid_t loc_id=*ifile;
  int i;
  hsize_t dimsr_t[10];
  if (rank > 10) {
      sprintf(AbortString,"%s: rank > 10 on object '%.72s'.",nomsub,namp);
      xabort_c(AbortString);
  }
  hid_t dataset = -1;
  H5E_BEGIN_TRY {
  dataset = H5Dopen(loc_id,namp,H5P_DEFAULT);
  } H5E_END_TRY;
  if (dataset  >= 0) {
    iretcd = H5Ldelete(loc_id,namp,H5P_DEFAULT);
    if (iretcd < 0) {
      sprintf(AbortString,"%s: HDF5 delete failure. iretcd=%d\n", nomsub, iretcd);
      xabort_c(AbortString);
    }
    H5Dclose(dataset);
  }
  for (i = 0; i < rank; ++i) dimsr_t[i] = dimsr[i];
  hid_t dataspace = H5Screate_simple(rank, dimsr_t, NULL); 
  hid_t datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
  dataset = H5Dcreate(loc_id, namp, datatype, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset < 0) {
      sprintf(AbortString,"%s: H5Dcreate failure on object '%.72s'.",nomsub,namp);
      xabort_c(AbortString);
  }
  iretcd = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata);
  if (iretcd != 0) {
      sprintf(AbortString,"%s: the object '%.72s' cannot be written.",nomsub,namp);
      xabort_c(AbortString);
  }
  H5Sclose(dataspace);
  H5Tclose(datatype);
  H5Dclose(dataset);
}

void hdf5_write_data_string_c(hid_t *ifile, const char *namp, int_32 rank, int_32 len, int_32 *dimsr, char *idata) {
/*
 *----------------------------------------------------------------------
 *
 * Copy a character array from memory into a hdf5 dataset
 *
 * input parameters:
 *   ifile : hdf5 file identificator.
 *    namp : character name of the dataset.
 *    rank : number of dimensions.
 *     len : length of a string element in the array (in bytes).
 *   dimsr : number of information along each dimension.
 *   idata : information elements.
 *
 *----------------------------------------------------------------------
 */
  char *nomsub="hdf5_write_data_string_c";
  hid_t loc_id=*ifile;
  int i;
  hsize_t dimsr_t[10];
  if (rank > 10) {
      sprintf(AbortString,"%s: rank > 10 on object '%.72s'.",nomsub,namp);
      xabort_c(AbortString);
  }
  hid_t dataset = -1;
  H5E_BEGIN_TRY {
  dataset = H5Dopen(loc_id,namp,H5P_DEFAULT);
  } H5E_END_TRY;
  if (dataset  >= 0) {
    iretcd = H5Ldelete(loc_id,namp,H5P_DEFAULT);
    if (iretcd < 0) {
      sprintf(AbortString,"%s: HDF5 delete failure. iretcd=%d\n", nomsub, iretcd);
      xabort_c(AbortString);
    }
    H5Dclose(dataset);
  }
  for (i = 0; i < rank; ++i) dimsr_t[i] = dimsr[i];
  hid_t dataspace = H5Screate_simple(rank, dimsr_t, NULL);
  hid_t datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, len * sizeof(char));
  dataset = H5Dcreate(loc_id, namp, datatype , dataspace, H5P_DEFAULT,
      H5P_DEFAULT, H5P_DEFAULT);
  if (dataset < 0) {
      sprintf(AbortString,"%s: H5Dcreate failure on object '%.72s'.",nomsub,namp);
      xabort_c(AbortString);
  }
  iretcd = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, idata);
  if (iretcd != 0) {
      sprintf(AbortString,"%s: the object '%.72s' cannot be written.",nomsub,namp);
      xabort_c(AbortString);
  }
  H5Sclose(dataspace);
  H5Tclose(datatype);
  H5Dclose(dataset);
}

herr_t print_info(hid_t loc_id, const char *name, void *opdata) {
    H5G_stat_t statbuf;
    char* ctype[]={"","INTEGER","REAL","CHARACTER","DOUBLE PRECISION",
                   "LOGICAL","COMPLEX","UNDEFINED"};
    /*
     * Get type of the object and display its name and type.
     * The name of the object is passed to this function by 
     * the Library.
     */
    H5Gget_objinfo(loc_id, name, 0, &statbuf);
    int rank, type, nbyte;
    int_32 dimsr[6];
    switch (statbuf.type) {
    case H5G_GROUP: 
         printf(" '%-72s' GROUP \n", name);
         break;
    case H5G_DATASET: 
         hdf5_info_c(&loc_id, name, &rank, &type, &nbyte, dimsr);
         if (rank == 1) {
           printf(" '%-72s' %-16s   %-10d %d \n", name,ctype[type],nbyte,dimsr[0]);
         } else if (rank == 2) {
           printf(" '%-72s' %-16s   %-10d %d  %d\n", name,ctype[type],nbyte,dimsr[0],dimsr[1]);
         } else if (rank == 3) {
           printf(" '%-72s' %-16s   %-10d %d  %d  %d\n", name,ctype[type],nbyte,dimsr[0],dimsr[1],dimsr[2]);
         } else if (rank == 4) {
           printf(" '%-72s' %-16s   %-10d %d  %d  %d  %d\n", name,ctype[type],nbyte,dimsr[0],dimsr[1],
           dimsr[2],dimsr[3]);
         } else if (rank == 5) {
           printf(" '%-72s' %-16s   %-10d %d  %d  %d  %d  %d\n", name,ctype[type],nbyte,dimsr[0],
           dimsr[1],dimsr[2],dimsr[3],dimsr[4]);
         } else if (rank == 6) {
           printf(" '%-72s' %-16s   %-10d %d  %d  %d  %d  %d  %d\n", name,ctype[type],nbyte,dimsr[0],
           dimsr[1],dimsr[2],dimsr[3],dimsr[4],dimsr[5]);
         }
         break;
    case H5G_TYPE: 
         printf(" '%-72s' NAMED DATATYPE \n", name);
         break;
    default:
         printf(" '%-72s' UNKNOWN \n", name);
    }
    fflush(stdout);
    return 0;
 }
#endif /* defined(HDF5_LIB) */
