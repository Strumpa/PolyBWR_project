
/**********************************/
/* C API for lcm object support   */
/* (auxiliary functions)          */
/* author: A. Hebert (30/04/2002) */
/**********************************/

/*
   Copyright (C) 2002 Ecole Polytechnique de Montreal

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
 */

#include <stdlib.h>
#include <string.h>
#include "lcm.h"

static char AbortString[132];

void lcmput_c(lcm **iplist,const char *namp,int_32 ilong,int_32 itype,int_32 *idata)
/*
 *----------------------------------------------------------------------
 *
 * copy a block of data from memory into a table.
 *
 * input parameters:
 *  iplist : address of the table.
 *    namp : character*12 name of the current block.
 *   ilong : number of information elements stored in the current block.
 *   itype : type of information elements stored in the current block.
 *           0: directory                1: integer
 *           2: single precision         3: character*4
 *           4: double precision         5: logical
 *           6: complex                 99: undefined
 *   idata : information elements.
 *
 *----------------------------------------------------------------------
 */
{
   if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      xsmput_c((xsm **)iplist,namp,ilong,itype,idata);
   } else {
      int_32 i, *iofdat;
      int_32 jlong = ilong;
      if (itype == 4 || itype == 6) jlong = 2*ilong;
      iofdat = (int_32 *)malloc(jlong*sizeof(int_32)); /* setara_c(jlong); */
      for (i = 0; i < jlong; ++i) iofdat[i] = idata[i];
      lcmppd_c(iplist,namp,ilong,itype,iofdat);
   }
}

void lcmpdl_c(lcm **iplist,int_32 iset,int_32 ilong,int_32 itype,int_32 *idata)
/*
 *----------------------------------------------------------------------
 *
 * copy a block of data from memory into a list.
 *
 * input parameters:
 *  iplist : address of the list.
 *    iset : position of the specific element.
 *   ilong : number of information elements stored in the current block.
 *   itype : type of information elements stored in the current block.
 *           0: directory                1: integer
 *           2: single precision         3: character*4
 *           4: double precision         5: logical
 *           6: complex                 99: undefined
 *   idata : information elements.
 *
 *----------------------------------------------------------------------
 */
{
   if ((*iplist)->header == 200) {
      /* USE A XSM FILE. */
      xsm *ipxsm = (xsm *)*iplist + iset;
      xsmput_c(&ipxsm," ",ilong,itype,idata);
   } else {
      int_32 i, *iofdat;
      int_32 jlong = ilong;
      if (itype == 4 || itype == 6) jlong = 2*ilong;
      iofdat = (int_32 *)malloc(jlong*sizeof(int_32)); /* setara_c(jlong); */
      for (i = 0; i < jlong; ++i) iofdat[i] = idata[i];
      lcmppl_c(iplist,iset,ilong,itype,iofdat);
   }
}

void lcmpcd_c(lcm **iplist,const char *namp,int_32 ilong,char *hdata[])
/*
 *----------------------------------------------------------------------
 *
 * copy an array of c string variables from memory into a table.
 *
 * input parameters:
 *  iplist : address of the table.
 *    namp : character*12 name of the block.
 *   ilong : dimension of the string array.
 *   hdata : array of ilong strings.
 *
 *----------------------------------------------------------------------
 */
{
   int_32 iset;
   lcm *jplist;
   jplist = lcmlid_c(iplist, namp, ilong);
   for (iset=0; iset<ilong; iset++) {
      int_32 i, ilen, *iofset;
      ilen = (strlen(hdata[iset]) + 4 ) / 4;
      iofset = (int_32 *)malloc(ilen*sizeof(int_32)); /* setara_c(ilen); */
      for (i=0; i<ilen; i++) strncpy ((char *)(iofset+i), hdata[iset]+4*i, 4);
      lcmppl_c(&jplist, iset, ilen, 3, iofset);
   }
}

void lcmgcd_c(lcm **iplist,const char *namp,char *hdata[])
/*
 *-----------------------------------------------------------------------
 *
 * copy an array of c string variables from a table into memory.
 *
 * input parameters:
 *  iplist : address of the table.
 *    namp : character*12 name of the existing block.
 *
 * output parameter:
 *   hdata : array of ilong strings (allocated by lcmgcd_c).
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmgcd_c";
   int_32 iset, ilong, itylcm;
   lcm *jplist;
   lcmlen_c(iplist, namp, &ilong, &itylcm);
   if (itylcm != 10) {
      sprintf(AbortString,"%s: LIST EXPECTED.",nomsub);
      xabort_c(AbortString);
   }
   jplist = lcmgid_c(iplist, namp);
   for (iset=0; iset<ilong; iset++) {
      int_32 j, ilcmg, *iass;
      lcmlel_c(&jplist, iset, &ilcmg, &itylcm);
      iass = (int_32 *)malloc(ilcmg*sizeof(int_32)); /* setara_c(ilcmg); */
      lcmgdl_c(&jplist, iset, iass);
      hdata[iset] = (char *)malloc((int)4*ilcmg+1);
      for (j=0; j<ilcmg; j++) strncpy ((hdata[iset]+4*j),(char *) (iass + j), 4);
      hdata[iset][4*ilcmg]=' ';
      free(iass); /* rlsara_c(iass); */
      for(j=4*ilcmg; j>0; j--) {
         if (hdata[iset][j] != ' ') break;
         hdata[iset][j]='\0';
      }
   }
}

void lcmpcl_c(lcm **iplist,int_32 iset,int_32 ilong,char *hdata[])
/*
 *----------------------------------------------------------------------
 *
 * copy an array of c string variables from memory into a list.
 *
 * input parameters:
 *  iplist : address of the table.
 *    iset : position of the block in the list.
 *   ilong : dimension of the character variable.
 *   hdata : array of ilong strings.
 *
 *----------------------------------------------------------------------
 */
{
   int_32 jset;
   lcm *jplist;
   jplist=lcmlil_c(iplist, iset, ilong);
   for (jset=0; jset<ilong; jset++) {
      int_32 i, ilen, *iofset;
      ilen = (strlen(hdata[jset]) + 4 ) / 4;
      iofset = (int_32 *)malloc(ilen*sizeof(int_32)); /* setara_c(ilen); */
      for (i=0; i<ilen; i++) strncpy ((char *)(iofset+i), hdata[jset]+4*i, 4);
      lcmppl_c(&jplist, jset, ilen, 3, iofset);
   }
}

void lcmgcl_c(lcm **iplist,int_32 iset,char *hdata[])
/*
 *-----------------------------------------------------------------------
 *
 * copy an array of c string variables from a list into memory.
 *
 * input parameters:
 *  iplist : address of the table.
 *    iset : position of the block in the list.
 *
 * output parameter:
 *   hdata : array of ilong strings (allocated by lcmgcl_c).
 *
 *-----------------------------------------------------------------------
 */
{
   char *nomsub="lcmgcl_c";
   int_32 jset, ilong, itylcm;
   lcm *jplist;
   lcmlel_c(iplist, iset, &ilong, &itylcm);
   if (itylcm != 10) {
      sprintf(AbortString,"%s: LIST EXPECTED.",nomsub);
      xabort_c(AbortString);
   }
   jplist = lcmgil_c(iplist, iset);
   for (jset=0; jset<ilong; jset++) {
      int_32 j, ilcmg, *iass;
      lcmlel_c(&jplist, jset, &ilcmg, &itylcm);
      iass = (int_32 *)malloc(ilcmg*sizeof(int_32)); /* setara_c(ilcmg); */
      lcmgdl_c(&jplist, jset, iass);
      hdata[jset] = (char *)malloc((int)4*ilcmg+1);
      for (j=0; j<ilcmg; j++) strncpy ((hdata[jset]+4*j),(char *) (iass + j), 4);
      hdata[jset][4*ilcmg]=' ';
      free(iass); /* rlsara_c(iass); */
      for(j=4*ilcmg; j>0; j--) {
         if (hdata[jset][j] != ' ') break;
         hdata[jset][j]='\0';
      }
   }
}

void lcmpsd_c(lcm **iplist,const char *namp,char *hdata)
/*
 *----------------------------------------------------------------------
 *
 * copy a single c string variable from memory into a table.
 *
 * input parameters:
 *  iplist : address of the table.
 *    namp : character*12 name of the block.
 *   hdata : c string.
 *
 *----------------------------------------------------------------------
 */
{
   int_32 i, ilong, *iofset;
   ilong = (strlen(hdata) + 4 ) / 4;
   iofset = (int_32 *)malloc(ilong*sizeof(int_32)); /* setara_c(ilong); */
   for (i=0; i<ilong; i++) strncpy ((char *)(iofset+i), hdata+4*i, 4);
   lcmppd_c(iplist, namp, ilong, 3, iofset);
}

char * lcmgsd_c(lcm **iplist,const char *namp)
/*
 *-----------------------------------------------------------------------
 *
 * copy a single c string variable from a table into memory.
 *
 * input parameters:
 *  iplist : address of the table.
 *    namp : character*12 name of the existing block.
 *
 * output parameter:
 *   lcmgsd_c : c string.
 *
 *-----------------------------------------------------------------------
 */
{
   static char nomstatic[133];
   char *nomsub="lcmgsd_c";
   int_32 i, ilong, itylcm, *iass;
   lcmlen_c(iplist, namp, &ilong, &itylcm);
   if (itylcm != 3) {
      sprintf(AbortString,"%s: CHARACTER DATA EXPECTED.",nomsub);
      xabort_c(AbortString);
   } else if (ilong*4 > 132) {
      sprintf(AbortString,"%s: CHARACTER DATA OVERFLOW.",nomsub);
      xabort_c(AbortString);
   }

   iass = (int_32 *)malloc(ilong*sizeof(int_32)); /* setara_c(ilong); */
   lcmget_c(iplist, namp, iass);
   for (i=0; i<ilong; i++) strncpy ((nomstatic+4*i),(char *) (iass+i), 4);
   nomstatic[ilong*4] = ' ';
   free(iass); /* rlsara_c(iass); */
   for(i=ilong*4; i>0; i--) {
      if(nomstatic[i] != ' ') break;
      nomstatic[i] = '\0';
   }
   return nomstatic;
}

void lcmpsl_c(lcm **iplist,int_32 iset,char *hdata)
/*
 *----------------------------------------------------------------------
 *
 * copy a single c string variable from memory into a list.
 *
 * input parameters:
 *  iplist : address of the table.
 *    iset : position of the block in the list.
 *   hdata : c string.
 *
 *----------------------------------------------------------------------
 */
{
   int_32 i, ilong, *iofset;
   ilong = (strlen(hdata) + 4 ) / 4;
   iofset = (int_32 *)malloc(ilong*sizeof(int_32)); /* setara_c(ilong); */
   for (i=0; i<ilong; i++) strncpy ((char *)(iofset+i), hdata+4*i, 4);
   lcmppl_c(iplist, iset, ilong, 3, iofset);
}

char * lcmgsl_c(lcm **iplist,int_32 iset)
/*
 *-----------------------------------------------------------------------
 *
 * copy a single c string variable from a list into memory.
 *
 * input parameters:
 *  iplist : address of the table.
 *    iset : position of the block in the list.
 *
 * output parameter:
 *   lcmgsd_c : c string.
 *
 *-----------------------------------------------------------------------
 */
{
   static char nomstatic[133];
   char *nomsub="lcmgsl_c";
   int_32 i, ilong, itylcm, *iass;
   lcmlel_c(iplist, iset, &ilong, &itylcm);
   if (itylcm != 3) {
      sprintf(AbortString,"%s: CHARACTER DATA EXPECTED.",nomsub);
      xabort_c(AbortString);
   } else if (ilong*4 > 132) {
      sprintf(AbortString,"%s: CHARACTER DATA OVERFLOW.",nomsub);
      xabort_c(AbortString);
   }

   iass = (int_32 *)malloc(ilong*sizeof(int_32)); /* setara_c(ilong); */
   lcmgdl_c(iplist, iset, iass);
   for (i=0; i<ilong; i++) strncpy ((nomstatic+4*i),(char *) (iass+i), 4);
   nomstatic[ilong*4] = ' ';
   free(iass); /* rlsara_c(iass); */
   for(i=ilong*4; i>0; i--) {
      if(nomstatic[i] != ' ') break;
      nomstatic[i] = '\0';
   }
   return nomstatic;
}
