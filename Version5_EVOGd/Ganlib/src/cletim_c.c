/*****************************************/
/*             CLE-2000 API              */
/*     AUTHOR: A. Hebert ; 16/07/10      */
/*****************************************/

#include <stdlib.h>
#include <time.h>
#include "cle2000.h"

#ifdef _OPENMP
#include <omp.h>
void cletim_c(double *sec){
   *sec = omp_get_wtime();
}
#else
void cletim_c(double *sec){
   long value = (long) clock();
   *sec = (double) (value / CLOCKS_PER_SEC);
}
#endif
