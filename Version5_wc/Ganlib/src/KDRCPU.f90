subroutine KDRCPU(cpusec)
!
!-----------------------------------------------------------------------
!
!Purpose:
! system clock support.
!
!Copyright:
! Copyright (C) 2002 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s): A. Hebert
!
!Parameters: output
!  cpusec : number of seconds elapsed since the first call to KDRCPU.
!
!-----------------------------------------------------------------------
!
   integer :: itloc,irate
   double precision :: dtloc
   integer,save :: isave=0,itloc0
   double precision,save :: dtloc0
!
   if(isave==0) then
      call CLETIM(dtloc0)
      isave=1
   endif
   call CLETIM(dtloc)
   cpusec=real(dtloc-dtloc0)
   return
end subroutine KDRCPU
