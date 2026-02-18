!
!-----------------------------------------------------------------------
!
!Purpose:
! Fortran-2003 bindings for lcm -- part 3.
! Support of character arrays.
!
!Copyright:
! Copyright (C) 2009 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s): A. Hebert
!
!Parameters: input
! iplcm   pointer to the LCM object.
! name    character name of the LCM node.
! ipos    heterogeneous list index.
! leng    length of each character string in the array carr.
! nlin    dimension of array carr.
!
!Parameters: input/output
! carr    array of character strings.
!
!-----------------------------------------------------------------------
!
module LCMTLC
   use LCMMOD
   private
   public :: LCMGTC, LCMPTC, LCMGLC, LCMPLC
   interface LCMGTC
      ! recover a string array from an associative table
      MODULE PROCEDURE LCMGTC_S0, LCMGTC_S1, LCMGTC_S2
   end interface
   interface LCMPTC
      ! store a string array from an associative table
      MODULE PROCEDURE LCMPTC_S0, LCMPTC_S1, LCMPTC_S2
   end interface
   interface LCMGLC
      ! recover a string array from an heterogeneous list
      MODULE PROCEDURE LCMGLC_S0, LCMGLC_S1, LCMGLC_S2
   end interface
   interface LCMPLC
      ! store a string array from an heterogeneous list
      MODULE PROCEDURE LCMPLC_S0, LCMPLC_S1, LCMPLC_S2
   end interface
contains
subroutine LCMGTC_S0(iplcm,name,leng,carr)
   use, intrinsic :: iso_c_binding
   use LCMAUX
   !----
   !  Subroutine arguments
   !----
   type(c_ptr), intent(in) :: iplcm
   integer, intent(in) :: leng
   character(len=*), intent(in) :: name
   character(len=*), intent(out) :: carr
   !----
   !  Local variables
   !----
   character(len=13)  :: fmt
   character(len=12) :: text12
   integer,allocatable :: ibase(:)
   !
   fmt='(           )'
   allocate(ibase((leng+3)/4))
   !----
   !  Read from LCM object
   !----
   call LCMLEN(iplcm,name,ilong,itylcm)
   if(ilong == 0) then
      call LCMLIB(iplcm)
      text12=name
      call XABORT('LCMGTC: record '//text12//' not found.')
   endif
   call LCMGET(iplcm,name,ibase)
   !----
   !  Define format and dimensions for conversion
   !----
   n=leng/4
   m=mod(leng,4)
   if(n > 0)write(fmt(2:9),'(i6,''a4'')') n
   if(m > 0)then
      if(n.gt.0)fmt(10:10)=','
      write(fmt(11:12),'(''a'',i1)') m
      n=n+1
   endif
   !----
   !  Convert integers to string
   !----
   write(carr(1:leng),fmt) (ibase(j),j=1,n)
   deallocate(ibase)
end subroutine LCMGTC_S0
!
subroutine LCMGTC_S1(iplcm,name,leng,nlin,carr)
   use, intrinsic :: iso_c_binding
   use LCMAUX
   !----
   !  Subroutine arguments
   !----
   type(c_ptr), intent(in) :: iplcm
   integer, intent(in) :: leng,nlin
   character(len=*), intent(in) :: name
   character(len=*), dimension(nlin), intent(out) :: carr
   !----
   !  Local variables
   !----
   character(len=13)  :: fmt
   character(len=12) :: text12
   integer,allocatable :: ibase(:)
   !
   fmt='(           )'
   allocate(ibase(nlin*(leng+3)/4))
   !----
   !  Read from LCM object
   !----
   call LCMLEN(iplcm,name,ilong,itylcm)
   if(ilong == 0) then
      call LCMLIB(iplcm)
      text12=name
      call XABORT('LCMGTC: record '//text12//' not found.')
   endif
   call LCMGET(iplcm,name,ibase)
   !----
   !  Define format and dimensions for conversion
   !----
   n=leng/4
   m=mod(leng,4)
   if(n > 0)write(fmt(2:9),'(i6,''a4'')') n
   if(m > 0)then
      if(n.gt.0)fmt(10:10)=','
      write(fmt(11:12),'(''a'',i1)') m
      n=n+1
   endif
   !----
   !  Convert integers to strings
   !----
   do i=1,nlin
      write(carr(i)(1:leng),fmt) (ibase((i-1)*n+j),j=1,n)
   enddo
   deallocate(ibase)
end subroutine LCMGTC_S1
!
subroutine LCMGTC_S2(iplcm,name,leng,nlin,carr)
   use, intrinsic :: iso_c_binding
   use LCMAUX
   !----
   !  Subroutine arguments
   !----
   type(c_ptr), intent(in) :: iplcm
   integer, intent(in) :: leng,nlin
   character(len=*), intent(in) :: name
   character(len=*), dimension(:,:), intent(out) :: carr
   !----
   !  Local variables
   !----
   character(len=13)  :: fmt
   character(len=12) :: text12
   character(len=131) :: hsmg
   integer,allocatable :: ibase(:)
   !
   fmt='(           )'
   allocate(ibase(nlin*(leng+3)/4))
   !----
   !  Read from LCM object
   !----
   call LCMLEN(iplcm,name,ilong,itylcm)
   if(ilong == 0) then
      call LCMLIB(iplcm)
      text12=name
      call XABORT('LCMGTC: record '//text12//' not found.')
   endif
   call LCMGET(iplcm,name,ibase)
   !----
   !  Define format and dimensions for conversion
   !----
   n=leng/4
   m=mod(leng,4)
   if(n > 0)write(fmt(2:9),'(i6,''a4'')') n
   if(m > 0)then
      if(n.gt.0)fmt(10:10)=','
      write(fmt(11:12),'(''a'',i1)') m
      n=n+1
   endif
   !----
   !  Convert integers to strings
   !----
   nlin1=size(carr,1)
   nlin2=size(carr,2)
   if(nlin1*nlin2.ne.nlin) then
      write(hsmg,'(29hLCMGTC_S2: allocated length (,i5,17h) is not equal to, &
      & 16h argument size (,i5,2h).)') nlin1*nlin2,nlin
      call xabort(hsmg)
   endif
   do i2=1,nlin2
      do i1=1,nlin1
         i=(i2-1)*nlin1+i1
         write(carr(i1,i2)(1:leng),fmt) (ibase((i-1)*n+j),j=1,n)
      enddo
   enddo
   deallocate(ibase)
end subroutine LCMGTC_S2
!
subroutine LCMPTC_S0(iplcm,name,leng,carr)
   use, intrinsic :: iso_c_binding
   use LCMMOD
   !----
   !  Subroutine arguments
   !----
   type(c_ptr), intent(in) :: iplcm
   integer, intent(in) :: leng
   character(len=*), intent(in) :: name
   character(len=*), intent(in) :: carr
   !----
   !  Local variables
   !----
   character(len=13) :: fmt
   integer,allocatable :: ibase(:)
   !
   fmt='(           )'
   !----
   !  Define format and dimensions for conversion
   !----
   allocate(ibase((leng+3)/4))
   n=leng/4
   m=mod(leng,4)
   if(n > 0)write(fmt(2:9),'(i6,''a4'')') n
   if(m > 0)then
      if(n.gt.0)fmt(10:10)=','
      write(fmt(11:12),'(''a'',i1)') m
      n=n+1
   endif
   !----
   !  Convert strings to integers
   !----
   read(carr(1:leng),fmt) (ibase(j),j=1,n)
   !----
   !  Write to LCM object
   !----
   call LCMPUT(iplcm,name,n,3,ibase)
   deallocate(ibase)
end subroutine LCMPTC_S0
!
subroutine LCMPTC_S1(iplcm,name,leng,nlin,carr)
   use, intrinsic :: iso_c_binding
   use LCMMOD
   !----
   !  Subroutine arguments
   !----
   type(c_ptr), intent(in) :: iplcm
   integer, intent(in) :: leng,nlin
   character(len=*), intent(in) :: name
   character(len=*), dimension(nlin), intent(in) :: carr
   !----
   !  Local variables
   !----
   character(len=13) :: fmt
   integer,allocatable :: ibase(:)
   !
   fmt='(           )'
   !----
   !  Define format and dimensions for conversion
   !----
   allocate(ibase(nlin*(leng+3)/4))
   n=leng/4
   m=mod(leng,4)
   if(n > 0)write(fmt(2:9),'(i6,''a4'')') n
   if(m > 0)then
      if(n.gt.0)fmt(10:10)=','
      write(fmt(11:12),'(''a'',i1)') m
      n=n+1
   endif
   !----
   !  Convert strings to integers
   !----
   do i=1,nlin
      read(carr(i)(1:leng),fmt) (ibase((i-1)*n+j),j=1,n)
   enddo
   !----
   !  Write to LCM object
   !----
   call LCMPUT(iplcm,name,n*nlin,3,ibase)
   deallocate(ibase)
end subroutine LCMPTC_S1
!
subroutine LCMPTC_S2(iplcm,name,leng,nlin,carr)
   use, intrinsic :: iso_c_binding
   use LCMMOD
   !----
   !  Subroutine arguments
   !----
   type(c_ptr), intent(in) :: iplcm
   integer, intent(in) :: leng,nlin
   character(len=*), intent(in) :: name
   character(len=*), dimension(:,:), intent(in) :: carr
   !----
   !  Local variables
   !----
   character(len=13) :: fmt
   character(len=131) :: hsmg
   integer,allocatable :: ibase(:)
   !
   fmt='(           )'
   !----
   !  Define format and dimensions for conversion
   !----
   allocate(ibase(nlin*(leng+3)/4))
   n=leng/4
   m=mod(leng,4)
   if(n > 0)write(fmt(2:9),'(i6,''a4'')') n
   if(m > 0)then
      if(n.gt.0)fmt(10:10)=','
      write(fmt(11:12),'(''a'',i1)') m
      n=n+1
   endif
   !----
   !  Convert strings to integers
   !----
   nlin1=size(carr,1)
   nlin2=size(carr,2)
   if(nlin1*nlin2.ne.nlin) then
      write(hsmg,'(29hLCMPTC_S2: allocated length (,i5,17h) is not equal to, &
      & 16h argument size (,i5,2h).)') nlin1*nlin2,nlin
      call xabort(hsmg)
   endif
   do i2=1,nlin2
      do i1=1,nlin1
         i=(i2-1)*nlin1+i1
         read(carr(i1,i2)(1:leng),fmt) (ibase((i-1)*n+j),j=1,n)
      enddo
   enddo
   !----
   !  Write to LCM object
   !----
   call LCMPUT(iplcm,name,n*nlin,3,ibase)
   deallocate(ibase)
end subroutine LCMPTC_S2
!
subroutine LCMGLC_S0(iplcm,ipos,leng,carr)
   use, intrinsic :: iso_c_binding
   use LCMAUX
   !----
   !  Subroutine arguments
   !----
   type(c_ptr), intent(in) :: iplcm
   integer, intent(in) :: ipos,leng
   character(len=*), intent(out) :: carr
   !----
   !  Local variables
   !----
   character(len=13) :: fmt
   character(len=131) :: hsmg
   integer,allocatable :: ibase(:)
   !
   fmt='(           )'
   allocate(ibase((leng+3)/4))
   !----
   !  Read from LCM object
   !----
   call lcmlel(iplcm,ipos,ilong,itylcm)
   if(ilong == 0) then
      call LCMLIB(iplcm)
      write(hsmg,'(8hLCMGLC: ,i5,21h-th record not found.)') ipos
      call XABORT(hsmg)
   endif
   call LCMGDL(iplcm,ipos,ibase)
   !----
   !  Define format and dimensions for conversion
   !----
   n=leng/4
   m=mod(leng,4)
   if(n > 0)write(fmt(2:9),'(i6,''a4'')') n
   if(m > 0)then
      if(n.gt.0)fmt(10:10)=','
      write(fmt(11:12),'(''a'',i1)') m
      n=n+1
   endif
   !----
   !  Convert integer to strings
   !----
   write(carr(1:leng),fmt) (ibase(j),j=1,n)
   deallocate(ibase)
end subroutine LCMGLC_S0
!
subroutine LCMGLC_S1(iplcm,ipos,leng,nlin,carr)
   use, intrinsic :: iso_c_binding
   use LCMAUX
   !----
   !  Subroutine arguments
   !----
   type(c_ptr), intent(in) :: iplcm
   integer, intent(in) :: ipos,leng,nlin
   character(len=*), dimension(nlin), intent(out) :: carr
   !----
   !  Local variables
   !----
   character(len=13) :: fmt
   character(len=131) :: hsmg
   integer,allocatable :: ibase(:)
   !
   fmt='(           )'
   allocate(ibase(nlin*(leng+3)/4))
   !----
   !  Read from LCM object
   !----
   call lcmlel(iplcm,ipos,ilong,itylcm)
   if(ilong == 0) then
      call LCMLIB(iplcm)
      write(hsmg,'(8hLCMGLC: ,i5,21h-th record not found.)') ipos
      call XABORT(hsmg)
   endif
   call LCMGDL(iplcm,ipos,ibase)
   !----
   !  Define format and dimensions for conversion
   !----
   n=leng/4
   m=mod(leng,4)
   if(n > 0)write(fmt(2:9),'(i6,''a4'')') n
   if(m > 0)then
      if(n.gt.0)fmt(10:10)=','
      write(fmt(11:12),'(''a'',i1)') m
      n=n+1
   endif
   !----
   !  Convert integers to strings
   !----
   do i=1,nlin
      write(carr(i)(1:leng),fmt) (ibase((i-1)*n+j),j=1,n)
   enddo
   deallocate(ibase)
end subroutine LCMGLC_S1
!
subroutine LCMGLC_S2(iplcm,ipos,leng,nlin,carr)
   use, intrinsic :: iso_c_binding
   use LCMAUX
   !----
   !  Subroutine arguments
   !----
   type(c_ptr), intent(in) :: iplcm
   integer, intent(in) :: ipos,leng,nlin
   character(len=*), dimension(:,:), intent(out) :: carr
   !----
   !  Local variables
   !----
   character(len=13) :: fmt
   character(len=131) :: hsmg
   integer,allocatable :: ibase(:)
   !
   fmt='(           )'
   allocate(ibase(nlin*(leng+3)/4))
   !----
   !  Read from LCM object
   !----
   call lcmlel(iplcm,ipos,ilong,itylcm)
   if(ilong == 0) then
      call LCMLIB(iplcm)
      write(hsmg,'(8hLCMGLC: ,i5,21h-th record not found.)') ipos
      call XABORT(hsmg)
   endif
   call LCMGDL(iplcm,ipos,ibase)
   !----
   !  Define format and dimensions for conversion
   !----
   n=leng/4
   m=mod(leng,4)
   if(n > 0)write(fmt(2:9),'(i6,''a4'')') n
   if(m > 0)then
      if(n.gt.0)fmt(10:10)=','
      write(fmt(11:12),'(''a'',i1)') m
      n=n+1
   endif
   !----
   !  Convert integers to strings
   !----
   nlin1=size(carr,1)
   nlin2=size(carr,2)
   if(nlin1*nlin2.ne.nlin) then
      write(hsmg,'(29hLCMGLC_S2: allocated length (,i5,17h) is not equal to, &
      & 16h argument size (,i5,2h).)') nlin1*nlin2,nlin
      call xabort(hsmg)
   endif
   do i2=1,nlin2
      do i1=1,nlin1
         i=(i2-1)*nlin1+i1
         write(carr(i1,i2)(1:leng),fmt) (ibase((i-1)*n+j),j=1,n)
      enddo
   enddo
   deallocate(ibase)
end subroutine LCMGLC_S2
!
subroutine LCMPLC_S0(iplcm,ipos,leng,carr)
   use, intrinsic :: iso_c_binding
   use LCMMOD
   !----
   !  Subroutine arguments
   !----
   type(c_ptr), intent(in) :: iplcm
   integer, intent(in) :: ipos,leng
   character(len=*), intent(in) :: carr
   !----
   !  Local variables
   !----
   character(len=13) :: fmt
   integer,allocatable :: ibase(:)
   !
   fmt='(           )'
   !----
   !  Define format and dimensions for conversion
   !----
   allocate(ibase((leng+3)/4))
   n=leng/4
   m=mod(leng,4)
   if(n > 0)write(fmt(2:9),'(i6,''a4'')') n
   if(m > 0)then
      if(n.gt.0)fmt(10:10)=','
      write(fmt(11:12),'(''a'',i1)') m
      n=n+1
   endif
   !----
   !  Convert strings to integers
   !----
   read(carr(1:leng),fmt) (ibase(j),j=1,n)
   !----
   !  Write to LCM object
   !----
   call LCMPDL(iplcm,ipos,n,3,ibase)
   deallocate(ibase)
end subroutine LCMPLC_S0
!
subroutine LCMPLC_S1(iplcm,ipos,leng,nlin,carr)
   use, intrinsic :: iso_c_binding
   use LCMMOD
   !----
   !  Subroutine arguments
   !----
   type(c_ptr), intent(in) :: iplcm
   integer, intent(in) :: ipos,leng,nlin
   character(len=*), dimension(nlin), intent(in) :: carr
   !----
   !  Local variables
   !----
   character(len=13) :: fmt
   integer,allocatable :: ibase(:)
   !
   fmt='(           )'
   !----
   !  Define format and dimensions for conversion
   !----
   allocate(ibase(nlin*(leng+3)/4))
   n=leng/4
   m=mod(leng,4)
   if(n > 0)write(fmt(2:9),'(i6,''a4'')') n
   if(m > 0)then
      if(n.gt.0)fmt(10:10)=','
      write(fmt(11:12),'(''a'',i1)') m
      n=n+1
   endif
   !----
   !  Convert strings to integers
   !----
   do i=1,nlin
      read(carr(i)(1:leng),fmt) (ibase((i-1)*n+j),j=1,n)
   enddo
   !----
   !  Write to LCM object
   !----
   call LCMPDL(iplcm,ipos,n*nlin,3,ibase)
   deallocate(ibase)
end subroutine LCMPLC_S1
!
subroutine LCMPLC_S2(iplcm,ipos,leng,nlin,carr)
   use, intrinsic :: iso_c_binding
   use LCMMOD
   !----
   !  Subroutine arguments
   !----
   type(c_ptr), intent(in) :: iplcm
   integer, intent(in) :: ipos,leng,nlin
   character(len=*), dimension(:,:), intent(in) :: carr
   !----
   !  Local variables
   !----
   character(len=13) :: fmt
   character(len=131) :: hsmg
   integer,allocatable :: ibase(:)
   !
   fmt='(           )'
   !----
   !  Define format and dimensions for conversion
   !----
   allocate(ibase(nlin*(leng+3)/4))
   n=leng/4
   m=mod(leng,4)
   if(n > 0)write(fmt(2:9),'(i6,''a4'')') n
   if(m > 0)then
      if(n.gt.0)fmt(10:10)=','
      write(fmt(11:12),'(''a'',i1)') m
      n=n+1
   endif
   !----
   !  Convert strings to integers
   !----
   nlin1=size(carr,1)
   nlin2=size(carr,2)
   if(nlin1*nlin2.ne.nlin) then
      write(hsmg,'(29hLCMPLC_S2: allocated length (,i5,17h) is not equal to, &
      & 16h argument size (,i5,2h).)') nlin1*nlin2,nlin
      call xabort(hsmg)
   endif
   do i2=1,nlin2
      do i1=1,nlin1
         i=(i2-1)*nlin1+i1
         read(carr(i1,i2)(1:leng),fmt) (ibase((i-1)*n+j),j=1,n)
      enddo
   enddo
   !----
   !  Write to LCM object
   !----
   call LCMPDL(iplcm,ipos,n*nlin,3,ibase)
   deallocate(ibase)
end subroutine LCMPLC_S2
end module LCMTLC
