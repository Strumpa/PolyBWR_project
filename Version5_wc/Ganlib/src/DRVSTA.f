*DECK DRVSTA
      SUBROUTINE DRVSTA(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
* STANDARD COMPARE MODULE.
*
* INPUT/OUTPUT PARAMETERS:
*  NENTRY : NUMBER OF LINKED LISTS AND FILES USED BY THE MODULE.
*  HENTRY : CHARACTER*12 NAME OF EACH LINKED LIST OR FILE.
*  IENTRY : =0 CLE-2000 VARIABLE; =1 LINKED LIST; =2 XSM FILE;
*           =3 SEQUENTIAL BINARY FILE; =4 SEQUENTIAL ASCII FILE;
*           =5 DIRECT ACCESS FILE.
*  JENTRY : =0 THE LINKED LIST OR FILE IS CREATED.
*           =1 THE LINKED LIST OR FILE IS OPEN FOR MODIFICATIONS;
*           =2 THE LINKED LIST OR FILE IS OPEN IN READ-ONLY MODE.
*  KENTRY : =FILE UNIT NUMBER; =LINKED LIST STARESS OTHERWISE.
*           DIMENSION HENTRY(NENTRY),IENTRY(NENTRY),JENTRY(NENTRY),
*           KENTRY(NENTRY)
*
*-------------------------------------- AUTHOR: A. HEBERT ; 21/12/93 ---
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR) KENTRY(NENTRY)
      CHARACTER HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      CHARACTER TEXT12*12
      TYPE(C_PTR) IPLIST1,IPLIST2
*----
*  PARAMETER VALIDATION.
*----
      IF(NENTRY.LE.1) CALL XABORT('DRVSTA: TWO PARAMETER EXPECTED.')
      TEXT12=HENTRY(1)
      IF((JENTRY(1).NE.2).OR.(IENTRY(1).GT.2)) CALL XABORT('DRVSTA: LIN'
     1 //'KED LIST OR XSM FILE IN READ-ONLY MODE EXPECTED AT RHS ('
     2 //TEXT12//').')
      IF((JENTRY(2).NE.2).OR.(IENTRY(2).GT.2)) CALL XABORT('DRVSTA: LIN'
     1 //'KED LIST OR XSM FILE IN READ-ONLY MODE EXPECTED AT RHS ('
     2 //TEXT12//').')
*----
*  PERFORM THE COMPARISON.
*----
      IPLIST1=KENTRY(1)
      IPLIST2=KENTRY(2)
      CALL LCMSTA(IPLIST2,IPLIST1)
      RETURN
      END