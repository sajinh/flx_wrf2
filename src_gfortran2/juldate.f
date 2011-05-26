      FUNCTION JULDATE(YYYYMMDD,HHMISS)
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
*                                                                             *
*     Calculates the Julian date                                              *
*                                                                             *
*     AUTHOR: Andreas Stohl (15 October 1993)                                 *
*                                                                             *
*     Variables:                                                              *
*     DD             Day                                                      *
*     HH             Hour                                                     *
*     HHMISS         Hour, minute + second                                    *
*     JA,JM,JY       help variables                                           *
*     JULDATE        Julian Date                                              *
*     JULDAY         help variable                                            *
*     MI             Minute                                                   *
*     MM             Month                                                    *
*     SS             Second                                                   *
*     YYYY           Year                                                     *
*     YYYYMMDDHH     Date and Time                                            *
*                                                                             *
*     Constants:                                                              *
*     IGREG          help constant                                            *
*                                                                             *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

      IMPLICIT NONE

      INTEGER YYYYMMDD,YYYY,MM,DD,HH,MI,SS,HHMISS
      INTEGER JULDAY,JY,JM,JA,IGREG
      DOUBLE PRECISION JULDATE
      PARAMETER (IGREG=15+31*(10+12*1582))

      YYYY=YYYYMMDD/10000
      MM=(YYYYMMDD-10000*YYYY)/100
      DD=YYYYMMDD-10000*YYYY-100*MM
      HH=HHMISS/10000
      MI=(HHMISS-10000*HH)/100
      SS=HHMISS-10000*HH-100*MI

c------modified 2011/03/29
      IF (YYYY.EQ.0) then
         write(0,*) 'There is no Year Zero.'
         return
      endif
c------------------------
      IF (YYYY.LT.0) YYYY=YYYY+1
      IF (MM.GT.2) THEN
        JY=YYYY
        JM=MM+1
      ELSE
        JY=YYYY-1
        JM=MM+13
      ENDIF
      JULDAY=INT(365.25*JY)+INT(30.6001*JM)+DD+1720995
      IF (DD+31*(MM+12*YYYY).GE.IGREG) THEN
        JA=INT(0.01*JY)
        JULDAY=JULDAY+2-JA+INT(0.25*JA)
      ENDIF

      JULDATE=DBLE(FLOAT(JULDAY))+DBLE(FLOAT(HH)/24.)+
     +DBLE(FLOAT(MI)/1440.)+DBLE(FLOAT(SS)/86400.)

      END
