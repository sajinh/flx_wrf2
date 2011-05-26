      SUBROUTINE CALDATE(JULDATE,YYYYMMDD,HHMISS)
c                           i       o       o
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
*                                                                             *
*     Calculates the Gregorian date from the Julian date                      *
*                                                                             *
*     AUTHOR: Andreas Stohl (21 January 1994), adapted from Numerical Recipes *
*                                                                             *
*     Variables:                                                              *
*     DD             Day                                                      *
*     HH             Hour                                                     *
*     HHMISS         Hour, Minute, Second                                     *
*     JA,JB,JC,JD,JE help variables                                           *
*     JALPHA         help variable                                            *
*     JULDATE        Julian Date                                              *
*     JULDAY         help variable                                            *
*     MI             Minute                                                   *
*     MM             Month                                                    *
*     SS             Seconds                                                  *
*     YYYY           Year                                                     *
*     YYYYMMDD       Year, Month, Day                                         *
*                                                                             *
*     Constants:                                                              *
*     IGREG          help constant                                            *
*                                                                             *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

      IMPLICIT NONE

      INTEGER YYYYMMDD,YYYY,MM,DD,HHMISS,HH,MI,SS
      INTEGER JULDAY,JA,JB,JC,JD,JE,IGREG,JALPHA
      DOUBLE PRECISION JULDATE
      PARAMETER (IGREG=2299161)

      JULDAY=INT(JULDATE)
      IF(JULDAY.GE.IGREG)THEN
        JALPHA=INT(((JULDAY-1867216)-0.25)/36524.25)
        JA=JULDAY+1+JALPHA-INT(0.25*JALPHA)
      ELSE
        JA=JULDAY
      ENDIF
      JB=JA+1524
      JC=INT(6680.+((JB-2439870)-122.1)/365.25)
      JD=365*JC+INT(0.25*JC)
      JE=INT((JB-JD)/30.6001)
      DD=JB-JD-INT(30.6001*JE)
      MM=JE-1
      IF (MM.GT.12) MM=MM-12
      YYYY=JC-4715
      IF (MM.GT.2) YYYY=YYYY-1
      IF (YYYY.LE.0) YYYY=YYYY-1

      YYYYMMDD=10000*YYYY+100*MM+DD
      HH=INT(24.*(JULDATE-FLOAT(JULDAY)))
      MI=INT(1440.*(JULDATE-FLOAT(JULDAY))-60.*FLOAT(HH))
      SS=NINT(86400.*(JULDATE-FLOAT(JULDAY))-3600.*FLOAT(HH))
     +-60.*FLOAT(MI)
      IF (SS.EQ.60) THEN  ! 60 seconds = 1 minute
        SS=0
        MI=MI+1
      ENDIF
      IF (MI.EQ.60) THEN
        MI=0
        HH=HH+1
      ENDIF
      HHMISS=10000*HH+100*MI+SS

      END
