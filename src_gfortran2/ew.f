      REAL FUNCTION EW(X)
C     ****************************************************************
C     SAETTIGUNGSDAMPFDRUCK UEBER WASSER IN PA. X IN KELVIN.
C     NACH DER GOFF-GRATCH-FORMEL.
C     ****************************************************************
      EW=0.
      IF(X.LE.0.) STOP 'SORRY: T NOT IN [K]'
      Y=373.16/X
      A=-7.90298*(Y-1.)
      A=A+(5.02808*0.43429*ALOG(Y))
      C=(1.-(1./Y))*11.344
      C=-1.+(10.**C)
      C=-1.3816*C/(10.**7)
      D=(1.-Y)*3.49149
      D=-1.+(10.**D)
      D=8.1328*D/(10.**3)
      Y=A+C+D
      EW=101324.6*(10.**Y)       ! Saettigungsdampfdruck in Pa

      END
