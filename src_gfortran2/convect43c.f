***************************************************************************
*****                       SUBROUTINE CONVECT                        *****
*****                          VERSION 4.3c                           *****
*****                          20 May, 2002                           *****
*****                          Kerry Emanuel                          *****
***************************************************************************
C
        SUBROUTINE CONVECT
     *    (ND,  NL,   DELT, IFLAG,  
     *     PRECIP, WD,   TPRIME, QPRIME, CBMF    )
c
c-cv ***************************************************************************
c-cv C. Forster, November 2003 - May 2004:
c-cv
c-cv The subroutine has been downloaded from Kerry Emanuel's homepage, 
c-cv where further infos on the convection scheme can be found
c-cv http://www-paoc.mit.edu/~emanuel/home.html
c-cv
c-cv The following changes have been made to integrate this subroutine
c-cv into FLEXPART
c-cv
c-cv Putting most of the variables in a new common block
c-cv renaming eps to eps0 because there is some eps already in includepar
c-cv
c-cv removing the arrays U,V,TRA and related arrays
c-cv
c-cv renaming the original arrays T,Q,QS,P,PH to
c-cv TCONV,QCONV,QSCONV,PCONV_HPA,PHCONV_HPA
c-cv
c-cv Initialization of variables has been put into parameter statements instead
c-cv of assignment of values at each call, in order to save computation time. 
c***************************************************************************
C
C-----------------------------------------------------------------------------
C    *** On input:      ***
C
C     T:   Array of absolute temperature (K) of dimension ND, with first
C           index corresponding to lowest model level. Note that this array
C           will be altered by the subroutine if dry convective adjustment
C           occurs and if IPBL is not equal to 0.
C
C     Q:   Array of specific humidity (gm/gm) of dimension ND, with first
C            index corresponding to lowest model level. Must be defined
C            at same grid levels as T. Note that this array will be altered
C            if dry convective adjustment occurs and if IPBL is not equal to 0.
C
C     QS:  Array of saturation specific humidity of dimension ND, with first
C            index corresponding to lowest model level. Must be defined
C            at same grid levels as T. Note that this array will be altered
C            if dry convective adjustment occurs and if IPBL is not equal to 0.
C
C     U:   Array of zonal wind velocity (m/s) of dimension ND, witth first
C            index corresponding with the lowest model level. Defined at
C            same levels as T. Note that this array will be altered if
C            dry convective adjustment occurs and if IPBL is not equal to 0.
C
C     V:   Same as U but for meridional velocity.
C
C     TRA: Array of passive tracer mixing ratio, of dimensions (ND,NTRA),
C            where NTRA is the number of different tracers. If no
C            convective tracer transport is needed, define a dummy
C            input array of dimension (ND,1). Tracers are defined at
C            same vertical levels as T. Note that this array will be altered
C            if dry convective adjustment occurs and if IPBL is not equal to 0.
C
C     P:   Array of pressure (mb) of dimension ND, with first
C            index corresponding to lowest model level. Must be defined
C            at same grid levels as T.
C
C     PH:  Array of pressure (mb) of dimension ND+1, with first index
C            corresponding to lowest level. These pressures are defined at
C            levels intermediate between those of P, T, Q and QS. The first
C            value of PH should be greater than (i.e. at a lower level than)
C            the first value of the array P.
C
C     ND:  The dimension of the arrays T,Q,QS,P,PH,FT and FQ
C
C     NL:  The maximum number of levels to which convection can
C            penetrate, plus 1.
C            NL MUST be less than or equal to ND-1.
C
C     NTRA:The number of different tracers. If no tracer transport
C            is needed, set this equal to 1. (On most compilers, setting
C            NTRA to 0 will bypass tracer calculation, saving some CPU.)  
C
C     DELT: The model time step (sec) between calls to CONVECT
C
C----------------------------------------------------------------------------
C    ***   On Output:         ***
C
C     IFLAG: An output integer whose value denotes the following:
C
C                VALUE                        INTERPRETATION
C                -----                        --------------
C                  0               No moist convection; atmosphere is not
C                                  unstable, or surface temperature is less
C                                  than 250 K or surface specific humidity
C                                  is non-positive.
C
C                  1               Moist convection occurs.
C
C                  2               No moist convection: lifted condensation
C                                  level is above the 200 mb level.
C
C                  3               No moist convection: cloud base is higher
C                                  then the level NL-1.
C
C                  4               Moist convection occurs, but a CFL condition
C                                  on the subsidence warming is violated. This
C                                  does not cause the scheme to terminate.
C
C     FT:   Array of temperature tendency (K/s) of dimension ND, defined at same
C             grid levels as T, Q, QS and P.
C
C     FQ:   Array of specific humidity tendencies ((gm/gm)/s) of dimension ND,
C             defined at same grid levels as T, Q, QS and P.
C
C     FU:   Array of forcing of zonal velocity (m/s^2) of dimension ND,
C             defined at same grid levels as T.
C
C     FV:   Same as FU, but for forcing of meridional velocity.
C
C     FTRA: Array of forcing of tracer content, in tracer mixing ratio per
C             second, defined at same levels as T. Dimensioned (ND,NTRA).
C
C     PRECIP: Scalar convective precipitation rate (mm/day).
C
C     WD:    A convective downdraft velocity scale. For use in surface
C             flux parameterizations. See convect.ps file for details.
C
C     TPRIME: A convective downdraft temperature perturbation scale (K).
C              For use in surface flux parameterizations. See convect.ps
C              file for details.
C
C     QPRIME: A convective downdraft specific humidity
C              perturbation scale (gm/gm).
C              For use in surface flux parameterizations. See convect.ps
C              file for details.
C
C     CBMF:   The cloud base mass flux ((kg/m**2)/s). THIS SCALAR VALUE MUST
C              BE STORED BY THE CALLING PROGRAM AND RETURNED TO CONVECT AT
C              ITS NEXT CALL. That is, the value of CBMF must be "remembered"
C              by the calling program between calls to CONVECT.
C
C------------------------------------------------------------------------------
C
C    ***  THE PARAMETER NA SHOULD IN GENERAL BE GREATER THAN   ***
C    ***                OR EQUAL TO  ND + 1                    ***
C
c
      include 'includepar'
      include 'includeconv'
C
C-cv====>Begin Module CONVECT    File convect.f      Undeclared variables
C
C     Argument variables
C
      integer iflag, nd, nl
C
      real cbmf, delt, precip, qprime, tprime, wd
C
C     Local variables
C
      integer i, icb, ihmin, inb, inb1, ipbl, j, jtt, k, minorig
      integer nk
C
      real ad, afac, ahmax, ahmin, alpha, alt, altem
      real am, amp1, anum, asij, awat, b6, beta, bf2, bsum, by
      real byp, c6, cape, capem, cbmfold, chi, cl, coeff, coeffr, coeffs
      real cpd, cpinv, cpv, cpvmcl, cwat, damp, damps, dbo, dbosum
      real defrac, dei, delm, delp, delt0, delti, denom, dhdp
      real dpinv, dtma, dtmax, dtmin, dtpbl, elacrit, elcrit, entp, ents
      real epmax, eps0, epsi, fac, fqold, frac, ftold, g, ginv
      real omtrain, omtsnow, plcl, qp1, qsm, qstm, qti, rat, rd
      real rdcp, revap, rh, rowl, rv, scrit, sigd, sigs, sigt, sjmax
      real sjmin, smid, smin, stemp, tca, tlcrit
      real tvaplcl, tvpplcl, tvx, tvy, wdtrain
      real cu

c     integer jc,jn
c     real alvnew,a2,ahm,alv,rm,sum,qnew,dphinv,tc,thbar,tnew,x

      real FUP(NA),FDOWN(NA)
C
C-cv====>End Module   CONVECT    File convect.f

      INTEGER NENT(NA)
      REAL M(NA),MP(NA),MENT(NA,NA),QENT(NA,NA),ELIJ(NA,NA)
      REAL SIJ(NA,NA),TVP(NA),TV(NA),WATER(NA)
      REAL QP(NA),EP(NA),TH(NA),WT(NA),EVAP(NA),CLW(NA)
      REAL SIGP(NA),TP(NA),CPN(NA)
      REAL LV(NA),LVCP(NA),LV0,H(NA),HP(NA),GZ(NA),HM(NA)
      REAL EPSILON
C     REAL TOLD(NA)
C
C -----------------------------------------------------------------------
C
C   ***                     Specify Switches                         ***
C
C   ***   IPBL: Set to zero to bypass dry adiabatic adjustment       ***
C   ***    Any other value results in dry adiabatic adjustment       ***
C   ***     (Zero value recommended for use in models with           ***
C   ***                   boundary layer schemes)                    ***
C
C   ***   MINORIG: Lowest level from which convection may originate  ***
C   ***     (Should be first model level at which T is defined       ***
C   ***      for models using bulk PBL schemes; otherwise, it should ***
C   ***      be the first model level at which T is defined above    ***
C   ***                      the surface layer)                      ***
C
        PARAMETER(IPBL=0)
        PARAMETER(MINORIG=1)
C
C------------------------------------------------------------------------------
C
C   ***                    SPECIFY PARAMETERS                        ***
C
C   *** ELCRIT IS THE AUTOCONVERSION THERSHOLD WATER CONTENT (gm/gm) ***
C   ***  TLCRIT IS CRITICAL TEMPERATURE BELOW WHICH THE AUTO-        ***
C   ***       CONVERSION THRESHOLD IS ASSUMED TO BE ZERO             ***
C   ***     (THE AUTOCONVERSION THRESHOLD VARIES LINEARLY            ***
C   ***               BETWEEN 0 C AND TLCRIT)                        ***
C   ***   ENTP IS THE COEFFICIENT OF MIXING IN THE ENTRAINMENT       ***
C   ***                       FORMULATION                            ***
C   ***  SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT  ***
C   ***  SIGS IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE       ***
C   ***                        OF CLOUD                              ***
C   ***        OMTRAIN IS THE ASSUMED FALL SPEED (P/s) OF RAIN       ***
C   ***     OMTSNOW IS THE ASSUMED FALL SPEED (P/s) OF SNOW          ***
C   ***  COEFFR IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
C   ***                          OF RAIN                             ***
C   ***  COEFFS IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
C   ***                          OF SNOW                             ***
C   ***     CU IS THE COEFFICIENT GOVERNING CONVECTIVE MOMENTUM      ***
C   ***                         TRANSPORT                            ***
C   ***    DTMAX IS THE MAXIMUM NEGATIVE TEMPERATURE PERTURBATION    ***
C   ***        A LIFTED PARCEL IS ALLOWED TO HAVE BELOW ITS LFC      ***
C   ***    ALPHA AND DAMP ARE PARAMETERS THAT CONTROL THE RATE OF    ***
C   ***                 APPROACH TO QUASI-EQUILIBRIUM                ***
C   ***   (THEIR STANDARD VALUES ARE  0.20 AND 0.1, RESPECTIVELY)    ***
C   ***                   (DAMP MUST BE LESS THAN 1)                 ***
C
        PARAMETER(ELCRIT=.0011)
        PARAMETER(TLCRIT=-55.0)
        PARAMETER(ENTP=1.5)
        PARAMETER(SIGD=0.05)
        PARAMETER(SIGS=0.12)
        PARAMETER(OMTRAIN=50.0)
        PARAMETER(OMTSNOW=5.5) 
        PARAMETER(COEFFR=1.0)
        PARAMETER(COEFFS=0.8)
        PARAMETER(CU=0.7)
        PARAMETER(BETA=10.0)
        PARAMETER(DTMAX=0.9) 
        PARAMETER(ALPHA=0.025)  !original 0.2
        PARAMETER(DAMP=0.1)
C
C   ***        ASSIGN VALUES OF THERMODYNAMIC CONSTANTS,        ***
C   ***            GRAVITY, AND LIQUID WATER DENSITY.           ***
C   ***             THESE SHOULD BE CONSISTENT WITH             ***
C   ***              THOSE USED IN CALLING PROGRAM              ***
C   ***     NOTE: THESE ARE ALSO SPECIFIED IN SUBROUTINE TLIFT  ***
C
      PARAMETER(CPD=1005.7)
      PARAMETER(CPV=1870.0)
      PARAMETER(CL=2500.0) 
      PARAMETER(RV=461.5)
      PARAMETER(RD=287.04)
      PARAMETER(LV0=2.501E6)
      PARAMETER(G=9.81)  
      PARAMETER(ROWL=1000.0)
C
      PARAMETER(CPVMCL=CL-CPV) 
      PARAMETER(EPS0=RD/RV)
      PARAMETER(EPSI=1./EPS0)
      PARAMETER(GINV=1.0/G)
      PARAMETER(EPSILON=1.e-20)

C EPSILON IS A SMALL NUMBER USED TO EXCLUDE MASS FLUXES OF ZERO
c
      DELTI=1.0/DELT
C
C           ***  INITIALIZE OUTPUT ARRAYS AND PARAMETERS  ***
C

        DO 5 I=1,NL+1
         FT(I)=0.0
         FQ(I)=0.0
         FDOWN(I)=0.0
         SUB(I)=0.0
         FUP(I)=0.0
         M(I)=0.0
         MP(I)=0.0
        DO 5 J=1,NL+1
         FMASS(I,J)=0.0
         MENT(I,J)=0.0
    5   CONTINUE
        DO 7 I=1,NL+1
         RDCP=(RD*(1.-QCONV(I))+QCONV(I)*RV)/
     1    (CPD*(1.-QCONV(I))+QCONV(I)*CPV)
         TH(I)=TCONV(I)*(1000.0/PCONV_HPA(I))**RDCP
    7   CONTINUE
        PRECIP=0.0
        WD=0.0
        TPRIME=0.0
        QPRIME=0.0
        IFLAG=0
C
c       IF(IPBL.NE.0)THEN
C
C     ***            PERFORM DRY ADIABATIC ADJUSTMENT            ***
C
c       JC=0
c       DO 30 I=NL-1,1,-1
c        JN=0
c         SUM=TH(I)*(1.+QCONV(I)*EPSI-QCONV(I))
c        DO 10 J=I+1,NL
c         SUM=SUM+TH(J)*(1.+QCONV(J)*EPSI-QCONV(J))
c         THBAR=SUM/FLOAT(J+1-I)
c         IF((TH(J)*(1.+QCONV(J)*EPSI-QCONV(J))).LT.THBAR)JN=J
c  10    CONTINUE
c        IF(I.EQ.1)JN=MAX(JN,2)
c        IF(JN.EQ.0)GOTO 30
c  12    CONTINUE
c        AHM=0.0
c        RM=0.0
c        DO 15 J=I,JN
c         AHM=AHM+(CPD*(1.-QCONV(J))+QCONV(J)*CPV)*TCONV(J)*
c    +   (PHCONV_HPA(J)-PHCONV_HPA(J+1))
c         RM=RM+QCONV(J)*(PHCONV_HPA(J)-PHCONV_HPA(J+1))
c  15    CONTINUE
c        DPHINV=1./(PHCONV_HPA(I)-PHCONV_HPA(JN+1))
c        RM=RM*DPHINV
c        A2=0.0
c        DO 20 J=I,JN
c         QCONV(J)=RM
c         RDCP=(RD*(1.-QCONV(J))+QCONV(J)*RV)/
c    1     (CPD*(1.-QCONV(J))+QCONV(J)*CPV)  
c         X=(0.001*PCONV_HPA(J))**RDCP
c         TOLD(J)=TCONV(J)
c         TCONV(J)=X
c         A2=A2+(CPD*(1.-QCONV(J))+QCONV(J)*CPV)*X*
c    1    (PHCONV_HPA(J)-PHCONV_HPA(J+1))
c  20    CONTINUE
c        DO 25 J=I,JN
c         TH(J)=AHM/A2
c         TCONV(J)=TCONV(J)*TH(J)
c         TC=TOLD(J)-273.15
c         ALV=LV0-CPVMCL*TC
c         QSCONV(J)=QSCONV(J)+QSCONV(J)*(1.+QSCONV(J)*(EPSI-1.))*ALV*
c    1    (TCONV(J)- TOLD(J))/(RV*TOLD(J)*TOLD(J))
c      if (qslev(j) .lt. 0.) then
c        write(*,*) 'qslev.lt.0 ',j,qslev
c      endif
c  25    CONTINUE
c        IF((TH(JN+1)*(1.+QCONV(JN+1)*EPSI-QCONV(JN+1))).LT.
c    1    (TH(JN)*(1.+QCONV(JN)*EPSI-QCONV(JN))))THEN
c         JN=JN+1
c         GOTO 12
c        END IF
c        IF(I.EQ.1)JC=JN 
c  30   CONTINUE
C
C   ***   Remove any supersaturation that results from adjustment ***
C
c     IF(JC.GT.1)THEN
c      DO 38 J=1,JC
c         IF(QSCONV(J).LT.QCONV(J))THEN 
c          ALV=LV0-CPVMCL*(TCONV(J)-273.15)  
c          TNEW=TCONV(J)+ALV*(QCONV(J)-QSCONV(J))/(CPD*(1.-QCONV(J))+
c    1      CL*QCONV(J)+QSCONV(J)*(CPV-CL+ALV*ALV/(RV*TCONV(J)*TCONV(J))))
c          ALVNEW=LV0-CPVMCL*(TNEW-273.15)
c          QNEW=(ALV*QCONV(J)-(TNEW-TCONV(J))*(CPD*(1.-QCONV(J))
c    1     +CL*QCONV(J)))/ALVNEW
c          PRECIP=PRECIP+24.*3600.*1.0E5*(PHCONV_HPA(J)-PHCONV_HPA(J+1))*
c    1      (QCONV(J)-QNEW)/(G*DELT*ROWL)
c          TCONV(J)=TNEW
c          QCONV(J)=QNEW
c          QSCONV(J)=QNEW
c         END IF     
c  38  CONTINUE  
c     END IF
C
c     END IF
C
C  *** CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY AND STATIC ENERGY
C  
        GZ(1)=0.0
        CPN(1)=CPD*(1.-QCONV(1))+QCONV(1)*CPV
        H(1)=TCONV(1)*CPN(1)
        LV(1)=LV0-CPVMCL*(TCONV(1)-273.15)
        HM(1)=LV(1)*QCONV(1)
        TV(1)=TCONV(1)*(1.+QCONV(1)*EPSI-QCONV(1))
        AHMIN=1.0E12
        IHMIN=NL
        DO 40 I=2,NL+1
          TVX=TCONV(I)*(1.+QCONV(I)*EPSI-QCONV(I))
          TVY=TCONV(I-1)*(1.+QCONV(I-1)*EPSI-QCONV(I-1))
          GZ(I)=GZ(I-1)+0.5*RD*(TVX+TVY)*(PCONV_HPA(I-1)-PCONV_HPA(I))/
     1    PHCONV_HPA(I)
          CPN(I)=CPD*(1.-QCONV(I))+CPV*QCONV(I)
          H(I)=TCONV(I)*CPN(I)+GZ(I)
          LV(I)=LV0-CPVMCL*(TCONV(I)-273.15)
          HM(I)=(CPD*(1.-QCONV(I))+CL*QCONV(I))*(TCONV(I)-TCONV(1))+
     1     LV(I)*QCONV(I)+GZ(I)
          TV(I)=TCONV(I)*(1.+QCONV(I)*EPSI-QCONV(I))
C
C  ***  Find level of minimum moist static energy    ***
C
          IF(I.GE.MINORIG.AND.HM(I).LT.AHMIN.AND.HM(I).LT.HM(I-1))THEN
           AHMIN=HM(I)
           IHMIN=I
          END IF
   40   CONTINUE
        IHMIN=MIN(IHMIN, NL-1)
C
C  ***     Find that model level below the level of minimum moist       ***
C  ***  static energy that has the maximum value of moist static energy ***
C
        AHMAX=0.0
        DO 42 I=MINORIG,IHMIN
         IF(HM(I).GT.AHMAX)THEN
          NK=I
          AHMAX=HM(I)
         END IF
   42   CONTINUE
C
C  ***  CHECK WHETHER PARCEL LEVEL TEMPERATURE AND SPECIFIC HUMIDITY   ***
C  ***                          ARE REASONABLE                         ***
C  ***      Skip convection if HM increases monotonically upward       ***
C
        IF(TCONV(NK).LT.250.0.OR.QCONV(NK).LE.0.0.OR.IHMIN.EQ.(NL-1))
     1  THEN
         IFLAG=0
         CBMF=0.0
         RETURN
        END IF
C
C   ***  CALCULATE LIFTED CONDENSATION LEVEL OF AIR AT PARCEL ORIGIN LEVEL ***
C   ***       (WITHIN 0.2% OF FORMULA OF BOLTON, MON. WEA. REV.,1980)      ***
C
        RH=QCONV(NK)/QSCONV(NK)
        CHI=TCONV(NK)/(1669.0-122.0*RH-TCONV(NK))
        PLCL=PCONV_HPA(NK)*(RH**CHI)
        IF(PLCL.LT.200.0.OR.PLCL.GE.2000.0)THEN
         IFLAG=2
         CBMF=0.0
         RETURN
        END IF
C
C   ***  CALCULATE FIRST LEVEL ABOVE LCL (=ICB)  ***
C
        ICB=NL-1
        DO 50 I=NK+1,NL
         IF(PCONV_HPA(I).LT.PLCL)THEN
          ICB=MIN(ICB,I)
         END IF
   50   CONTINUE
        IF(ICB.GE.(NL-1))THEN
         IFLAG=3
         CBMF=0.0
         RETURN
        END IF
C
C   *** FIND TEMPERATURE UP THROUGH ICB AND TEST FOR INSTABILITY           ***
C
C   *** SUBROUTINE TLIFT CALCULATES PART OF THE LIFTED PARCEL VIRTUAL      ***
C   ***  TEMPERATURE, THE ACTUAL TEMPERATURE AND THE ADIABATIC             ***
C   ***                   LIQUID WATER CONTENT                             ***
C
        CALL TLIFT(GZ,ICB,NK,TVP,TP,CLW,ND,NL,1)
        DO 54 I=NK,ICB
         TVP(I)=TVP(I)-TP(I)*QCONV(NK)
   54   CONTINUE
C
C   ***  If there was no convection at last time step and parcel    ***
C   ***       is stable at ICB then skip rest of calculation        ***
C
        IF(CBMF.EQ.0.0.AND.TVP(ICB).LE.(TV(ICB)-DTMAX))THEN
         IFLAG=0
         RETURN
        END IF
C
C   ***  IF THIS POINT IS REACHED, MOIST CONVECTIVE ADJUSTMENT IS NECESSARY ***
C
        IF(IFLAG.NE.4)IFLAG=1
C
C   ***  FIND THE REST OF THE LIFTED PARCEL TEMPERATURES          ***
C
        CALL TLIFT(GZ,ICB,NK,TVP,TP,CLW,ND,NL,2)
C
C   ***  SET THE PRECIPITATION EFFICIENCIES AND THE FRACTION OF   ***
C   ***          PRECIPITATION FALLING OUTSIDE OF CLOUD           ***
C   ***      THESE MAY BE FUNCTIONS OF TP(I), PCONV_HPA(I) AND CLW(I)     ***
C                 
        DO 57 I=1,NK
         EP(I)=0.0
         SIGP(I)=SIGS
   57   CONTINUE
        DO 60 I=NK+1,NL
         TCA=TP(I)-273.15
         IF(TCA.GE.0.0)THEN
          ELACRIT=ELCRIT
         ELSE
          ELACRIT=ELCRIT*(1.0-TCA/TLCRIT)
         END IF
         ELACRIT=MAX(ELACRIT,0.0)
	   EPMAX=0.999
         EP(I)=EPMAX*(1.0-ELACRIT/MAX(CLW(I),1.0E-8))
         EP(I)=MAX(EP(I),0.0)
         EP(I)=MIN(EP(I),EPMAX)
         SIGP(I)=SIGS
   60   CONTINUE
C
C   ***       CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL     ***
C   ***                    VIRTUAL TEMPERATURE                    ***
C
        DO 64 I=ICB+1,NL
         TVP(I)=TVP(I)-TP(I)*QCONV(NK)
   64   CONTINUE
        TVP(NL+1)=TVP(NL)-(GZ(NL+1)-GZ(NL))/CPD
C
C   ***        NOW INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS       ***
C
        DO 70 I=1,NL+1
         HP(I)=H(I)
         NENT(I)=0
         WATER(I)=0.0
         EVAP(I)=0.0
         WT(I)=OMTSNOW
         LVCP(I)=LV(I)/CPN(I)
         DO 70 J=1,NL+1
          QENT(I,J)=QCONV(J)
          ELIJ(I,J)=0.0
          SIJ(I,J)=0.0
   70   CONTINUE
        QP(1)=QCONV(1)
        DO 72 I=2,NL+1
         QP(I)=QCONV(I-1)
   72	   CONTINUE
C
C  ***  FIND THE FIRST MODEL LEVEL (INB1) ABOVE THE PARCEL'S      ***
C  ***          HIGHEST LEVEL OF NEUTRAL BUOYANCY                 ***
C  ***     AND THE HIGHEST LEVEL OF POSITIVE CAPE (INB)           ***
C
        CAPE=0.0
        CAPEM=0.0
        INB=ICB+1
        INB1=INB
	  BYP=0.0
        DO 82 I=ICB+1,NL-1
         BY=(TVP(I)-TV(I))*(PHCONV_HPA(I)-PHCONV_HPA(I+1))/PCONV_HPA(I)
         CAPE=CAPE+BY
         IF(BY.GE.0.0)INB1=I+1
         IF(CAPE.GT.0.0)THEN
          INB=I+1
          BYP=(TVP(I+1)-TV(I+1))*(PHCONV_HPA(I+1)-PHCONV_HPA(I+2))/
     1    PCONV_HPA(I+1)
          CAPEM=CAPE
         END IF
   82 	CONTINUE
        INB=MAX(INB,INB1)
        CAPE=CAPEM+BYP
        DEFRAC=CAPEM-CAPE
        DEFRAC=MAX(DEFRAC,0.001)
        FRAC=-CAPE/DEFRAC
        FRAC=MIN(FRAC,1.0)
        FRAC=MAX(FRAC,0.0)
C
C   ***   CALCULATE LIQUID WATER STATIC ENERGY OF LIFTED PARCEL   ***
C
        DO 95 I=ICB,INB
         HP(I)=H(NK)+(LV(I)+(CPD-CPV)*TCONV(I))*EP(I)*CLW(I)
   95   CONTINUE                  
C
C   ***  CALCULATE CLOUD BASE MASS FLUX AND RATES OF MIXING, M(I),  ***
c   ***                   AT EACH MODEL LEVEL                       ***
C
        DBOSUM=0.0
C   
C   ***     INTERPOLATE DIFFERENCE BETWEEN LIFTED PARCEL AND      ***
C   ***  ENVIRONMENTAL TEMPERATURES TO LIFTED CONDENSATION LEVEL  ***
C	
        TVPPLCL=TVP(ICB-1)-RD*TVP(ICB-1)*(PCONV_HPA(ICB-1)-PLCL)/
     1    (CPN(ICB-1)*PCONV_HPA(ICB-1))
        TVAPLCL=TV(ICB)+(TVP(ICB)-TVP(ICB+1))*(PLCL-PCONV_HPA(ICB))/
     1    (PCONV_HPA(ICB)-PCONV_HPA(ICB+1))
        DTPBL=0.0
        DO 96 I=NK,ICB-1
         DTPBL=DTPBL+(TVP(I)-TV(I))*(PHCONV_HPA(I)-PHCONV_HPA(I+1))
   96   CONTINUE
        DTPBL=DTPBL/(PHCONV_HPA(NK)-PHCONV_HPA(ICB))
        DTMIN=TVPPLCL-TVAPLCL+DTMAX+DTPBL
        DTMA=DTMIN
C
C   ***  ADJUST CLOUD BASE MASS FLUX   ***
C
      CBMFOLD=CBMF
	DELT0=300.0
      DAMPS=DAMP*DELT/DELT0 
      CBMF=(1.-DAMPS)*CBMF+0.1*ALPHA*DTMA 
      CBMF=MAX(CBMF,0.0)
C
C   *** If cloud base mass flux is zero, skip rest of calculation  ***
C
      IF(CBMF.EQ.0.0.AND.CBMFOLD.EQ.0.0)THEN
       RETURN
      END IF

C
C   ***   CALCULATE RATES OF MIXING,  M(I)   ***
C
      M(ICB)=0.0
      DO 103 I=ICB+1,INB
       K=MIN(I,INB1)
       DBO=ABS(TV(K)-TVP(K))+
     1  ENTP*0.02*(PHCONV_HPA(K)-PHCONV_HPA(K+1))
       DBOSUM=DBOSUM+DBO
       M(I)=CBMF*DBO
  103 CONTINUE
      DO 110 I=ICB+1,INB
       M(I)=M(I)/DBOSUM  
  110 CONTINUE     
C
C   ***  CALCULATE ENTRAINED AIR MASS FLUX (MENT), TOTAL WATER MIXING  ***
C   ***     RATIO (QENT), TOTAL CONDENSED WATER (ELIJ), AND MIXING     ***
C   ***                        FRACTION (SIJ)                          ***
C
        DO 170 I=ICB+1,INB
         QTI=QCONV(NK)-EP(I)*CLW(I)
         DO 160 J=ICB,INB
          BF2=1.+LV(J)*LV(J)*QSCONV(J)/(RV*TCONV(J)*TCONV(J)*CPD)
          ANUM=H(J)-HP(I)+(CPV-CPD)*TCONV(J)*(QTI-QCONV(J))
          DENOM=H(I)-HP(I)+(CPD-CPV)*(QCONV(I)-QTI)*TCONV(J)
          DEI=DENOM
          IF(ABS(DEI).LT.0.01)DEI=0.01
          SIJ(I,J)=ANUM/DEI
          SIJ(I,I)=1.0
          ALTEM=SIJ(I,J)*QCONV(I)+(1.-SIJ(I,J))*QTI-QSCONV(J)
          ALTEM=ALTEM/BF2
          CWAT=CLW(J)*(1.-EP(J))
          STEMP=SIJ(I,J)
          IF((STEMP.LT.0.0.OR.STEMP.GT.1.0.OR.
     1      ALTEM.GT.CWAT).AND.J.GT.I)THEN
           ANUM=ANUM-LV(J)*(QTI-QSCONV(J)-CWAT*BF2)
           DENOM=DENOM+LV(J)*(QCONV(I)-QTI)
           IF(ABS(DENOM).LT.0.01)DENOM=0.01
           SIJ(I,J)=ANUM/DENOM
           ALTEM=SIJ(I,J)*QCONV(I)+(1.-SIJ(I,J))*QTI-QSCONV(J)
           ALTEM=ALTEM-(BF2-1.)*CWAT
          END IF
          IF(SIJ(I,J).GT.0.0.AND.SIJ(I,J).LT.0.9)THEN
           QENT(I,J)=SIJ(I,J)*QCONV(I)+(1.-SIJ(I,J))*QTI
           ELIJ(I,J)=ALTEM
           ELIJ(I,J)=MAX(0.0,ELIJ(I,J))
           MENT(I,J)=M(I)/(1.-SIJ(I,J))
           NENT(I)=NENT(I)+1
          END IF
          SIJ(I,J)=MAX(0.0,SIJ(I,J))
          SIJ(I,J)=MIN(1.0,SIJ(I,J))
  160    CONTINUE
C
C   ***   IF NO AIR CAN ENTRAIN AT LEVEL I ASSUME THAT UPDRAFT DETRAINS  ***
C   ***   AT THAT LEVEL AND CALCULATE DETRAINED AIR FLUX AND PROPERTIES  ***
C
         IF(NENT(I).EQ.0)THEN
          MENT(I,I)=M(I)
          QENT(I,I)=QCONV(NK)-EP(I)*CLW(I)
          ELIJ(I,I)=CLW(I)
          SIJ(I,I)=1.0
         END IF 
  170   CONTINUE
        SIJ(INB,INB)=1.0
C
C   ***  NORMALIZE ENTRAINED AIR MASS FLUXES TO REPRESENT EQUAL  ***
C   ***              PROBABILITIES OF MIXING                     ***
C
        DO 200 I=ICB+1,INB
        IF(NENT(I).NE.0)THEN
         QP1=QCONV(NK)-EP(I)*CLW(I)
         ANUM=H(I)-HP(I)-LV(I)*(QP1-QSCONV(I))
         DENOM=H(I)-HP(I)+LV(I)*(QCONV(I)-QP1)
         IF(ABS(DENOM).LT.0.01)DENOM=0.01
         SCRIT=ANUM/DENOM
         ALT=QP1-QSCONV(I)+SCRIT*(QCONV(I)-QP1)
         IF(ALT.LT.0.0)SCRIT=1.0
	   SCRIT=MAX(SCRIT,0.0)
         ASIJ=0.0
         SMIN=1.0
         DO 175 J=ICB,INB
          IF(SIJ(I,J).GT.0.0.AND.SIJ(I,J).LT.0.9)THEN
           IF(J.GT.I)THEN
            SMID=MIN(SIJ(I,J),SCRIT)
            SJMAX=SMID
            SJMIN=SMID
            IF(SMID.LT.SMIN.AND.SIJ(I,J+1).LT.SMID)THEN
             SMIN=SMID
             SJMAX=MIN(SIJ(I,J+1),SIJ(I,J),SCRIT)
             SJMIN=MAX(SIJ(I,J-1),SIJ(I,J))
             SJMIN=MIN(SJMIN,SCRIT)
            END IF
           ELSE
            SJMAX=MAX(SIJ(I,J+1),SCRIT)
            SMID=MAX(SIJ(I,J),SCRIT)
            SJMIN=0.0
            IF(J.GT.1)SJMIN=SIJ(I,J-1)
            SJMIN=MAX(SJMIN,SCRIT)
           END IF
           DELP=ABS(SJMAX-SMID)
           DELM=ABS(SJMIN-SMID)
           ASIJ=ASIJ+(DELP+DELM)*(PHCONV_HPA(J)-PHCONV_HPA(J+1))
           MENT(I,J)=MENT(I,J)*(DELP+DELM)*
     1     (PHCONV_HPA(J)-PHCONV_HPA(J+1))
          END IF
  175    CONTINUE
         ASIJ=MAX(1.0E-21,ASIJ)
         ASIJ=1.0/ASIJ
         DO 180 J=ICB,INB
          MENT(I,J)=MENT(I,J)*ASIJ
  180    CONTINUE
         BSUM=0.0
         DO 190 J=ICB,INB
          BSUM=BSUM+MENT(I,J)
  190    CONTINUE
         IF(BSUM.LT.1.0E-18)THEN
          NENT(I)=0
          MENT(I,I)=M(I)
          QENT(I,I)=QCONV(NK)-EP(I)*CLW(I)
          ELIJ(I,I)=CLW(I)
          SIJ(I,I)=1.0
         END IF
        END IF
  200   CONTINUE
C
C   ***  CHECK WHETHER EP(INB)=0, IF SO, SKIP PRECIPITATING    ***
C   ***             DOWNDRAFT CALCULATION                      ***
C
        IF(EP(INB).LT.0.0001)GOTO 405
C
C   ***  INTEGRATE LIQUID WATER EQUATION TO FIND CONDENSED WATER   ***
C   ***                AND CONDENSED WATER FLUX                    ***
C
        JTT=2
C
C    ***                    BEGIN DOWNDRAFT LOOP                    ***
C
        DO 400 I=INB,1,-1
C
C    ***              CALCULATE DETRAINED PRECIPITATION             ***
C
        WDTRAIN=G*EP(I)*M(I)*CLW(I)
        IF(I.GT.1)THEN
         DO 320 J=1,I-1
         AWAT=ELIJ(J,I)-(1.-EP(I))*CLW(I)
         AWAT=MAX(0.0,AWAT)
  320    WDTRAIN=WDTRAIN+G*AWAT*MENT(J,I)
        END IF
C
C    ***    FIND RAIN WATER AND EVAPORATION USING PROVISIONAL   ***
C    ***              ESTIMATES OF QP(I)AND QP(I-1)             ***
C     
c
c  ***  Value of terminal velocity and coefficient of evaporation for snow   ***
c 
        COEFF=COEFFS
        WT(I)=OMTSNOW
c      
c  ***  Value of terminal velocity and coefficient of evaporation for rain   ***
c
        IF(TCONV(I).GT.273.0)THEN
         COEFF=COEFFR
         WT(I)=OMTRAIN
        END IF
        QSM=0.5*(QCONV(I)+QP(I+1))
        AFAC=COEFF*PHCONV_HPA(I)*(QSCONV(I)-QSM)/
     1  (1.0E4+2.0E3*PHCONV_HPA(I)*QSCONV(I))
        AFAC=MAX(AFAC,0.0)
        SIGT=SIGP(I)
        SIGT=MAX(0.0,SIGT)
        SIGT=MIN(1.0,SIGT)
        B6=100.*(PHCONV_HPA(I)-PHCONV_HPA(I+1))*SIGT*AFAC/WT(I)
        C6=(WATER(I+1)*WT(I+1)+WDTRAIN/SIGD)/WT(I)
        REVAP=0.5*(-B6+SQRT(B6*B6+4.*C6))
        EVAP(I)=SIGT*AFAC*REVAP
        WATER(I)=REVAP*REVAP
C
C    ***  CALCULATE PRECIPITATING DOWNDRAFT MASS FLUX UNDER     ***
C    ***              HYDROSTATIC APPROXIMATION                 ***
C   
        IF(I.EQ.1)GOTO 360
        DHDP=(H(I)-H(I-1))/(PCONV_HPA(I-1)-PCONV_HPA(I))
        DHDP=MAX(DHDP,10.0)
        MP(I)=100.*GINV*LV(I)*SIGD*EVAP(I)/DHDP
        MP(I)=MAX(MP(I),0.0)
C
C   ***   ADD SMALL AMOUNT OF INERTIA TO DOWNDRAFT              ***
C
        FAC=20.0/(PHCONV_HPA(I-1)-PHCONV_HPA(I))
        MP(I)=(FAC*MP(I+1)+MP(I))/(1.+FAC)
C   
C    ***      FORCE MP TO DECREASE LINEARLY TO ZERO                 ***
C    ***      BETWEEN ABOUT 950 MB AND THE SURFACE                  ***
C
          IF(PCONV_HPA(I).GT.(0.949*PCONV_HPA(1)))THEN
           JTT=MAX(JTT,I)
           MP(I)=MP(JTT)*(PCONV_HPA(1)-PCONV_HPA(I))/(PCONV_HPA(1)-
     1     PCONV_HPA(JTT))
          END IF              
  360   CONTINUE
C
C    ***       FIND MIXING RATIO OF PRECIPITATING DOWNDRAFT     ***
C
        IF(I.EQ.INB)GOTO 400
        IF(I.EQ.1)THEN
         QSTM=QSCONV(1)
        ELSE
         QSTM=QSCONV(I-1)
        END IF
        IF(MP(I).GT.MP(I+1))THEN
          RAT=MP(I+1)/MP(I)
          QP(I)=QP(I+1)*RAT+QCONV(I)*(1.0-RAT)+100.*GINV*
     1       SIGD*(PHCONV_HPA(I)-PHCONV_HPA(I+1))*(EVAP(I)/MP(I))
         ELSE
          IF(MP(I+1).GT.0.0)THEN
            QP(I)=(GZ(I+1)-GZ(I)+QP(I+1)*(LV(I+1)+TCONV(I+1)*(CL-CPD))+
     1      CPD*(TCONV(I+1)-TCONV(I)))/(LV(I)+TCONV(I)*(CL-CPD))
          END IF
        END IF
        QP(I)=MIN(QP(I),QSTM)
        QP(I)=MAX(QP(I),0.0)
  400   CONTINUE
C
C   ***  CALCULATE SURFACE PRECIPITATION IN MM/DAY     ***
C
        PRECIP=PRECIP+WT(1)*SIGD*WATER(1)*3600.*24000./(ROWL*G)
C
  405   CONTINUE
C
C   ***  CALCULATE DOWNDRAFT VELOCITY SCALE AND SURFACE TEMPERATURE AND  ***
c   ***                    WATER VAPOR FLUCTUATIONS                      ***
C
      WD=BETA*ABS(MP(ICB))*0.01*RD*TCONV(ICB)/(SIGD*PCONV_HPA(ICB))
      QPRIME=0.5*(QP(1)-QCONV(1))
      TPRIME=LV0*QPRIME/CPD
C
C   ***  CALCULATE TENDENCIES OF LOWEST LEVEL POTENTIAL TEMPERATURE  ***
C   ***                      AND MIXING RATIO                        ***
C
        DPINV=0.01/(PHCONV_HPA(1)-PHCONV_HPA(2))
        AM=0.0
        IF(NK.EQ.1)THEN
         DO 410 K=2,INB
  410    AM=AM+M(K)
        END IF
c save saturated upward mass flux for first level
        FUP(1)=AM
        IF((2.*G*DPINV*AM).GE.DELTI)IFLAG=4
        FT(1)=FT(1)+G*DPINV*AM*(TCONV(2)-TCONV(1)+(GZ(2)-GZ(1))/CPN(1))
        FT(1)=FT(1)-LVCP(1)*SIGD*EVAP(1)
        FT(1)=FT(1)+SIGD*WT(2)*(CL-CPD)*WATER(2)*(TCONV(2)-
     1   TCONV(1))*DPINV/CPN(1)
        FQ(1)=FQ(1)+G*MP(2)*(QP(2)-QCONV(1))*
     1    DPINV+SIGD*EVAP(1)
        FQ(1)=FQ(1)+G*AM*(QCONV(2)-QCONV(1))*DPINV
        DO 415 J=2,INB
         FQ(1)=FQ(1)+G*DPINV*MENT(J,1)*(QENT(J,1)-QCONV(1))
  415   CONTINUE
C
C   ***  CALCULATE TENDENCIES OF POTENTIAL TEMPERATURE AND MIXING RATIO  ***
C   ***               AT LEVELS ABOVE THE LOWEST LEVEL                   ***
C
C   ***  FIRST FIND THE NET SATURATED UPDRAFT AND DOWNDRAFT MASS FLUXES  ***
C   ***                      THROUGH EACH LEVEL                          ***
C
        DO 500 I=2,INB
        DPINV=0.01/(PHCONV_HPA(I)-PHCONV_HPA(I+1))
        CPINV=1.0/CPN(I)
        AMP1=0.0
        AD=0.0
        IF(I.GE.NK)THEN
         DO 440 K=I+1,INB+1
  440    AMP1=AMP1+M(K)
        END IF
        DO 450 K=1,I
        DO 450 J=I+1,INB+1
         AMP1=AMP1+MENT(K,J)
  450   CONTINUE
c save saturated upward mass flux
        FUP(I)=AMP1
        IF((2.*G*DPINV*AMP1).GE.DELTI)IFLAG=4
        DO 470 K=1,I-1
        DO 470 J=I,INB
         AD=AD+MENT(J,K)
  470   CONTINUE
c save saturated downward mass flux
        FDOWN(I)=AD
        FT(I)=FT(I)+G*DPINV*(AMP1*(TCONV(I+1)-TCONV(I)+(GZ(I+1)-GZ(I))*
     1   CPINV)-AD*(TCONV(I)-TCONV(I-1)+(GZ(I)-GZ(I-1))*CPINV))
     2   -SIGD*LVCP(I)*EVAP(I)
        FT(I)=FT(I)+G*DPINV*MENT(I,I)*(HP(I)-H(I)+
     1    TCONV(I)*(CPV-CPD)*(QCONV(I)-QENT(I,I)))*CPINV
        FT(I)=FT(I)+SIGD*WT(I+1)*(CL-CPD)*WATER(I+1)*
     1    (TCONV(I+1)-TCONV(I))*DPINV*CPINV
        FQ(I)=FQ(I)+G*DPINV*(AMP1*(QCONV(I+1)-QCONV(I))-
     1    AD*(QCONV(I)-QCONV(I-1)))
        DO 480 K=1,I-1
         AWAT=ELIJ(K,I)-(1.-EP(I))*CLW(I)
         AWAT=MAX(AWAT,0.0)
         FQ(I)=FQ(I)+G*DPINV*MENT(K,I)*(QENT(K,I)-AWAT-QCONV(I))
  480   CONTINUE
        DO 490 K=I,INB
         FQ(I)=FQ(I)+G*DPINV*MENT(K,I)*(QENT(K,I)-QCONV(I))
  490   CONTINUE
        FQ(I)=FQ(I)+SIGD*EVAP(I)+G*(MP(I+1)*
     1    (QP(I+1)-QCONV(I))-MP(I)*(QP(I)-QCONV(I-1)))*DPINV
  500   CONTINUE
C
C   *** Adjust tendencies at top of convection layer to reflect  ***
C   ***       actual position of the level zero CAPE             ***
C
        FQOLD=FQ(INB)
        FQ(INB)=FQ(INB)*(1.-FRAC)
        FQ(INB-1)=FQ(INB-1)+FRAC*FQOLD*((PHCONV_HPA(INB)-
     1    PHCONV_HPA(INB+1))/
     2   (PHCONV_HPA(INB-1)-PHCONV_HPA(INB)))*LV(INB)/LV(INB-1)
        FTOLD=FT(INB)
        FT(INB)=FT(INB)*(1.-FRAC)
        FT(INB-1)=FT(INB-1)+FRAC*FTOLD*((PHCONV_HPA(INB)-
     1  PHCONV_HPA(INB+1))/
     2   (PHCONV_HPA(INB-1)-PHCONV_HPA(INB)))*CPN(INB)/CPN(INB-1)
C
C   ***   Very slightly adjust tendencies to force exact   ***
C   ***     enthalpy, momentum and tracer conservation     ***
C
        ENTS=0.0
        DO 680 I=1,INB
         ENTS=ENTS+(CPN(I)*FT(I)+LV(I)*FQ(I))*
     +   (PHCONV_HPA(I)-PHCONV_HPA(I+1))	
  680	CONTINUE
        ENTS=ENTS/(PHCONV_HPA(1)-PHCONV_HPA(INB+1))
        DO 640 I=1,INB
         FT(I)=FT(I)-ENTS/CPN(I)
  640	CONTINUE

c ************************************************
c **** DETERMINE MASS DISPLACEMENT MATRIX
c ***** AND COMPENSATING SUBSIDENCE
c ************************************************

c mass displacement matrix due to saturated up-and downdrafts
c inside the cloud and determine compensating subsidence
c FUP(I) (saturated updrafts), FDOWN(I) (saturated downdrafts) are assumed to be
c balanced by  compensating subsidence (SUB(I))
c FDOWN(I) and SUB(I) defined positive downwards

C NCONVTOP IS THE TOP LEVEL AT WHICH CONVECTIVE MASS FLUXES ARE DIAGNOSED
C EPSILON IS A SMALL NUMBER

       SUB(1)=0.
       NCONVTOP=1
       do i=1,INB+1
       do j=1,INB+1
        if (j.eq.NK) then
         FMASS(j,i)=FMASS(j,i)+M(i)
        endif
         FMASS(j,i)=FMASS(j,i)+MENT(j,i)
         IF (FMASS(J,I).GT.EPSILON) NCONVTOP=MAX(NCONVTOP,I,J)
       end do
       if (i.gt.1) then
        SUB(i)=FUP(i-1)-FDOWN(i)
       endif
       end do
       NCONVTOP=NCONVTOP+1

        RETURN
C
        END
C
C ---------------------------------------------------------------------------
C
        SUBROUTINE TLIFT(GZ,ICB,NK,TVP,TPK,CLW,ND,NL,KK)
c
c-cv
        include 'includepar'
        include 'includeconv'
c-cv
C====>Begin Module TLIFT      File convect.f      Undeclared variables
C
C     Argument variables
C
      integer icb, kk, nd, nk, nl
C
C     Local variables
C
      integer i, j, nsb, nst
C
      real ah0, ahg, alv, cl, cpd, cpinv, cpp, cpv, cpvmcl, denom, epsi
      real eps0, es, qg, rd, rg, rv, s, tc, tg
C
C====>End Module   TLIFT      File convect.f

        REAL GZ(ND),TPK(ND),CLW(ND)
        REAL TVP(ND),LV0
C
C   ***   ASSIGN VALUES OF THERMODYNAMIC CONSTANTS     ***
C
        PARAMETER(CPD=1005.7)
        PARAMETER(CPV=1870.0)
        PARAMETER(CL=2500.0)
        PARAMETER(RV=461.5)
        PARAMETER(RD=287.04)
        PARAMETER(LV0=2.501E6)
C
        PARAMETER(CPVMCL=CL-CPV)
        PARAMETER(EPS0=RD/RV)
        PARAMETER(EPSI=1./EPS0)
C
C   ***  CALCULATE CERTAIN PARCEL QUANTITIES, INCLUDING STATIC ENERGY   ***
C
        AH0=(CPD*(1.-QCONV(NK))+CL*QCONV(NK))*TCONV(NK)+QCONV(NK)*
     1  (LV0-CPVMCL*(
     1   TCONV(NK)-273.15))+GZ(NK)
        CPP=CPD*(1.-QCONV(NK))+QCONV(NK)*CPV
        CPINV=1./CPP
C
        IF(KK.EQ.1)THEN
C
C   ***   CALCULATE LIFTED PARCEL QUANTITIES BELOW CLOUD BASE   ***
C
        DO 50 I=1,ICB-1
         CLW(I)=0.0
   50   CONTINUE
        DO 100 I=NK,ICB-1
         TPK(I)=TCONV(NK)-(GZ(I)-GZ(NK))*CPINV
         TVP(I)=TPK(I)*(1.+QCONV(NK)*EPSI)
  100   CONTINUE
        END IF
C
C    ***  FIND LIFTED PARCEL QUANTITIES ABOVE CLOUD BASE    ***
C
        NST=ICB
        NSB=ICB
        IF(KK.EQ.2)THEN  
         NST=NL
         NSB=ICB+1
        END IF
        DO 300 I=NSB,NST
         TG=TCONV(I)
         QG=QSCONV(I)
         ALV=LV0-CPVMCL*(TCONV(I)-273.15)
         DO 200 J=1,2
          S=CPD+ALV*ALV*QG/(RV*TCONV(I)*TCONV(I))
          S=1./S
          AHG=CPD*TG+(CL-CPD)*QCONV(NK)*TCONV(I)+ALV*QG+GZ(I)
          TG=TG+S*(AH0-AHG)
          TG=MAX(TG,35.0)
          TC=TG-273.15
          DENOM=243.5+TC
          IF(TC.GE.0.0)THEN  
           ES=6.112*EXP(17.67*TC/DENOM)
          ELSE  
           ES=EXP(23.33086-6111.72784/TG+0.15215*LOG(TG))
          END IF  
          QG=EPS0*ES/(PCONV_HPA(I)-ES*(1.-EPS0))
  200    CONTINUE
         ALV=LV0-CPVMCL*(TCONV(I)-273.15)
         TPK(I)=(AH0-(CL-CPD)*QCONV(NK)*TCONV(I)-GZ(I)-ALV*QG)/CPD
         CLW(I)=QCONV(NK)-QG
         CLW(I)=MAX(0.0,CLW(I))
         RG=QG/(1.-QCONV(NK))
         TVP(I)=TPK(I)*(1.+RG*EPSI)
  300   CONTINUE
        RETURN
        END
