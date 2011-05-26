C  Taken from Press et al., Numerical Recipes

      FUNCTION RAN1(IDUM)
      INTEGER IDUM,IA,IM,IQ,IR,NTAB,NDIV
      REAL RAN1,AM,EPS,RNMX
      PARAMETER(IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     +NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2E-7,RNMX=1.-EPS)
      INTEGER J,K,IV(NTAB),IY
      SAVE IV,IY
      DATA IV/NTAB*0/,IY/0/
      IF (IDUM.LE.0.OR.IY.EQ.0) THEN
        IDUM=MAX(-IDUM,1)
        DO J=NTAB+8,1,-1
          K=IDUM/IQ
          IDUM=IA*(IDUM-K*IQ)-IR*K
          IF (IDUM.LT.0) IDUM=IDUM+IM
          IF (J.LE.NTAB) IV(J)=IDUM
        ENDDO
        IY=IV(1)
      ENDIF
      K=IDUM/IQ
      IDUM=IA*(IDUM-K*IQ)-IR*K
      IF (IDUM.LT.0) IDUM=IDUM+IM
      J=1+IY/NDIV
      IY=IV(J)
      IV(J)=IDUM
      RAN1=MIN(AM*IY,RNMX)
      RETURN
      END



      FUNCTION GASDEV(IDUM)
      INTEGER IDUM,ISET
      REAL GASDEV
      DATA ISET/0/
      SAVE ISET,GSET
      IF (ISET.EQ.0) THEN
1       V1=2.*RAN3(IDUM)-1.
        V2=2.*RAN3(IDUM)-1.
        R=V1**2+V2**2
        IF(R.GE.1.0 .OR. R.EQ.0.0) GO TO 1
        FAC=SQRT(-2.*LOG(R)/R)
        GSET=V1*FAC
        GASDEV=V2*FAC
        ISET=1
      ELSE
        GASDEV=GSET
        ISET=0
      ENDIF
      RETURN
      END

      SUBROUTINE GASDEV1(IDUM,RANDOM1,RANDOM2)
      INTEGER IDUM
1     V1=2.*RAN3(IDUM)-1.
      V2=2.*RAN3(IDUM)-1.
      R=V1**2+V2**2
      IF(R.GE.1.0 .OR. R.EQ.0.0) GO TO 1
      FAC=SQRT(-2.*LOG(R)/R)
      RANDOM1=V1*FAC
      RANDOM2=V2*FAC
      RETURN
      END


      FUNCTION ran3(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software us.
