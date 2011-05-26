      subroutine calcmatrix(lconv,delt,cbmf)
c                             o    i    o
C******************************************************************
C     This subroutine calculates the matrix describing convective 
C     redistribution of mass in a grid column, using the subroutine
C     convect43c.f provided by Kerry Emanuel.
c
C     Petra Seibert, Bernd C. Krueger, 2000-2001
c
c     changed by C. Forster, November 2003 - February 2004
c     array fmassfrac(nconvlevmax,nconvlevmax) represents
c     the convective redistribution matrix for the particles
c
c     20 Oct 2005 - R. Easter - added calc of pconv_hpa(nconvlev+1)
c     16 Nov 2005 - R. Easter - pconv,phconv are set in convmix
c                               using pph & pphn
c
C******************************************************************

C lconv        indicates whether there is convection in this cell, or not
C delt         time step for convection [s]
C cbmf         cloud base mass flux

      include 'includepar'
      include 'includecom'
      include 'includeconv'

      real rlevmass,summe

      integer iflag, k, kk, kuvz
     
C     1-d variables for convection
C     variables for redistribution matrix
      real cbmfold, precip, qprime
      real tprime, wd, f_qvsat
      real delt,cbmf
      logical lconv

      lconv = .false.

   
C calculate pressure at eta levels for use in convect
C and assign temp & spec. hum. to 1D workspace
C -------------------------------------------------------

C pconv(1) is the pressure at the first level above ground
C phconv(k) is the pressure between levels k-1 and k 
C dp(k) is the pressure difference "around" tconv(k)
C phconv(kmax) must also be defined 1/2 level above pconv(kmax)
C Therefore, we define k = kuvz-1 and let kuvz start from 2
C top layer cannot be used for convection because p at top of this layer is
C not given


c     phconv(1) = psconv
C Emanuel subroutine needs pressure in hPa, therefore convert all pressures       
      do kuvz = 2,nuvz
        k = kuvz-1
c       pconv(k) = (akz(kuvz) + bkz(kuvz)*psconv)
c       phconv(kuvz) = (akm(kuvz) + bkm(kuvz)*psconv)
        dp(k) = phconv(k) - phconv(kuvz)
        qsconv(k) = f_qvsat( pconv(k), tconv(k) )

c initialize mass fractions 
        do kk=1,nconvlev
          fmassfrac(k,kk)=0.
        enddo
      enddo


C     note that Emanuel says it is important
C     a. to set this =0. every grid point
C     b. to keep this value in the calling programme in the iteration

c CALL CONVECTION
c******************

        cbmfold = cbmf
C Convert pressures to hPa, as required by Emanuel scheme
*********************************************************
        do 10 k=1,nconvlev
          pconv_hpa(k)=pconv(k)/100.
10        phconv_hpa(k)=phconv(k)/100.
        pconv_hpa(nconvlev+1)=pconv(nconvlev+1)/100.
        phconv_hpa(nconvlev+1)=phconv(nconvlev+1)/100.
        call convect(nconvlevmax, nconvlev, delt, iflag,
     &            precip, wd, tprime, qprime, cbmf)

c      do not update fmassfrac and cloudbase massflux
c      if no convection takes place or
c      if a CFL criterion is violated in convect43c.f
       if (iflag .ne. 1 .and. iflag .ne. 4) then
         cbmf=cbmfold
         goto 200
       endif

c      do not update fmassfrac and cloudbase massflux
c      if the old and the new cloud base mass
c      fluxes are zero
       if (cbmf.le.0..and.cbmfold.le.0.) then
         cbmf=cbmfold
         goto 200
       endif

c      Update fmassfrac
c      account for mass displaced from level k to level k

       lconv = .true.
        do 60 k=1,nconvtop
          rlevmass = dp(k)/ga
          summe = 0.
          do 61 kk=1,nconvtop
            fmassfrac(k,kk) = delt*fmass(k,kk)
            summe = summe + fmassfrac(k,kk)
61        continue
          fmassfrac(k,k)=fmassfrac(k,k) + rlevmass - summe
60      continue

200     continue

      end
