      subroutine redist (ipart,ktop,ipconv)
c**************************************************************************
C Do the redistribution of particles due to convection
C This subroutine is called for each particle which is assigned
C a new vertical position randomly, based on the convective redistribution
C matrix
c**************************************************************************

C Petra Seibert, Feb 2001, Apr 2001, May 2001, Jan 2002, Nov 2002 and
C Andreas Frank, Nov 2002

c Caroline Forster:  November 2004 - February 2005

      include 'includepar'
      include 'includecom'
      include 'includeconv'

      real const
      parameter (const=r_air/ga) 
      integer ipart, ktop,ipconv
      integer iseed, k, kz, levnew, levold
      real uvzlev(nuvzmax)
      real wsub(nuvzmax)
      real totlevmass, wsubpart 
      real temp_levold,temp_levold1
      real sub_levold,sub_levold1
      real pint, pold, rn, tv, tvold, dlevfrac
      real ew,ran3, ztold,fraction
      real tv1, tv2, dlogp, dz, dz1, dz2
      save iseed, uvzlev
      data iseed /-88/

C ipart   ... number of particle to be treated

      ipconv=1

c determine height of the eta half-levels (uvzlev)
c do that only once for each grid column
c i.e. when ktop.eq.1
***************************************************************

      if (ktop .le. 1) then

        tvold=tt2conv*(1.+0.378*ew(td2conv)/psconv)
        pold=psconv
        uvzlev(1)=0.

        pint = phconv(2)
C       determine next virtual temperatures
        tv1 = tconv(1)*(1.+0.608*qconv(1))
        tv2 = tconv(2)*(1.+0.608*qconv(2))
C       interpolate virtual temperature to half-level
        tv = tv1 + (tv2-tv1)*(pconv(1)-phconv(2))/(pconv(1)-pconv(2))
        if (abs(tv-tvold).gt.0.2) then
          uvzlev(2) = uvzlev(1) + 
     &                  const*log(pold/pint)*
     &                 (tv-tvold)/log(tv/tvold)
        else
          uvzlev(2) = uvzlev(1)+
     &          const*log(pold/pint)*tv
        endif
        tvold=tv
        tv1=tv2
        pold=pint

C integrate profile (calculation of height agl of eta layers) as required
        do 10 kz = 3, nconvtop+1
C         note that variables defined in calcmatrix.f (pconv,tconv,qconv)
C         start at the first real ECMWF model level whereas kz and
C         thus uvzlev(kz) starts at the surface. uvzlev is defined at the
C         half-levels (between the tconv, qconv etc. values !)
C         Thus, uvzlev(kz) is the lower boundary of the tconv(kz) cell.
          pint = phconv(kz)
C         determine next virtual temperatures
          tv2 = tconv(kz)*(1.+0.608*qconv(kz))
C         interpolate virtual temperature to half-level
          tv = tv1 + (tv2-tv1)*(pconv(kz-1)-phconv(kz))/
     &         (pconv(kz-1)-pconv(kz))
          if (abs(tv-tvold).gt.0.2) then
            uvzlev(kz) = uvzlev(kz-1) + 
     &                  const*log(pold/pint)*
     &                 (tv-tvold)/log(tv/tvold)
          else
            uvzlev(kz) = uvzlev(kz-1)+
     &          const*log(pold/pint)*tv
          endif
          tvold=tv
          tv1=tv2
          pold=pint

10      continue

        ktop = 2

      endif ! (if ktop .le. 1) then

C  determine vertical grid position of particle in the eta system
C****************************************************************

      ztold = ztra1(abs(ipart))
C find old particle grid position
      do 20 kz = 2, nconvtop
        if (uvzlev(kz) .ge. ztold ) then
          levold = kz-1
          goto 30
        endif
20    continue

C Particle is above the potentially convective domain. Skip it.
      goto 90

30    continue

c now redistribute particles
c****************************

C  Choose a random number and find corresponding level of destination
C  Random numbers to be evenly distributed in [0,1]

      rn = ran3(iseed)

c initialize levnew

      levnew = levold

      fraction = 0.
      totlevmass=dp(levold)/ga
      do 35 k = 1,nconvtop
c for backward runs use the transposed matrix
       if (ldirect.eq.1) then
         fraction=fraction+fmassfrac(levold,k)
     +       /totlevmass
       else
         fraction=fraction+fmassfrac(k,levold)
     +       /totlevmass
       endif
       if (rn.le.fraction) then
         levnew=k
c avoid division by zero or a too small number
c if division by zero or a too small number happens the
c particle is assigned to the center of the grid cell
         if (fraction.gt.1.e-20) then
          if (ldirect.eq.1) then
            dlevfrac = (fraction-rn) / fmassfrac(levold,k) * totlevmass
          else
            dlevfrac = (fraction-rn) / fmassfrac(k,levold) * totlevmass 
          endif
         else
           dlevfrac = 0.5
         endif
         goto 40
       endif
 35   continue

40    continue

C now assign new position to particle

      if (levnew.le.nconvtop) then
       if (levnew.eq.levold) then
          ztra1(abs(ipart)) = ztold
       else
        dlogp = (1.-dlevfrac)*
     +        (log(phconv(levnew+1))-log(phconv(levnew)))
        pint = log(phconv(levnew))+dlogp
        dz1 = pint - log(phconv(levnew))
        dz2 = log(phconv(levnew+1)) - pint
        dz = dz1 + dz2
        ztra1(abs(ipart)) = (uvzlev(levnew)*dz2+uvzlev(levnew+1)*dz1)/dz
         if (ztra1(abs(ipart)).lt.0.) 
     +       ztra1(abs(ipart))=-1.*ztra1(abs(ipart))
         if (ipconv.gt.0) ipconv=-1
       endif
      endif

c displace particle according to compensating subsidence
c this is done to those particles, that were not redistributed 
c by the matrix
***************************************************************

      if (levnew.le.nconvtop.and.levnew.eq.levold) then

      ztold = ztra1(abs(ipart))

c determine compensating vertical velocity at the levels
c above and below the particel position
c increase compensating subsidence by the fraction that
c is displaced by convection to this level

        if (levold.gt.1) then
         temp_levold = tconv(levold-1) + 
     &           (tconv(levold)-tconv(levold-1))
     &           *(pconv(levold-1)-phconv(levold))/
     &         (pconv(levold-1)-pconv(levold))
         sub_levold = sub(levold)/(1.-sub(levold)/dp(levold)*ga)
         wsub(levold)=-1.*sub_levold*r_air*temp_levold/(phconv(levold))
        else
         wsub(levold)=0.
        endif

         temp_levold1 = tconv(levold) + 
     &           (tconv(levold+1)-tconv(levold))
     &           *(pconv(levold)-phconv(levold+1))/
     &         (pconv(levold)-pconv(levold+1))
         sub_levold1 = sub(levold+1)/(1.-sub(levold+1)/dp(levold+1)*ga)
         wsub(levold+1)=-1.*sub_levold1*r_air*temp_levold1/
     & (phconv(levold+1))

c interpolate wsub to the vertical particle position

      dz1 = ztold - uvzlev(levold)
      dz2 = uvzlev(levold+1) - ztold
      dz = dz1 + dz2

      wsubpart = (dz2*wsub(levold)+dz1*wsub(levold+1))/dz
      ztra1(abs(ipart)) = ztold+wsubpart*float(lsynctime)
      if (ztra1(abs(ipart)).lt.0.) then
         ztra1(abs(ipart))=-1.*ztra1(abs(ipart))
      endif

      endif      !(levnew.le.nconvtop.and.levnew.eq.levold)

C Maximum altitude .5 meter below uppermost model level
********************************************************

 90   continue

      if (ztra1(abs(ipart)) .gt. height(nz)-0.5)
     &    ztra1(abs(ipart)) = height(nz)-0.5

      end
