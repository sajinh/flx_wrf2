      subroutine calcpar(n,uuh,vvh,pvh)
C                        i  i   i   o
********************************************************************************
*                                                                              *
*     Computation of several boundary layer parameters needed for the          *
*     dispersion calculation and calculation of dry deposition velocities.     *
*     All parameters are calculated over the entire grid.                      *
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subroutine calcpar.           *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     21 May 1995                                                              *
*                                                                              *
* ------------------------------------------------------------------           *
*     Petra Seibert, Feb 2000:                                                 *
*     convection scheme:                                                       *
*     new variables in call to richardson                                      *
*                                                                              *
*     Changes, Bernd C. Krueger, Feb. 2001:
*        Variables tth and qvh (on eta coordinates) in common block
*                                                                              *
*     17 Oct 2005 - R. Easter - added ierr in call to richardson               *
*     18 Oct 2005 - J. Fast - limit ustar to < 5.0 m/s                         *
*     -- Oct 2005 - R. Easter - use xy_to_ll_wrf to get latitude               *
*             use pph for calculating zlev                                     *
*             pass level-2 pph directly to obukhov                             *
*     11 Nov 2005 - R. Easter - changed name of "xy to latlon" routine         *
*     15 Nov 2005 - R. Easter - pass pplev to richardson instead of akz,bkz    *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* n                  temporal index for meteorological fields (1 to 3)         *
*                                                                              *
* Constants:                                                                   *
*                                                                              *
*                                                                              *
* Functions:                                                                   *
* scalev             computation of ustar                                      *
* obukhov            computatio of Obukhov length                              *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      integer n,ix,jy,i,kz,lz,kzmin
      integer ierr
      integer ientry
      save ientry
      data ientry / 0 /
      real ttlev(nuvzmax),qvlev(nuvzmax),obukhov,scalev,ol,hmixplus
      real ulev(nuvzmax),vlev(nuvzmax),ew,rh,vd(maxspec),subsceff,ylat
      real xlon,dumx,dumy,dumxb,dumyb,pplev(nuvzmax)
      real altmin,tvold,pold,zold,pint,tv,zlev(nuvzmax),const
      real uuh(0:nxmax-1,0:nymax-1,nuvzmax)
      real vvh(0:nxmax-1,0:nymax-1,nuvzmax)
      real pvh(0:nxmax-1,0:nymax-1,nuvzmax)
      parameter(const=r_air/ga)


C Loop over entire grid
***********************
      ientry = ientry + 1

      do 10 jy=0,nymin1

        do 10 ix=0,nxmin1

C Set minimum height for tropopause
***********************************

c FLEXPART_WRF - use this routine to get lat,lon
c       ylat=ylat0+float(jy)*dy
        call xyindex_to_ll_wrf( 0, float(ix), float(jy), xlon, ylat )

c       if ( ((ix.eq.0) .or. (ix.eq.nxmin1) .or. 
c    &                       (ix.eq.nxmin1/2)) .and.
c    &       ((jy.eq.0) .or. (jy.eq.nymin1) .or. 
c    &                       (jy.eq.nymin1/2)) ) then
c           if (ientry .eq. 1) then
c               write(*,'(a,2i4,2f12.5)') 
c    &              'calcpar i,j, xlon,ylat', ix, jy, xlon, ylat
c               write(*,'(a, 8x,2f12.5)') 
c    &              '             dlon,dlat', 
c    &              (xlon-xlon2d(ix,jy)), (ylat-ylat2d(ix,jy))
c               call ll_to_xyindex_wrf(
c    &              xlon2d(ix,jy), ylat2d(ix,jy), dumx, dumy )
c               write(*,'(a, 8x,2f12.5)') 
c    &              '             dxkm,dykm', 
c    &              ((dumx-ix)*dx*1.0e-3), ((dumy-jy)*dy*1.0e-3) 
c
c               if ((ix .eq. 0) .and. (jy .eq. 0)) then
c                  dumxb = 2.33
c                  dumyb = 3.44
c                  call xyindex_to_ll_wrf( 0, dumxb, dumyb, dumx, dumy )
c                  call ll_to_xyindex_wrf( dumx, dumy, dumx, dumy )
c                  write(*,'(a,2f5.2,2f12.5)') 
c    &                'xi,yj,     dxkm,dykm', dumxb, dumyb,
c    &                ((dumx-dumxb)*dx*1.0e-3), ((dumy-dumyb)*dy*1.0e-3)
c                  dumxb = 4.55
c                  dumyb = 6.77
c                  call xyindex_to_ll_wrf( 0, dumxb, dumyb, dumx, dumy )
c                  call ll_to_xyindex_wrf( dumx, dumy, dumx, dumy )
c                  write(*,'(a,2f5.2,2f12.5)') 
c    &                'xi,yj,     dxkm,dykm', dumxb, dumyb,
c    &                ((dumx-dumxb)*dx*1.0e-3), ((dumy-dumyb)*dy*1.0e-3)
c               end if
c
c           end if
c       end if

        if ((ylat.ge.-20.).and.(ylat.le.20.)) then
          altmin = 5000.
        else
          if ((ylat.gt.20.).and.(ylat.lt.40.)) then
            altmin=2500.+(40.-ylat)*125.
          else if ((ylat.gt.-40.).and.(ylat.lt.-20.)) then
            altmin=2500.+(40.+ylat)*125.
          else
            altmin=2500.
          endif
        endif

C 1) Calculation of friction velocity
*************************************

          ustar(ix,jy,1,n)=scalev(ps(ix,jy,1,n),tt2(ix,jy,1,n),
     +    td2(ix,jy,1,n),surfstr(ix,jy,1,n))
          if (ustar(ix,jy,1,n).le.1.e-8) ustar(ix,jy,1,n)=1.e-8
c FLEXPART_WRF - limit ustar
          if (ustar(ix,jy,1,n).ge.5.0)   ustar(ix,jy,1,n)=5.0

C 2) Calculation of inverse Obukhov length scale
************************************************

c FLEXPART_WRF - pass k=2 pressure directly
c         ol=obukhov(ps(ix,jy,1,n),tt2(ix,jy,1,n),td2(ix,jy,1,n),
c    +    tth(ix,jy,2,n),ustar(ix,jy,1,n),sshf(ix,jy,1,n),akm,bkm)
          ol=obukhov(ps(ix,jy,1,n),tt2(ix,jy,1,n),td2(ix,jy,1,n),
     +    tth(ix,jy,2,n),ustar(ix,jy,1,n),sshf(ix,jy,1,n),
     +    pph(ix,jy,2,n) )
          if (ol.ne.0.) then
            oli(ix,jy,1,n)=1./ol
          else
            oli(ix,jy,1,n)=99999.
          endif


C 3) Calculation of convective velocity scale and mixing height
***************************************************************

          do 20 i=1,nuvz
            ulev(i) =uuh(ix,jy,i)
            vlev(i) =vvh(ix,jy,i)
            pplev(i)=pph(ix,jy,i,n)
            ttlev(i)=tth(ix,jy,i,n)
20          qvlev(i)=qvh(ix,jy,i,n)

c FLEXPART_WRF - use  & check ierr argument
c FLEXPART_WRF - pass pplev instead of akz,bkz
c         call richardson(ps(ix,jy,1,n),ustar(ix,jy,1,n),ttlev,qvlev,
c    +    ulev,vlev,nuvz,akz,bkz,sshf(ix,jy,1,n),tt2(ix,jy,1,n),
c    +    td2(ix,jy,1,n),hmix(ix,jy,1,n),wstar(ix,jy,1,n),hmixplus)
          call richardson(ps(ix,jy,1,n),ustar(ix,jy,1,n),ttlev,qvlev,
     +    ulev,vlev,nuvz,  pplev,sshf(ix,jy,1,n),tt2(ix,jy,1,n),
     +    td2(ix,jy,1,n),hmix(ix,jy,1,n),wstar(ix,jy,1,n),hmixplus,
     +    ierr,sfc_option )

          if (ierr .gt. 0) then
              write(*,9500) 'warning', ix, jy
          else if (ierr .lt. 0) then
              write(*,9500) 'failure', ix, jy
              stop
          end if
9500      format( 'calcpar - richardson ', a, ' - ix,jy=', 2i5 )

          if(lsubgrid.eq.1) then
            subsceff=min(excessoro(ix,jy),hmixplus)
          else
            subsceff=0
          endif
*
* CALCULATE HMIX EXCESS ACCORDING TO SUBGRIDSCALE VARIABILITY AND STABILITY
*
          hmix(ix,jy,1,n)=hmix(ix,jy,1,n)+subsceff
          hmix(ix,jy,1,n)=max(hmixmin,hmix(ix,jy,1,n)) ! set minimum PBL height
          hmix(ix,jy,1,n)=min(hmixmax,hmix(ix,jy,1,n)) ! set maximum PBL height

C 4) Calculation of dry deposition velocities
*********************************************

          if (DRYDEP) then
            z0(4)=0.016*ustar(ix,jy,1,n)*ustar(ix,jy,1,n)/ga
            z0(9)=0.016*ustar(ix,jy,1,n)*ustar(ix,jy,1,n)/ga

C Calculate relative humidity at surface
****************************************
            rh=ew(td2(ix,jy,1,n))/ew(tt2(ix,jy,1,n))

            call getvdep(n,ix,jy,ustar(ix,jy,1,n),
     +      tt2(ix,jy,1,n),ps(ix,jy,1,n),1./oli(ix,jy,1,n),
     +      ssr(ix,jy,1,n),rh,lsprec(ix,jy,1,n)+convprec(ix,jy,1,n),vd)

            do 30 i=1,nspec
30            vdep(ix,jy,i,n)=vd(i)

          endif

*******************************************************
C Calculate height of thermal tropopause (Hoinka, 1997)
*******************************************************

C 1) Calculate altitudes of ECMWF model levels
**********************************************

          tvold=tt2(ix,jy,1,n)*(1.+0.378*ew(td2(ix,jy,1,n))/
     +                                   ps(ix,jy,1,n))
          pold=ps(ix,jy,1,n)
          zold=0.
c FLEXPART_WRF - set zlev(1)
          zlev(1)=zold
          do 40 kz=2,nuvz
c FLEXPART_WRF - use pph for pressure
c           pint=akz(kz)+bkz(kz)*ps(ix,jy,1,n)  ! pressure on model layers
            pint=pph(ix,jy,kz,n)  ! pressure on model layers
            tv=tth(ix,jy,kz,n)*(1.+0.608*qvh(ix,jy,kz,n))

            if (abs(tv-tvold).gt.0.2) then
             zlev(kz)=zold+const*log(pold/pint)*(tv-tvold)/log(tv/tvold)
            else
              zlev(kz)=zold+const*log(pold/pint)*tv
            endif
            tvold=tv
            pold=pint
            zold=zlev(kz)
40          continue

C 2) Define a minimum level kzmin, from which upward the tropopause is
C    searched for. This is to avoid inversions in the lower troposphere
C    to be identified as the tropopause
*************************************************************************

          do 44 kz=1,nuvz
            if (zlev(kz).ge.altmin) then
              kzmin=kz
              goto 45
            endif
44          continue
45        continue

C 3) Search for first stable layer above minimum height that fulfills the
C    thermal tropopause criterion
*************************************************************************

          do 50 kz=kzmin,nuvz
            do 60 lz=kz+1,nuvz
              if ((zlev(lz)-zlev(kz)).gt.2000.) then
                if (((tth(ix,jy,kz,n)-tth(ix,jy,lz,n))/
     +          (zlev(lz)-zlev(kz))).lt.0.002) then
                  tropopause(ix,jy,1,n)=zlev(kz)
                  goto 51
                endif
                goto 50
              endif
60            continue
50          continue
51        continue


10        continue


C Calculation of potential vorticity on 3-d grid, if plume trajectory mode is used
**********************************************************************************

      if ((iout.eq.4).or.(iout.eq.5).or.(mdomainfill.eq.2)) then
        call calcpv(n,uuh,vvh,pvh)
      endif


      end
