      subroutine calcpar_nests(n,uuhn,vvhn,pvhn)
C                              i  i    i    o
********************************************************************************
*                                                                              *
*     Computation of several boundary layer parameters needed for the          *
*     dispersion calculation and calculation of dry deposition velocities.     *
*     All parameters are calculated over the entire grid.                      *
*     This routine is similar to calcpar, but is used for the nested grids.    *
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subroutine calcpar.           *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     8 February 1999                                                          *
*                                                                              *
* ------------------------------------------------------------------           *
*     Petra Seibert, Feb 2000:                                                 *
*     convection scheme:                                                       *
*     new variables in call to richardson                                      *
*                                                                              *
*     Changes, Bernd C. Krueger, Feb. 2001:                                    *
*        Variables tth and qvh (on eta coordinates) in common block            *
*                                                                              *
*     14 Nov 2005 - R. Easter -                                                *
*          use xyindex_to_ll_wrf to get latitude                               *
*          limit ustar to < 5.0 m/s                                            *
*          added ierr in call to richardson                                    *
*          use pph for calculating zlev                                        *
*          pass level-2 pph directly to obukhov                                *
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

      integer n,ix,jy,i,l,kz,lz,kzmin
      integer ierr
      integer ientry
      save ientry
      data ientry / 0 /
      real ttlev(nuvzmax),qvlev(nuvzmax),obukhov,scalev,ol,hmixplus
      real ulev(nuvzmax),vlev(nuvzmax),ew,rh,vd(maxspec),subsceff,ylat
      real xlon,pplev(nuvzmax)
      real altmin,tvold,pold,zold,pint,tv,zlev(nuvzmax),const
      real uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
      real vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
      real pvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
      parameter(const=r_air/ga)


C Loop over all nests
*********************
      ientry = ientry + 1

      do 80 l=1,numbnests

C Loop over entire grid
***********************

      do 10 jy=0,nyn(l)-1

        do 10 ix=0,nxn(l)-1

C Set minimum height for tropopause
***********************************

c FLEXPART_WRF - use this routine to get lat,lon
c       ylat=ylat0n(l)+float(jy)*dyn(l)
        call xyindex_to_ll_wrf( l, float(ix), float(jy), xlon, ylat )

c       if ( ((ix.eq.0) .or. (ix.eq.(nxn(l)-1)) .or. 
c    &                       (ix.eq.(nxn(l)-1)/2)) .and.
c    &       ((jy.eq.0) .or. (jy.eq.(nyn(l)-1)) .or. 
c    &                       (jy.eq.(nyn(l)-1)/2)) ) then
c           if (ientry .eq. 1) then
c               write(*,'(a,3i4,2f12.5)') 
c    &              'calcpar_nests l,i,j, xlon,ylat', 
c    &              l, ix, jy, xlon, ylat
c               write(*,'(a,12x,2f12.5)') 
c    &              '                     dlon,dlat', 
c    &              (xlon-xlon2dn(ix,jy,l)), (ylat-ylat2dn(ix,jy,l))
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

          ustarn(ix,jy,1,n,l)=scalev(psn(ix,jy,1,n,l),tt2n(ix,jy,1,n,l),
     +    td2n(ix,jy,1,n,l),surfstrn(ix,jy,1,n,l))
c FLEXPART_WRF - limit ustar
          if (ustarn(ix,jy,1,n,l).le.1.e-8) ustarn(ix,jy,1,n,l)=1.e-8
          if (ustarn(ix,jy,1,n,l).ge.5.0)   ustarn(ix,jy,1,n,l)=5.0

C 2) Calculation of inverse Obukhov length scale
************************************************

c FLEXPART_WRF - pass k=2 pressure directly
c         ol=obukhov(psn(ix,jy,1,n,l),tt2n(ix,jy,1,n,l),
c    +    td2n(ix,jy,1,n,l),tthn(ix,jy,2,n,l),ustarn(ix,jy,1,n,l),
c    +    sshfn(ix,jy,1,n,l),akm,bkm)
          ol=obukhov(psn(ix,jy,1,n,l),tt2n(ix,jy,1,n,l),
     +    td2n(ix,jy,1,n,l),tthn(ix,jy,2,n,l),ustarn(ix,jy,1,n,l),
     +    sshfn(ix,jy,1,n,l),pphn(ix,jy,2,n,l) )
          if (ol.ne.0.) then
            olin(ix,jy,1,n,l)=1./ol
          else
            olin(ix,jy,1,n,l)=99999.
          endif


C 3) Calculation of convective velocity scale and mixing height
***************************************************************

          do 20 i=1,nuvz
            ulev(i) =uuhn(ix,jy,i,l)
            vlev(i) =vvhn(ix,jy,i,l)
            pplev(i)=pphn(ix,jy,i,n,l)
            ttlev(i)=tthn(ix,jy,i,n,l)
20          qvlev(i)=qvhn(ix,jy,i,n,l)

c FLEXPART_WRF - use  & check ierr argument
c FLEXPART_WRF - pass pplev instead of akz,bkz
c         call richardson(psn(ix,jy,1,n,l),ustarn(ix,jy,1,n,l),ttlev,
c    +    qvlev,ulev,vlev,nuvz,akz,bkz,sshfn(ix,jy,1,n,l),
c    +    tt2n(ix,jy,1,n,l),td2n(ix,jy,1,n,l),hmixn(ix,jy,1,n,l),
c    +    wstarn(ix,jy,1,n,l),hmixplus)
          call richardson(psn(ix,jy,1,n,l),ustarn(ix,jy,1,n,l),ttlev,
     +    qvlev,ulev,vlev,nuvz,  pplev,sshfn(ix,jy,1,n,l),
     +    tt2n(ix,jy,1,n,l),td2n(ix,jy,1,n,l),hmixn(ix,jy,1,n,l),
     +    wstarn(ix,jy,1,n,l),hmixplus,ierr,sfc_option)

          if (ierr .gt. 0) then
              write(*,9500) 'warning', l, ix, jy
          else if (ierr .lt. 0) then
              write(*,9500) 'failure', l, ix, jy
              stop
          end if
9500      format( 'calcpar_nests - richardson ', a, ' - l,ix,jy=', 3i5 )

          if(lsubgrid.eq.1) then
            subsceff=min(excessoron(ix,jy,l),hmixplus)
          else
            subsceff=0
          endif
*
* CALCULATE HMIX EXCESS ACCORDING TO SUBGRIDSCALE VARIABILITY AND STABILITY
*
          hmixn(ix,jy,1,n,l)=hmixn(ix,jy,1,n,l)+subsceff
          hmixn(ix,jy,1,n,l)=max(hmixmin,hmixn(ix,jy,1,n,l)) ! minim PBL height
          hmixn(ix,jy,1,n,l)=min(hmixmax,hmixn(ix,jy,1,n,l)) ! maxim PBL height


C 4) Calculation of dry deposition velocities
*********************************************

          if (DRYDEP) then
            z0(4)=0.016*ustarn(ix,jy,1,n,l)*ustarn(ix,jy,1,n,l)/ga
            z0(9)=0.016*ustarn(ix,jy,1,n,l)*ustarn(ix,jy,1,n,l)/ga

C Calculate relative humidity at surface
****************************************
            rh=ew(td2n(ix,jy,1,n,l))/ew(tt2n(ix,jy,1,n,l))

            call getvdep(n,ix,jy,ustarn(ix,jy,1,n,l),
     +      tt2n(ix,jy,1,n,l),psn(ix,jy,1,n,l),1./olin(ix,jy,1,n,l),
     +      ssrn(ix,jy,1,n,l),rh,lsprecn(ix,jy,1,n,l)+
     +      convprecn(ix,jy,1,n,l),vd)

            do 30 i=1,nspec
30            vdepn(ix,jy,i,n,l)=vd(i)

          endif

*******************************************************
C Calculate height of thermal tropopause (Hoinka, 1997)
*******************************************************

C 1) Calculate altitudes of ECMWF model levels
**********************************************

          tvold=tt2n(ix,jy,1,n,l)*(1.+0.378*ew(td2n(ix,jy,1,n,l))/
     +                                   psn(ix,jy,1,n,l))
          pold=psn(ix,jy,1,n,l)
          zold=0.
c FLEXPART_WRF - set zlev(1)
          zlev(1)=zold
          do 40 kz=2,nuvz
c FLEXPART_WRF - use pph for pressure
c           pint=akz(kz)+bkz(kz)*psn(ix,jy,1,n,l)  ! pressure on model layers
            pint=pphn(ix,jy,kz,n,l)  ! pressure on model layers
            tv=tthn(ix,jy,kz,n,l)*(1.+0.608*qvhn(ix,jy,kz,n,l))

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
                if (((tthn(ix,jy,kz,n,l)-tthn(ix,jy,lz,n,l))/
     +          (zlev(lz)-zlev(kz))).lt.0.002) then
                  tropopausen(ix,jy,1,n,l)=zlev(kz)
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

        if ((iout.eq.4).or.(iout.eq.5)) then
          call calcpv_nests(l,n,uuhn,vvhn,pvhn)
        endif

80      continue


      return
      end
