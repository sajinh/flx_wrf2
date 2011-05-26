      subroutine verttransform(n,uuh,vvh,wwh,pvh)
C                              i  i   i   i   i
********************************************************************************
*                                                                              *
*     This subroutine transforms temperature, dew point temperature and        *
*     wind components from eta to meter coordinates.                           *
*     The vertical wind component is transformed from Pa/s to m/s using        *
*     the conversion factor pinmconv.                                          *
*     In addition, this routine calculates vertical density gradients          *
*     needed for the parameterization of the turbulent velocities.             *
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subroutine assignland.        *
*            The computational grid is the WRF x-y grid rather than lat-lon.   *
*                                                                              *
*     Author: A. Stohl, G. Wotawa                                              *
*                                                                              *
*     12 August 1996                                                           *
*     Update: 16 January 1998                                                  *
*                                                                              *
*     Major update: 17 February 1999                                           *
*     by G. Wotawa                                                             *
*                                                                              *
*     - Vertical levels for u, v and w are put together                        *
*     - Slope correction for vertical velocity: Modification of calculation    *
*       procedure                                                              *
*                                                                              *
*     Changes, Bernd C. Krueger, Feb. 2001:                                    *
*        Variables tth and qvh (on eta coordinates) from common block          *
*                                                                              *
*     Oct-Nov 2005 - R. Easter - conversion to wrf                             *
*     17 Nov 2005 - R. Easter - terrain correction applied to ww.  There are   *
*            now 3 options, controlled by "method_w_terrain_correction"        *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* nx,ny,nz                        field dimensions in x,y and z direction      *
* uu(0:nxmax,0:nymax,nzmax,2)     wind components in x-direction [m/s]         *
* vv(0:nxmax,0:nymax,nzmax,2)     wind components in y-direction [m/s]         *
* ww(0:nxmax,0:nymax,nzmax,2)     wind components in z-direction [deltaeta/s]  *
* tt(0:nxmax,0:nymax,nzmax,2)     temperature [K]                              *
* pv(0:nxmax,0:nymax,nzmax,2)     potential voriticity (pvu)                   *
* ps(0:nxmax,0:nymax,2)           surface pressure [Pa]                        *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      integer ix,jy,kz,iz,n,kmin,kl,klp,ix1,jy1,ixp,jyp,ixm,jym
      save ixm,jym
      integer method_z_compute
      real uvzlev(nuvzmax),rhoh(nuvzmax),pinmconv(nzmax)
      real ew,pint,tv,tvold,pold,const,dz1,dz2,dz,ui,vi
      real xlon,ylat,xlonr,dzdx,dzdy
      real dzdx1,dzdx2,dzdy1,dzdy2
      real uuaux,vvaux,uupolaux,vvpolaux,ddpol,ffpol,wdummy
      real uuh(0:nxmax-1,0:nymax-1,nuvzmax)
      real vvh(0:nxmax-1,0:nymax-1,nuvzmax)
      real pvh(0:nxmax-1,0:nymax-1,nuvzmax)
      real wwh(0:nxmax-1,0:nymax-1,nwzmax)
      real wzlev(nwzmax),uvwzlev(0:nxmax-1,0:nymax-1,nzmax)
      real wwh_svaa(nwzmax), wtc_stat(4,nzmax)
      parameter(const=r_air/ga)

      logical init
      save init
      data init/.true./



c set method_w_terrain_correction  & method_z_compute
      method_w_terrain_correction = 11
      method_w_terrain_correction = 20
      method_z_compute = 10

      do iz = 1, nz
      do ix = 1, 4
          wtc_stat(ix,iz) = 0.0
      end do
      end do


**************************************************************************
* If verttransform is called the first time, initialize heights of the   *
* z levels in meter. The heights are the heights of model levels, where  *
* u,v,T and qv are given, and of the interfaces, where w is given. So,   *
* the vertical resolution in the z system is doubled. As reference point,*
* the lower left corner of the grid is used.                             *
* Unlike in the eta system, no difference between heights for u,v and    *
* heights for w exists.                                                  *
**************************************************************************

      if (init) then

C Search for a point with high surface pressure (i.e. not above significant topography)
C Then, use this point to construct a reference z profile, to be used at all times
C
C FLEXPART_WRF - use grid point with highest surface pressure
***************************************************************************************

        pint = -1.0
        ixm = -999888777
        jym = -999888777
        do 4 jy=0,nymin1
          do 4 ix=0,nxmin1
c           if (ps(ix,jy,1,n).gt.100000.) then
            if (ps(ix,jy,1,n).gt.pint) then
              pint = ps(ix,jy,1,n)
              ixm=ix
              jym=jy
c             goto 3
            endif
4           continue
3       continue
c       write(*,'(/a,2i4,1pe11.2)') 
c    &          'verttransform -- ixm,jym,ps() =', ixm, jym, pint


        tvold=tt2(ixm,jym,1,n)*(1.+0.378*ew(td2(ixm,jym,1,n))/
     +  ps(ixm,jym,1,n))
        pold=ps(ixm,jym,1,n)
        height(1)=0.

        do 5 kz=2,nuvz
c use pressure from wrf met file
c         pint=akz(kz)+bkz(kz)*ps(ixm,jym,1,n)
          pint=pph(ixm,jym,kz,n)
          tv=tth(ixm,jym,kz,n)*(1.+0.608*qvh(ixm,jym,kz,n))


C NOTE: In FLEXPART versions up to 4.0, the number of model levels was doubled
C upon the transformation to z levels. In order to save computer memory, this is
C not done anymore in the standard version. However, this option can still be
C switched on by replacing the following lines with those below, that are
C currently commented out.
C Note that two more changes are necessary in this subroutine below.
C One change is also necessary in gridcheck.f, and another one in verttransform_nests.
**************************************************************************************

          if (abs(tv-tvold).gt.0.2) then
            height(kz)=
     +      height(kz-1)+const*log(pold/pint)*
     +      (tv-tvold)/log(tv/tvold)
          else
            height(kz)=height(kz-1)+
     +      const*log(pold/pint)*tv
          endif

C 
C *** NOTE -- the doubled vertical resolution has not been tested in FLEXPART_WRF
C 
C Switch on following lines to use doubled vertical resolution
**************************************************************
c         if (abs(tv-tvold).gt.0.2) then
c           height((kz-1)*2)=
c    +      height(max((kz-2)*2,1))+const*log(pold/pint)*
c    +      (tv-tvold)/log(tv/tvold)
c         else
c           height((kz-1)*2)=height(max((kz-2)*2,1))+
c    +      const*log(pold/pint)*tv
c         endif
C End doubled vertical resolution
 
c FLEXPART_WRF - get height from zzh
          if (method_z_compute .eq. 10) then
             if ((add_sfc_level .eq. 1) .and. (kz .eq. 2)) then
                height(kz) = 0.5*(zzh(ixm,jym,   3,n)+zzh(ixm,jym, 1,n))
     &                     - zzh(ixm,jym,1,n)
             else
                height(kz) = 0.5*(zzh(ixm,jym,kz+1,n)+zzh(ixm,jym,kz,n))
     &                     - zzh(ixm,jym,1,n)
             end if
          end if

          tvold=tv
5         pold=pint

C 
C *** NOTE -- the doubled vertical resolution has not been tested in FLEXPART_WRF
C 
C Switch on following lines to use doubled vertical resolution
**************************************************************
c       do 7 kz=3,nz-1,2
c         height(kz)=0.5*(height(kz-1)+height(kz+1))
c       height(nz)=height(nz-1)+height(nz-1)-height(nz-2)
C End doubled vertical resolution


C Determine highest levels that can be within PBL
*************************************************

        do 8 kz=1,nz
          if (height(kz).gt.hmixmax) then
            nmixz=kz
            goto 9
          endif
8         continue
9       continue

C Do not repeat initialization of the Cartesian z grid
******************************************************

        init=.false.

      endif


C Loop over the whole grid
**************************

      do 10 jy=0,nymin1
        do 10 ix=0,nxmin1
          tvold=tt2(ix,jy,1,n)*(1.+0.378*ew(td2(ix,jy,1,n))/
     &                                   ps(ix,jy,1,n))
          pold=ps(ix,jy,1,n)
          uvzlev(1)=0.
          wzlev(1)=0.
          rhoh(1)=pold/(r_air*tvold)


C Compute heights of eta levels
*******************************

          do 20 kz=2,nuvz
c use pressure from wrf met file
c           pint=akz(kz)+bkz(kz)*ps(ix,jy,1,n)
            pint=pph(ix,jy,kz,n)
            tv=tth(ix,jy,kz,n)*(1.+0.608*qvh(ix,jy,kz,n))
            rhoh(kz)=pint/(r_air*tv)

            if (abs(tv-tvold).gt.0.2) then
              uvzlev(kz)=uvzlev(kz-1)+const*log(pold/pint)*
     +        (tv-tvold)/log(tv/tvold)
            else
              uvzlev(kz)=uvzlev(kz-1)+const*log(pold/pint)*tv
            endif
           
            tvold=tv
20          pold=pint

          do 21 kz=2,nwz-1
21          wzlev(kz)=(uvzlev(kz+1)+uvzlev(kz))/2.
          wzlev(nwz)=wzlev(nwz-1)+
     +    uvzlev(nuvz)-uvzlev(nuvz-1)

c FLEXPART_WRF - get uvzlev & wzlev from zzh
          if (method_z_compute .eq. 10) then
            do kz = 2, nuvz
              if ((add_sfc_level .eq. 1) .and. (kz .eq. 2)) then
                uvzlev(kz) = 0.5*(zzh(ix,jy,   3,n) + zzh(ix,jy, 1,n))
     &                     - zzh(ix,jy,1,n)
              else
                uvzlev(kz) = 0.5*(zzh(ix,jy,kz+1,n) + zzh(ix,jy,kz,n))
     &                     - zzh(ix,jy,1,n)
              end if
            end do
            do kz = 2, nwz
              wzlev(kz) = zzh(ix,jy,kz+add_sfc_level,n) 
     &                  - zzh(ix,jy,1,n)
            end do
          end if

          uvwzlev(ix,jy,1)=0.0
          do 22 kz=2,nuvz
22          uvwzlev(ix,jy,kz)=uvzlev(kz)

c     if ((ix .eq. ixm) .and. (jy .eq. jym)) then
c        write(*,'(/a)') 
c    &     'kz, height, uvzlev, wzlev, zzh-zzh(1) at ixm,jym  (in km)'
c        write(*,'(i3,4f8.3)') (kz, height(kz)*1.0e-3, 
c    &     uvzlev(kz)*1.0e-3, wzlev(kz)*1.0e-3,
c    &     (zzh(ix,jy,kz,n)-zzh(ix,jy,1,n))*1.0e-3, kz=nz,1,-1)
c        ixm = -9
c     end if

C Switch on following lines to use doubled vertical resolution
C Switch off the three lines above.
**************************************************************
c22          uvwzlev(ix,jy,(kz-1)*2)=uvzlev(kz)
c          do 23 kz=2,nwz
c23          uvwzlev(ix,jy,(kz-1)*2+1)=wzlev(kz)
C End doubled vertical resolution

C pinmconv=(h2-h1)/(p2-p1)
C
C in flexpart_ecmwf, pinmconv is used to convert etadot to w
C in FLEXPART_WRF, vertical velocity is already m/s, so pinmconv=1.0

c         pinmconv(1)=(uvwzlev(ix,jy,2)-uvwzlev(ix,jy,1))/
c    +    ((aknew(2)+bknew(2)*ps(ix,jy,1,n))-
c    +    (aknew(1)+bknew(1)*ps(ix,jy,1,n)))
          pinmconv(1)=1.0
          do 24 kz=2,nz-1
c           pinmconv(kz)=(uvwzlev(ix,jy,kz+1)-uvwzlev(ix,jy,kz-1))/
c    +      ((aknew(kz+1)+bknew(kz+1)*ps(ix,jy,1,n))-
c    +      (aknew(kz-1)+bknew(kz-1)*ps(ix,jy,1,n)))
            pinmconv(kz)=1.0
24        continue
c         pinmconv(nz)=(uvwzlev(ix,jy,nz)-uvwzlev(ix,jy,nz-1))/
c    +    ((aknew(nz)+bknew(nz)*ps(ix,jy,1,n))-
c    +    (aknew(nz-1)+bknew(nz-1)*ps(ix,jy,1,n)))
          pinmconv(nz)=1.0

C Levels, where u,v,t and q are given
*************************************

          uu(ix,jy,1,n)=uuh(ix,jy,1)
          vv(ix,jy,1,n)=vvh(ix,jy,1)
          tt(ix,jy,1,n)=tth(ix,jy,1,n)
          qv(ix,jy,1,n)=qvh(ix,jy,1,n)
          pv(ix,jy,1,n)=pvh(ix,jy,1)
          rho(ix,jy,1,n)=rhoh(1)
          uu(ix,jy,nz,n)=uuh(ix,jy,nuvz)
          vv(ix,jy,nz,n)=vvh(ix,jy,nuvz)
          tt(ix,jy,nz,n)=tth(ix,jy,nuvz,n)
          qv(ix,jy,nz,n)=qvh(ix,jy,nuvz,n)
          pv(ix,jy,nz,n)=pvh(ix,jy,nuvz)
          rho(ix,jy,nz,n)=rhoh(nuvz)
          kmin=2
          do 30 iz=2,nz-1
            do 35 kz=kmin,nuvz
              if(height(iz).gt.uvzlev(nuvz)) then
                uu(ix,jy,iz,n)=uu(ix,jy,nz,n)
                vv(ix,jy,iz,n)=vv(ix,jy,nz,n)
                tt(ix,jy,iz,n)=tt(ix,jy,nz,n)
                qv(ix,jy,iz,n)=qv(ix,jy,nz,n)
                pv(ix,jy,iz,n)=pv(ix,jy,nz,n)
                rho(ix,jy,iz,n)=rho(ix,jy,nz,n)
                goto 30
              endif
              if ((height(iz).gt.uvzlev(kz-1)).and.
     +            (height(iz).le.uvzlev(kz))) then
               dz1=height(iz)-uvzlev(kz-1)
               dz2=uvzlev(kz)-height(iz)
               dz=dz1+dz2
               uu(ix,jy,iz,n)=(uuh(ix,jy,kz-1)*dz2+uuh(ix,jy,kz)*dz1)/dz
               vv(ix,jy,iz,n)=(vvh(ix,jy,kz-1)*dz2+vvh(ix,jy,kz)*dz1)/dz
               tt(ix,jy,iz,n)=(tth(ix,jy,kz-1,n)*dz2
     $              +tth(ix,jy,kz,n)*dz1)/dz
               qv(ix,jy,iz,n)=(qvh(ix,jy,kz-1,n)*dz2
     $              +qvh(ix,jy,kz,n)*dz1)/dz
               pv(ix,jy,iz,n)=(pvh(ix,jy,kz-1)*dz2+pvh(ix,jy,kz)*dz1)/dz
               rho(ix,jy,iz,n)=(rhoh(kz-1)*dz2+rhoh(kz)*dz1)/dz
               kmin=kz
               goto 30
              endif
35            continue
30          continue


C Levels, where w is given
**************************

          ww(ix,jy,1,n)=wwh(ix,jy,1)*pinmconv(1)
          ww(ix,jy,nz,n)=wwh(ix,jy,nwz)*pinmconv(nz)
          kmin=2
          do 40 iz=2,nz
            do 45 kz=kmin,nwz
              if ((height(iz).gt.wzlev(kz-1)).and.
     +            (height(iz).le.wzlev(kz))) then
               dz1=height(iz)-wzlev(kz-1)
               dz2=wzlev(kz)-height(iz)
               dz=dz1+dz2
               ww(ix,jy,iz,n)=(wwh(ix,jy,kz-1)*dz2+wwh(ix,jy,kz)*dz1)/dz
     +         *pinmconv(iz)
               kmin=kz
               goto 40
              endif
45            continue
40          continue

          if (method_w_terrain_correction .eq. 20) then
c apply w correction assuming that the WRF w is "absolute w";
c apply it here to wwh; set wwh=0 at iz=1
             do iz = 1, nz
                wtc_stat(1,iz) = wtc_stat(1,iz) + ww(ix,jy,iz,n)
                wtc_stat(2,iz) = wtc_stat(2,iz) + abs(ww(ix,jy,iz,n))
             end do

             if ((ix.eq.0) .and. (jy.eq.0)) write(*,*) 
     &            'verttransform doing method_w_terrain_correction =', 
     &            method_w_terrain_correction
             ix1 = max( ix-1, 0 )
             jy1 = max( jy-1, 0 )
             ixp = min( ix+1, nx-1 )
             jyp = min( jy+1, ny-1 )
             dzdx=(oro(ixp,jy) - oro(ix1,jy))/(dx*(ixp-ix1))
             dzdy=(oro(ix,jyp) - oro(ix,jy1))/(dy*(jyp-jy1))

             do kz = 1, nwz
                wwh_svaa(kz) = wwh(ix,jy,kz)
                wwh(ix,jy,kz) = wwh(ix,jy,kz)
     &               - (uuh(ix,jy,kz)*dzdx + vvh(ix,jy,kz)*dzdy)  
                if (kz .eq. 1) wwh(ix,jy,kz) = 0.0
             end do

             ww(ix,jy,1,n)=wwh(ix,jy,1)*pinmconv(1)
             ww(ix,jy,nz,n)=wwh(ix,jy,nwz)*pinmconv(nz)
             kmin=2
             do 4000 iz=2,nz
               do 4500 kz=kmin,nwz
                 if ((height(iz).gt.wzlev(kz-1)).and.
     +               (height(iz).le.wzlev(kz))) then
                  dz1=height(iz)-wzlev(kz-1)
                  dz2=wzlev(kz)-height(iz)
                  dz=dz1+dz2
                  ww(ix,jy,iz,n)=(wwh(ix,jy,kz-1)*dz2+wwh(ix,jy,kz)*dz1)
     +              /dz*pinmconv(iz)
                  kmin=kz
                  goto 4000
                 endif
4500           continue
4000         continue

             do kz = 1, nwz
                wwh(ix,jy,kz) = wwh_svaa(kz)
             end do

             do iz = 1, nz
                wtc_stat(3,iz) = wtc_stat(3,iz) + ww(ix,jy,iz,n)
                wtc_stat(4,iz) = wtc_stat(4,iz) + abs(ww(ix,jy,iz,n))
             end do
          end if

C Compute density gradients at intermediate levels
**************************************************

          drhodz(ix,jy,1,n)=(rho(ix,jy,2,n)-rho(ix,jy,1,n))/
     +      (height(2)-height(1))
          do 50 kz=2,nz-1
50          drhodz(ix,jy,kz,n)=(rho(ix,jy,kz+1,n)-rho(ix,jy,kz-1,n))/
     +      (height(kz+1)-height(kz-1))
          drhodz(ix,jy,nz,n)=drhodz(ix,jy,nz-1,n)

10        continue


*****************************************************************
C Compute slope of eta levels in windward direction and resulting
C vertical wind correction
c
c The ECMWF model uses a hybrid-pressure vertical coordinate, "eta"
c    The "eta" coordinate transitions from terrain-following near
c    the surface to constant pressure in the stratosphere.
c    The vertical velocities in the ECMWF grib files are "eta_dot"
c FLEXPART uses a "height above ground" vertical coordinate
c    which we will call "hag".
c    The vertical velocity is uses (in ww array) is "hag_dot".
c Converting from eta_dot to hag_dot involves
c    >> multiplying by pinmconv = [d(hag)/d(eta)]
c    >> adding a term that accounts for the fact that
c       "eta" varies on constant "hag" surfaces.
c       This term is [u*d(hag)/dx + v*d(hag)/dy], with the
c       partial derivatives taken with "eta" being constant
c
c The WRF model uses a similar (to ECMWF) vertical coordinate.
c    HOWEVER, the vertical velocities in the WRF output files
c    are the "true/absolute w" in m/s.  (Is this true?)
c Converting from "absolute w" to hag_dot involves
c    adding a term that accounts for the fact that
c    "absolute z" varies on constant "hag" surfaces.
c    This term is [- u*d(oro)/dx - v*d(oro)/dy]
c
c The FLEXPART code did not apply the terrain corrections
c    at jy=0 & ny-1; ix=0 & nx-1; iz=1 & nz.
c FLEXPART_WRF applies the correction at all grid points
*****************************************************************

      if (method_w_terrain_correction .eq. 1) then
c
c apply w correction as it is done with ecmwf met
c
      write(*,*) 'verttransform doing method_w_terrain_correction =', 
     &           method_w_terrain_correction

      do 60 jy=1,ny-2
        do 60 ix=1,nx-2
c     do 60 jy=0,ny-1
c       do 60 ix=0,nx-1

          kmin=2
          do 41 iz=2,nz-1

c following converts from (m/s) to ("grid index coord. units"/s)
c for FLEXPART_WRF, dxconst=1/dx, dyconst=1/dy, and no cos(lat) is needed
c           ui=uu(ix,jy,iz,n)*dxconst/cos((float(jy)*dy+ylat0)*pi180)
            ui=uu(ix,jy,iz,n)*dxconst
            vi=vv(ix,jy,iz,n)*dyconst

            do 46 kz=kmin,nz
              if ((height(iz).gt.uvwzlev(ix,jy,kz-1)).and.
     +            (height(iz).le.uvwzlev(ix,jy,kz))) then
                dz1=height(iz)-uvwzlev(ix,jy,kz-1)
                dz2=uvwzlev(ix,jy,kz)-height(iz)
                dz=dz1+dz2
                kl=kz-1
                klp=kz
                kmin=kz
                goto 47
              endif
46            continue

c47         ix1=ix-1
c           jy1=jy-1
c           ixp=ix+1
c           jyp=jy+1
47          ix1 = max( ix-1, 0 )
            jy1 = max( jy-1, 0 )
            ixp = min( ix+1, nx-1 )
            jyp = min( jy+1, ny-1 )

c           dzdx1=(uvwzlev(ixp,jy,kl)-uvwzlev(ix1,jy,kl))/2.
c           dzdx2=(uvwzlev(ixp,jy,klp)-uvwzlev(ix1,jy,klp))/2.
            dzdx1=(uvwzlev(ixp,jy,kl )-uvwzlev(ix1,jy,kl ))/(ixp-ix1)
            dzdx2=(uvwzlev(ixp,jy,klp)-uvwzlev(ix1,jy,klp))/(ixp-ix1)
            dzdx=(dzdx1*dz2+dzdx2*dz1)/dz

c           dzdy1=(uvwzlev(ix,jyp,kl)-uvwzlev(ix,jy1,kl))/2.
c           dzdy2=(uvwzlev(ix,jyp,klp)-uvwzlev(ix,jy1,klp))/2.
            dzdy1=(uvwzlev(ix,jyp,kl )-uvwzlev(ix,jy1,kl ))/(jyp-jy1)
            dzdy2=(uvwzlev(ix,jyp,klp)-uvwzlev(ix,jy1,klp))/(jyp-jy1)
            dzdy=(dzdy1*dz2+dzdy2*dz1)/dz

            ww(ix,jy,iz,n)=ww(ix,jy,iz,n)+(dzdx*ui+dzdy*vi)
 
41          continue
 
60        continue


      else if ((method_w_terrain_correction .eq. 10) .or.
     &         (method_w_terrain_correction .eq. 11)) then
c
c apply w correction assuming that the WRF w is "absolute w"
c
c also, set ww=0 at iz=1 when method_w_terrain_correction=11
c
      write(*,*) 'verttransform doing method_w_terrain_correction =', 
     &           method_w_terrain_correction

      do 65 jy=0,ny-1
      do 65 ix=0,nx-1
          ix1 = max( ix-1, 0 )
          jy1 = max( jy-1, 0 )
          ixp = min( ix+1, nx-1 )
          jyp = min( jy+1, ny-1 )
          dzdx=(oro(ixp,jy) - oro(ix1,jy))/(dx*(ixp-ix1))
          dzdy=(oro(ix,jyp) - oro(ix,jy1))/(dy*(jyp-jy1))

          do iz=1,nz
            wtc_stat(1,iz) = wtc_stat(1,iz) + ww(ix,jy,iz,n)
            wtc_stat(2,iz) = wtc_stat(2,iz) + abs(ww(ix,jy,iz,n))

            ww(ix,jy,iz,n) = ww(ix,jy,iz,n)
     &           - (uu(ix,jy,iz,n)*dzdx + vv(ix,jy,iz,n)*dzdy)  
            if ((method_w_terrain_correction .eq. 11) .and.
     &          (iz .eq. 1)) ww(ix,jy,iz,n) = 0.0

            wtc_stat(3,iz) = wtc_stat(3,iz) + ww(ix,jy,iz,n)
            wtc_stat(4,iz) = wtc_stat(4,iz) + abs(ww(ix,jy,iz,n))
          end do
65    continue

      else if (method_w_terrain_correction .eq. 20) then
          continue
      else
      write(*,*) 'verttransform SKIPPING method_w_terrain_correction'
      end if

c turn off this output for now
      if (nz .le. -123) then
      write(*,'(/a/a)') 'verttransform - ww statistics',
     &    '  k    avg & abs before        avg & abs after '
      do iz = 1, nz
        if ((iz .le. 5) .or. (iz .ge. nwz-4)
     *                  .or. (mod(iz,2) .eq. 1)) then
          do ix = 1, 4
            wtc_stat(ix,iz) = wtc_stat(ix,iz)/(nx*ny)
          end do
          write(*,'(i3,1p,2e10.2,4x,2e10.2)') iz,
     &        (wtc_stat(ix,iz), ix=1,4)
        end if
      end do
      end if



C If north pole is in the domain, calculate wind velocities in polar
C stereographic coordinates
********************************************************************
 
      if (nglobal) then
        write(*,*)
        write(*,*) '*** stopping in verttransform ***'
        write(*,*) '    the nglobal code section should not be active'
        write(*,*)
        stop
c        do 74 jy=int(switchnorthg)-2,nymin1
c          ylat=ylat0+float(jy)*dy
c          do 74 ix=0,nxmin1
c            xlon=xlon0+float(ix)*dx
c            do 74 iz=1,nz
c74            call cc2gll(northpolemap,ylat,xlon,uu(ix,jy,iz,n),
c     +        vv(ix,jy,iz,n),uupol(ix,jy,iz,n),
c     +        vvpol(ix,jy,iz,n))
c
c
c        do 76 iz=1,nz
c
c* CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
c          xlon=xlon0+float(nx/2-1)*dx
c          xlonr=xlon*pi/180.
c          ffpol=sqrt(uu(nx/2-1,nymin1,iz,n)**2+
c     &               vv(nx/2-1,nymin1,iz,n)**2)
c          if(vv(nx/2-1,nymin1,iz,n).lt.0.) then
c            ddpol=atan(uu(nx/2-1,nymin1,iz,n)/
c     &                 vv(nx/2-1,nymin1,iz,n))-xlonr
c          else
c            ddpol=pi+atan(uu(nx/2-1,nymin1,iz,n)/
c     &                    vv(nx/2-1,nymin1,iz,n))-xlonr
c          endif
c          if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
c          if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi
c 
c* CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
c          xlon=180.0
c          xlonr=xlon*pi/180.
c          ylat=90.0
c          uuaux=-ffpol*sin(xlonr+ddpol)
c          vvaux=-ffpol*cos(xlonr+ddpol)
c          call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux,
c     +      vvpolaux)
c      
c          jy=nymin1
c          do 76 ix=0,nxmin1
c            uupol(ix,jy,iz,n)=uupolaux
c            vvpol(ix,jy,iz,n)=vvpolaux
c76      continue
c 
c 
c* Fix: Set W at pole to the zonally averaged W of the next equator-
c* ward parallel of latitude
c 
c      do 85 iz=1,nz
c          wdummy=0.
c          jy=ny-2
c          do 80 ix=0,nxmin1
c80          wdummy=wdummy+ww(ix,jy,iz,n)
c          wdummy=wdummy/float(nx)
c          jy=nymin1
c          do 85 ix=0,nxmin1
c85          ww(ix,jy,iz,n)=wdummy
 
      endif 

 
C If south pole is in the domain, calculate wind velocities in polar
C stereographic coordinates
********************************************************************
 
      if (sglobal) then
        write(*,*)
        write(*,*) '*** stopping in verttransform ***'
        write(*,*) '    the sglobal code section should not be active'
        write(*,*)
        stop
c        do 77 jy=0,int(switchsouthg)+3
c          ylat=ylat0+float(jy)*dy
c          do 77 ix=0,nxmin1
c            xlon=xlon0+float(ix)*dx
c            do 77 iz=1,nz
c77            call cc2gll(southpolemap,ylat,xlon,uu(ix,jy,iz,n),
c     +        vv(ix,jy,iz,n),uupol(ix,jy,iz,n),
c     +        vvpol(ix,jy,iz,n))
c      
c        do 79 iz=1,nz
c 
c* CALCULATE FFPOL, DDPOL FOR CENTRAL GRID POINT
c          xlon=xlon0+float(nx/2-1)*dx
c          xlonr=xlon*pi/180.
c          ffpol=sqrt(uu(nx/2-1,0,iz,n)**2+
c     &               vv(nx/2-1,0,iz,n)**2)
c          if(vv(nx/2-1,0,iz,n).lt.0.) then
c            ddpol=atan(uu(nx/2-1,0,iz,n)/
c     &                 vv(nx/2-1,0,iz,n))+xlonr
c          else
c            ddpol=pi+atan(uu(nx/2-1,0,iz,n)/
c     &                    vv(nx/2-1,0,iz,n))+xlonr
c          endif
c          if(ddpol.lt.0.) ddpol=2.0*pi+ddpol
c          if(ddpol.gt.2.0*pi) ddpol=ddpol-2.0*pi
c 
c* CALCULATE U,V FOR 180 DEG, TRANSFORM TO POLAR STEREOGRAPHIC GRID
c          xlon=180.0
c          xlonr=xlon*pi/180.
c          ylat=-90.0
c          uuaux=+ffpol*sin(xlonr-ddpol)
c          vvaux=-ffpol*cos(xlonr-ddpol)
c          call cc2gll(northpolemap,ylat,xlon,uuaux,vvaux,uupolaux,
c     +      vvpolaux)
c      
c          jy=0
c          do 79 ix=0,nxmin1
c            uupol(ix,jy,iz,n)=uupolaux
c79          vvpol(ix,jy,iz,n)=vvpolaux
c 
c 
c* Fix: Set W at pole to the zonally averaged W of the next equator-
c* ward parallel of latitude
c 
c        do 95 iz=1,nz
c          wdummy=0.
c          jy=1
c          do 90 ix=0,nxmin1
c90          wdummy=wdummy+ww(ix,jy,iz,n)
c          wdummy=wdummy/float(nx)
c          jy=0
c          do 95 ix=0,nxmin1
c95          ww(ix,jy,iz,n)=wdummy
      endif



      end
