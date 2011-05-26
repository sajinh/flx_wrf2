      subroutine verttransform_nests(n,uuhn,vvhn,wwhn,pvhn)
C                                    i   i    i    i   i
********************************************************************************
*                                                                              *
* Note:  This is the FLEXPART_WRF version of subroutine verttransform_nests.   *
*     The computational grid is the WRF x-y grid rather than lat-lon.          *
*                                                                              *
*     This subroutine transforms temperature, dew point temperature and        *
*     wind components from eta to meter coordinates.                           *
*     The vertical wind component is transformed from Pa/s to m/s using        *
*     the conversion factor pinmconv.                                          *
*     In addition, this routine calculates vertical density gradients          *
*     needed for the parameterization of the turbulent velocities.             *
*     It is similar to verttransform, but makes the transformations for        *
*     the nested grids.                                                        *
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
*     Changes, Bernd C. Krueger, Feb. 2001:       (marked "C-cv")              *
*        Variables tthn and qvhn (on eta coordinates) from common block        *
*                                                                              *
*     16 Nov 2005, R. Easter - changes for FLEXPART_WRF                        *
*     17 Nov 2005 - R. Easter - terrain correction applied to ww.  There are   *
*            now 3 options, controlled by "method_w_terrain_correction"        *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* nxn,nyn,nuvz,nwz                field dimensions in x,y and z direction      *
* uun                             wind components in x-direction [m/s]         *
* vvn                             wind components in y-direction [m/s]         *
* wwn                             wind components in z-direction [deltaeta/s]  *
* ttn                             temperature [K]                              *
* pvn                             potential vorticity (pvu)                    *
* psn                             surface pressure [Pa]                        *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      integer ix,jy,kz,iz,n,l,kmin,kl,klp,ix1,jy1,ixp,jyp
      integer method_z_compute
      real uvzlev(nuvzmax),wzlev(nwzmax),rhoh(nuvzmax),pinmconv(nzmax)
      real uvwzlev(0:nxmaxn-1,0:nymaxn-1,nzmax)
      real ew,pint,tv,tvold,pold,const,dz1,dz2,dz,ui,vi
      real dzdx,dzdy
      real dzdx1,dzdx2,dzdy1,dzdy2
      real uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
      real vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
      real pvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
      real wwhn(0:nxmaxn-1,0:nymaxn-1,nwzmax,maxnests)
      real wwhn_svaa(nwzmax)
      parameter(const=r_air/ga)



c set method_z_compute
      method_z_compute = 10


C Loop over all nests
*********************

      do 100 l=1,numbnests

C Loop over the whole grid
**************************

      do 10 jy=0,nyn(l)-1
        do 10 ix=0,nxn(l)-1

          tvold=tt2n(ix,jy,1,n,l)*(1.+0.378*ew(td2n(ix,jy,1,n,l))/
     &                                   psn(ix,jy,1,n,l))
          pold=psn(ix,jy,1,n,l)
          uvzlev(1)=0.
          wzlev(1)=0.
          rhoh(1)=pold/(r_air*tvold)


C Compute heights of eta levels
*******************************

          do 20 kz=2,nuvz
c FLEXPART_WRF - pphn hold pressure
c           pint=akz(kz)+bkz(kz)*psn(ix,jy,1,n,l)
            pint=pphn(ix,jy,kz,n,l)
            tv=tthn(ix,jy,kz,n,l)*(1.+0.608*qvhn(ix,jy,kz,n,l))
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
                uvzlev(kz) = 0.5*(zzhn(ix,jy,3,n,l) + 
     &                            zzhn(ix,jy,1,n,l))
     &                     - zzhn(ix,jy,1,n,l)
              else
                uvzlev(kz) = 0.5*(zzhn(ix,jy,kz+1,n,l) + 
     &                            zzhn(ix,jy,kz  ,n,l))
     &                     - zzhn(ix,jy,1,n,l)
              end if
            end do
            do kz = 2, nwz
              wzlev(kz) = zzhn(ix,jy,kz+add_sfc_level,n,l) 
     &                  - zzhn(ix,jy,1,n,l)
            end do
          end if

C NOTE: In FLEXPART versions up to 4.0, the number of model levels was doubled
C upon the transformation to z levels. In order to save computer memory, this is
C not done anymore in the standard version. However, this option can still be
C switched on by replacing the following lines with those below, that are
C currently commented out.  
C Note that one change is also necessary in gridcheck.f,
C and three changes in verttransform.f
C
C *** NOTE -- the doubled vertical resolution has not been tested in FLEXPART_WRF
********************************************************************************
          uvwzlev(ix,jy,1)=0.0
          do 22 kz=2,nuvz
22          uvwzlev(ix,jy,kz)=uvzlev(kz)

C Switch on following lines to use doubled vertical resolution
C Switch off the three lines above.
C
C *** NOTE -- the doubled vertical resolution has not been tested in FLEXPART_WRF
**************************************************************
c22          uvwzlev(ix,jy,(kz-1)*2)=uvzlev(kz)
c          do 23 kz=2,nwz
c23          uvwzlev(ix,jy,(kz-1)*2+1)=wzlev(kz)
C End doubled vertical resolution

C pinmconv=(h2-h1)/(p2-p1)
C
C in flexpart_ecmwf, pinmconv is used to convert etadot to w
C in FLEXPART_WRF, vertical velocity is already m/s, so pinmconv=1.0
C
c         pinmconv(1)=(uvwzlev(ix,jy,2)-uvwzlev(ix,jy,1))/
c    +    ((aknew(2)+bknew(2)*psn(ix,jy,1,n,l))-
c    +    (aknew(1)+bknew(1)*psn(ix,jy,1,n,l)))
          pinmconv(1)=1.0
          do 24 kz=2,nz-1
c           pinmconv(kz)=(uvwzlev(ix,jy,kz+1)-uvwzlev(ix,jy,kz-1))/
c    +      ((aknew(kz+1)+bknew(kz+1)*psn(ix,jy,1,n,l))-
c    +      (aknew(kz-1)+bknew(kz-1)*psn(ix,jy,1,n,l)))
            pinmconv(kz)=1.0
24        continue
c         pinmconv(nz)=(uvwzlev(ix,jy,nz)-uvwzlev(ix,jy,nz-1))/
c    +    ((aknew(nz)+bknew(nz)*psn(ix,jy,1,n,l))-
c    +    (aknew(nz-1)+bknew(nz-1)*psn(ix,jy,1,n,l)))
          pinmconv(nz)=1.0


C Levels, where u,v,t and q are given
*************************************

          uun(ix,jy,1,n,l)=uuhn(ix,jy,1,l)
          vvn(ix,jy,1,n,l)=vvhn(ix,jy,1,l)
          ttn(ix,jy,1,n,l)=tthn(ix,jy,1,n,l)
          qvn(ix,jy,1,n,l)=qvhn(ix,jy,1,n,l)
          pvn(ix,jy,1,n,l)=pvhn(ix,jy,1,l)
          rhon(ix,jy,1,n,l)=rhoh(1)
          uun(ix,jy,nz,n,l)=uuhn(ix,jy,nuvz,l)
          vvn(ix,jy,nz,n,l)=vvhn(ix,jy,nuvz,l)
          ttn(ix,jy,nz,n,l)=tthn(ix,jy,nuvz,n,l)
          qvn(ix,jy,nz,n,l)=qvhn(ix,jy,nuvz,n,l)
          pvn(ix,jy,nz,n,l)=pvhn(ix,jy,nuvz,l)
          rhon(ix,jy,nz,n,l)=rhoh(nuvz)
          kmin=2
          do 30 iz=2,nz-1
            do 35 kz=kmin,nuvz
              if(height(iz).gt.uvzlev(nuvz)) then
                uun(ix,jy,iz,n,l)=uun(ix,jy,nz,n,l)
                vvn(ix,jy,iz,n,l)=vvn(ix,jy,nz,n,l)
                ttn(ix,jy,iz,n,l)=ttn(ix,jy,nz,n,l)
                qvn(ix,jy,iz,n,l)=qvn(ix,jy,nz,n,l)
                pvn(ix,jy,iz,n,l)=pvn(ix,jy,nz,n,l)
                rhon(ix,jy,iz,n,l)=rhon(ix,jy,nz,n,l)
                goto 30
              endif
              if ((height(iz).gt.uvzlev(kz-1)).and.
     +            (height(iz).le.uvzlev(kz))) then
               dz1=height(iz)-uvzlev(kz-1)
               dz2=uvzlev(kz)-height(iz)
               dz=dz1+dz2
               uun(ix,jy,iz,n,l)=(uuhn(ix,jy,kz-1,l)*dz2+
     +         uuhn(ix,jy,kz,l)*dz1)/dz
               vvn(ix,jy,iz,n,l)=(vvhn(ix,jy,kz-1,l)*dz2+
     +         vvhn(ix,jy,kz,l)*dz1)/dz
               ttn(ix,jy,iz,n,l)=(tthn(ix,jy,kz-1,n,l)*dz2+
     +         tthn(ix,jy,kz,n,l)*dz1)/dz
               qvn(ix,jy,iz,n,l)=(qvhn(ix,jy,kz-1,n,l)*dz2+
     +         qvhn(ix,jy,kz,n,l)*dz1)/dz
               pvn(ix,jy,iz,n,l)=(pvhn(ix,jy,kz-1,l)*dz2+
     +         pvhn(ix,jy,kz,l)*dz1)/dz
               rhon(ix,jy,iz,n,l)=(rhoh(kz-1)*dz2+rhoh(kz)*dz1)/dz
               kmin=kz
               goto 30
              endif
35            continue
30          continue


C Levels, where w is given
**************************

          if (method_w_terrain_correction .eq. 20) then
c apply w correction assuming that the WRF w is "absolute w";
c apply it here to wwh; set wwh=0 at iz=1
             ix1 = max( ix-1, 0 )
             jy1 = max( jy-1, 0 )
             ixp = min( ix+1, nxn(l)-1 )
             jyp = min( jy+1, nyn(l)-1 )
             dzdx=(oron(ixp,jy,l) - oron(ix1,jy,l))/(dxn(l)*(ixp-ix1))
             dzdy=(oron(ix,jyp,l) - oron(ix,jy1,l))/(dyn(l)*(jyp-jy1))

             do kz = 1, nwz
                wwhn_svaa(kz) = wwhn(ix,jy,kz,l)
                wwhn(ix,jy,kz,l) = wwhn(ix,jy,kz,l)
     &               - (uuhn(ix,jy,kz,l)*dzdx + vvhn(ix,jy,kz,l)*dzdy)  
                if (kz .eq. 1) wwhn(ix,jy,kz,l) = 0.0
             end do
          end if

          wwn(ix,jy,1,n,l)=wwhn(ix,jy,1,l)*pinmconv(1)
          wwn(ix,jy,nz,n,l)=wwhn(ix,jy,nwz,l)*pinmconv(nz)
          kmin=2
          do 40 iz=2,nz
            do 45 kz=kmin,nwz
              if ((height(iz).gt.wzlev(kz-1)).and.
     +            (height(iz).le.wzlev(kz))) then
               dz1=height(iz)-wzlev(kz-1)
               dz2=wzlev(kz)-height(iz)
               dz=dz1+dz2
               wwn(ix,jy,iz,n,l)=(wwhn(ix,jy,kz-1,l)*dz2+
     +         wwhn(ix,jy,kz,l)*dz1)/dz*pinmconv(iz)
               kmin=kz
               goto 40
              endif
45            continue
40          continue

          if (method_w_terrain_correction .eq. 20) then
             do kz = 1, nwz
                wwhn(ix,jy,kz,l) = wwhn_svaa(kz)
             end do
          end if

C Compute density gradients at intermediate levels
**************************************************

          drhodzn(ix,jy,1,n,l)=(rhon(ix,jy,2,n,l)-rhon(ix,jy,1,n,l))/
     +      (height(2)-height(1))
          do 50 kz=2,nz-1
50          drhodzn(ix,jy,kz,n,l)=(rhon(ix,jy,kz+1,n,l)-
     +      rhon(ix,jy,kz-1,n,l))/(height(kz+1)-height(kz-1))
          drhodzn(ix,jy,nz,n,l)=drhodzn(ix,jy,nz-1,n,l)

10        continue


*****************************************************************
C Compute slope of eta levels in windward direction and resulting
C vertical wind correction
c
c See notes in verttransform.f about the w correction done here.
*****************************************************************

      if (method_w_terrain_correction .eq. 1) then
c
c apply w correction as it is done with ecmwf met
c
      do 60 jy=1,nyn(l)-2
        do 60 ix=1,nxn(l)-2
c     do 60 jy=0,nyn(l)-1
c       do 60 ix=0,nxn(l)-1

          kmin=2
          do 41 iz=2,nz-1

c for FLEXPART_WRF, dx & dy are in meters,
c dxconst=1/dx, dyconst=1/dy, and no cos(lat) is needed
c           ui=uun(ix,jy,iz,n,l)*dxconst*xresoln(l)/
c    +      cos((float(jy)*dyn(l)+ylat0n(l))*pi180)
            ui=uun(ix,jy,iz,n,l)*dxconst*xresoln(l)
            vi=vvn(ix,jy,iz,n,l)*dyconst*yresoln(l)

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
            ixp = min( ix+1, nxn(l)-1 )
            jyp = min( jy+1, nyn(l)-1 )

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

            wwn(ix,jy,iz,n,l)=wwn(ix,jy,iz,n,l)+(dzdx*ui+dzdy*vi)
 
41          continue
 
60        continue

      else if ((method_w_terrain_correction .eq. 10) .or.
     &         (method_w_terrain_correction .eq. 11)) then
c
c apply w correction assuming that the WRF w is "absolute w"
c
c also, set ww=0 at iz=1 when method_w_terrain_correction=11
c
      do 65 jy=0,nyn(l)-1
      do 65 ix=0,nxn(l)-1
          ix1 = max( ix-1, 0 )
          jy1 = max( jy-1, 0 )
          ixp = min( ix+1, nxn(l)-1 )
          jyp = min( jy+1, nyn(l)-1 )
          dzdx=(oron(ixp,jy,l) - oron(ix1,jy,l))/(dxn(l)*(ixp-ix1))
          dzdy=(oron(ix,jyp,l) - oron(ix,jy1,l))/(dyn(l)*(jyp-jy1))

          do iz=1,nz
            wwn(ix,jy,iz,n,l) = wwn(ix,jy,iz,n,l)
     &           - (uun(ix,jy,iz,n,l)*dzdx + vvn(ix,jy,iz,n,l)*dzdy)  
            if ((method_w_terrain_correction .eq. 11) .and.
     &          (iz .eq. 1)) wwn(ix,jy,iz,n,l) = 0.0

          end do
65    continue

c
c "else" -- apply no w correction
c
      end if



100   continue

      return
      end
