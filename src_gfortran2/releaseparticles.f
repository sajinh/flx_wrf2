      subroutine releaseparticles(itime)
C                                   o
********************************************************************************
*                                                                              *
*     This subroutine releases particles from the release locations.           *
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subroutine gridcheck.         *
*            The computational grid is the WRF x-y grid rather than lat-lon.   *
*                                                                              *
*     It searches for a "vacant" storage space and assigns all particle        *
*     information to that space. A space is vacant either when no particle     *
*     is yet assigned to it, or when it's particle is expired and, thus,       *
*     the storage space is made available to a new particle.                   *
*                                                                              *
*     Author: A. Stohl                                                         *
*     29 June 2002                                                             *
*                                                                              *
*     14 Nov 2005, R. Easter - use xyindex_to_ll_wrf to get lat,lon            *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* itime [s]            current time                                            *
* ireleasestart, ireleaseend          start and end times of all releases      *
* npart(maxpoint)      number of particles to be released in total             *
* numrel               number of particles to be released during this time step*
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      real xaux,yaux,zaux,ran1,fraction,xmasssave(maxpoint)
      real topo,rhoaux(2),r,t,rhoout,ddx,ddy,rddx,rddy,p1,p2,p3,p4,eps
      real dz1,dz2,dz,xtn,ytn,xlonav,timecorrect(maxspec),press,pressold
      real presspart,average_timecorrect
      real dumx, dumy, ylatav
      parameter(eps=1.e-6)
      integer itime,numrel,idummy,i,j,k,n,ix,jy,ixp,jyp,ipart,minpart,ii
      integer indz,indzp,kz,ngrid
      integer nweeks,ndayofweek,nhour,jjjjmmdd,ihmmss,mm
      double precision juldate,julmonday,jul,jullocal,juldiff

      save idummy,xmasssave
      data idummy/-7/,xmasssave/maxpoint*0./



C Determine the actual date and time in Greenwich (i.e., UTC + correction for daylight savings time)
****************************************************************************************************

      julmonday=juldate(19000101,0)          ! this is a Monday
      jul=bdate+dble(float(itime)/86400.)    ! this is the current day
      call caldate(jul,jjjjmmdd,ihmmss)
      mm=(jjjjmmdd-10000*(jjjjmmdd/10000))/100
      if ((mm.ge.4).and.(mm.le.9)) jul=jul+dble(1./24.)   ! daylight savings time in summer
          

C For every release point, check whether we are in the release time interval
****************************************************************************

      minpart=1
      do 10 i=1,numpoint
        if ((itime.ge.ireleasestart(i)).and.             ! are we within release interval?
     +  (itime.le.ireleaseend(i))) then

C Determine the local day and time
**********************************

c FLEXPART_WRF - use this routine to get lat,lon
c         xlonav=xlon0+(xpoint2(i)+xpoint1(i))/2.*dx  ! longitude needed to determine local time
          dumx = (xpoint2(i)+xpoint1(i))*0.5
          dumy = (ypoint2(i)+ypoint1(i))*0.5
          call xyindex_to_ll_wrf( 0, dumx, dumy, xlonav, ylatav )

          if (xlonav.lt.-180.) xlonav=xlonav+360.
          if (xlonav.gt.180.) xlonav=xlonav-360.
          jullocal=jul+dble(xlonav/360.)   ! correct approximately for time zone to obtain local time

          juldiff=jullocal-julmonday
          nweeks=int(juldiff/7.)
          juldiff=juldiff-dble(float(nweeks)*7.)
          ndayofweek=int(juldiff)+1              ! this is the current day of week, starting with Monday
          nhour=nint((juldiff-dble(float(ndayofweek-1)))*24.)    ! this is the current hour
          if (nhour.eq.0) then
            nhour=24
            ndayofweek=ndayofweek-1
            if (ndayofweek.eq.0) ndayofweek=7
          endif

C Calculate a species- and time-dependent correction factor, distinguishing between
C area (those with release starting at surface) and point (release starting above surface) sources
C Also, calculate an average time correction factor (species independent)
**************************************************************************************************
          average_timecorrect=0.
          do 25 k=1,nspec
            if (zpoint1(i).gt.0.5) then      ! point source
              timecorrect(k)=point_hour(k,nhour)*point_dow(k,ndayofweek)
            else                             ! area source
              timecorrect(k)=area_hour(k,nhour)*area_dow(k,ndayofweek)
            endif
            average_timecorrect=average_timecorrect+timecorrect(k)
25          continue
          average_timecorrect=average_timecorrect/float(nspec)

C Determine number of particles to be released this time; at start and at end of release,
C only half the particles are released
*****************************************************************************************

          if (ireleasestart(i).ne.ireleaseend(i)) then
            fraction=abs(float(npart(i))*float(lsynctime)/
     +      float(ireleaseend(i)-ireleasestart(i)))
            if ((itime.eq.ireleasestart(i)).or.
     +      (itime.eq.ireleaseend(i))) fraction=fraction/2.

C Take the species-average time correction factor in order to scale the
C number of particles released this time
***********************************************************************
            fraction=fraction*average_timecorrect

            fraction=fraction+xmasssave(i)  ! number to be released at this time
            numrel=int(fraction)
            xmasssave(i)=fraction-float(numrel)
          else
            numrel=npart(i)
          endif

          xaux=xpoint2(i)-xpoint1(i)
          yaux=ypoint2(i)-ypoint1(i)
          zaux=zpoint2(i)-zpoint1(i)
          do 30 j=1,numrel                       ! loop over particles to be released this time
            do 33 ipart=minpart,maxpart          ! search for free storage space

C If a free storage space is found, attribute everything to this array element
******************************************************************************

              if (itra1(ipart).ne.itime) then   

C Particle coordinates are determined by using a random position within the release volume
******************************************************************************************

C Determine horizontal particle position
****************************************

                xtra1(ipart)=xpoint1(i)+ran1(idummy)*xaux
                if (xglobal) then
                  if (xtra1(ipart).gt.float(nxmin1)) xtra1(ipart)=
     +            xtra1(ipart)-float(nxmin1)
                  if (xtra1(ipart).lt.0.) xtra1(ipart)=
     +            xtra1(ipart)+float(nxmin1)
                endif
                ytra1(ipart)=ypoint1(i)+ran1(idummy)*yaux

C Assign mass to particle: Total mass divided by total number of particles.
C Time variation has partly been taken into account already by a species-average
C correction factor, by which the number of particles released this time has been
C scaled. Adjust the mass per particle by the species-dependent time correction factor
C divided by the species-average one  
**************************************************************************************

                do 20 k=1,nspec
20                xmass1(ipart,k)=xmass(i,k)/float(npart(i))
     +            *timecorrect(k)/average_timecorrect

C Assign certain properties to particle
***************************************

                nclass(ipart)=min(int(ran1(idummy)*float(nclassunc))+1,
     +          nclassunc)
                numparticlecount=numparticlecount+1
                if (mquasilag.eq.0) then
                  npoint(ipart)=i
                else
                  npoint(ipart)=numparticlecount
                endif
                idt(ipart)=mintime               ! first time step
                itra1(ipart)=itime
                itramem(ipart)=itra1(ipart)
                itrasplit(ipart)=itra1(ipart)+ldirect*itsplit


C Determine vertical particle position
**************************************

                ztra1(ipart)=zpoint1(i)+ran1(idummy)*zaux

C Interpolation of topography and density
*****************************************

C Determine the nest we are in
******************************

                ngrid=0
                do 42 k=numbnests,1,-1
                  if ((xtra1(ipart).gt.xln(k)).and.
     +            (xtra1(ipart).lt.xrn(k)).and.
     +            (ytra1(ipart).gt.yln(k)).and.
     +            (ytra1(ipart).lt.yrn(k))) then
                    ngrid=k
                    goto 43
                  endif
42                continue
43              continue

C Determine (nested) grid coordinates and auxiliary parameters used for interpolation
*************************************************************************************

                if (ngrid.gt.0) then
                  xtn=(xtra1(ipart)-xln(ngrid))*xresoln(ngrid)
                  ytn=(ytra1(ipart)-yln(ngrid))*yresoln(ngrid)
                  ix=int(xtn)
                  jy=int(ytn)
                  ddy=ytn-float(jy)
                  ddx=xtn-float(ix)
                else
                  ix=int(xtra1(ipart))
                  jy=int(ytra1(ipart))
                  ddy=ytra1(ipart)-float(jy)
                  ddx=xtra1(ipart)-float(ix)
                endif
                ixp=ix+1
                jyp=jy+1
                rddx=1.-ddx
                rddy=1.-ddy
                p1=rddx*rddy
                p2=ddx*rddy
                p3=rddx*ddy
                p4=ddx*ddy

                if (ngrid.gt.0) then
                  topo=p1*oron(ix ,jy ,ngrid)
     +               + p2*oron(ixp,jy ,ngrid)
     +               + p3*oron(ix ,jyp,ngrid)
     +               + p4*oron(ixp,jyp,ngrid)
                else
                  topo=p1*oro(ix ,jy)
     +               + p2*oro(ixp,jy)
     +               + p3*oro(ix ,jyp)
     +               + p4*oro(ixp,jyp)
                endif

C If starting height is in pressure coordinates, retrieve pressure profile and convert zpart1 to meters
*******************************************************************************************************
                if (kindz(i).eq.3) then
                  presspart=ztra1(ipart)
                  do 70 kz=1,nz
                    if (ngrid.gt.0) then
                      r=p1*rhon(ix ,jy ,kz,2,ngrid)
     +                 +p2*rhon(ixp,jy ,kz,2,ngrid)
     +                 +p3*rhon(ix ,jyp,kz,2,ngrid)
     +                 +p4*rhon(ixp,jyp,kz,2,ngrid)
                      t=p1*ttn(ix ,jy ,kz,2,ngrid)
     +                 +p2*ttn(ixp,jy ,kz,2,ngrid)
     +                 +p3*ttn(ix ,jyp,kz,2,ngrid)
     +                 +p4*ttn(ixp,jyp,kz,2,ngrid)
                    else
                      r=p1*rho(ix ,jy ,kz,2)
     +                 +p2*rho(ixp,jy ,kz,2)
     +                 +p3*rho(ix ,jyp,kz,2)
     +                 +p4*rho(ixp,jyp,kz,2)
                      t=p1*tt(ix ,jy ,kz,2)
     +                 +p2*tt(ixp,jy ,kz,2)
     +                 +p3*tt(ix ,jyp,kz,2)
     +                 +p4*tt(ixp,jyp,kz,2)
                    endif
                    press=r*r_air*t/100.

                    if (press.lt.presspart) then
                      if (kz.eq.1) then
                        ztra1(ipart)=height(1)/2.
                      else
                        dz1=pressold-presspart
                        dz2=presspart-press
                        ztra1(ipart)=(height(kz-1)*dz2+height(kz)*dz1)
     +                  /(dz1+dz2)
                      endif
                      goto 71
                    endif
                    pressold=press
70                  continue
71                continue
                endif

C If release positions are given in meters above sea level, subtract the
C topography from the starting height
************************************************************************

                if (kindz(i).eq.2) ztra1(ipart)=ztra1(ipart)-topo
                if (ztra1(ipart).lt.eps) ztra1(ipart)=eps   ! Minimum starting height is eps
                if (ztra1(ipart).gt.height(nz)-0.5) ztra1(ipart)=
     +          height(nz)-0.5 ! Maximum starting height is uppermost level - 0.5 meters



C For special simulations, multiply particle concentration air density;
C Simply take the 2nd field in memory to do this (accurate enough)
************************************************************************
cAF IND_SOURCE switches between different units for concentrations at the source
cAf    NOTE that in backward simulations the release of particles takes place at the 
cAf         receptor and the sampling at the source.
cAf          1="mass" 
cAf          2="mass mixing ratio" 
cAf IND_RECEPTOR switches between different units for concentrations at the receptor
cAf          1="mass" 
cAf          2="mass mixing ratio" 

cAf switches for the releasefile:
cAf IND_REL =  1 : xmass * rho
cAf IND_REL =  0 : xmass * 1

cAf ind_rel is defined in readcommand.f

                if (ind_rel .eq. 1) then

C Interpolate the air density
*****************************

                  do 5 ii=2,nz
                    if (height(ii).gt.ztra1(ipart)) then
                      indz=ii-1
                      indzp=ii
                      goto 6
                    endif
5                   continue
6                 continue

                  dz1=ztra1(ipart)-height(indz)
                  dz2=height(indzp)-ztra1(ipart)
                  dz=1./(dz1+dz2)

                  if (ngrid.gt.0) then
                    do 44 n=1,2
44                    rhoaux(n)=p1*rhon(ix ,jy ,indz+n-1,2,ngrid)
     +                         +p2*rhon(ixp,jy ,indz+n-1,2,ngrid)
     +                         +p3*rhon(ix ,jyp,indz+n-1,2,ngrid)
     +                         +p4*rhon(ixp,jyp,indz+n-1,2,ngrid)
                  else
                    do 45 n=1,2
45                    rhoaux(n)=p1*rho(ix ,jy ,indz+n-1,2)
     +                         +p2*rho(ixp,jy ,indz+n-1,2)
     +                         +p3*rho(ix ,jyp,indz+n-1,2)
     +                         +p4*rho(ixp,jyp,indz+n-1,2)
                  endif
                  rhoout=(dz2*rhoaux(1)+dz1*rhoaux(2))*dz


C Multiply "mass" (i.e., mass mixing ratio in forward runs) with density
*********************************************************************

                  do 46 k=1,nspec
46                  xmass1(ipart,k)=xmass1(ipart,k)*rhoout
                endif



                numpart=max(numpart,ipart)
                goto 34      ! Storage space has been found, stop searching
              endif
33            continue       ! end search loop
            if (ipart.gt.maxpart) goto 996
34          minpart=ipart+1
30          continue         ! end loop over particles to be released
          endif
10        continue


      return

996   continue
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE RELEASEPARTICLES: ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - TOTAL NUMBER OF PARTICLES REQUIRED  ####'
      write(*,*) '#### EXCEEDS THE MAXIMUM ALLOWED NUMBER. REDUCE  ####'
      write(*,*) '#### EITHER NUMBER OF PARTICLES PER RELEASE POINT####'
      write(*,*) '#### OR REDUCE NUMBER OF RELEASE POINTS.         ####'
      write(*,*) '#####################################################'
      stop

      end
