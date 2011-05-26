      subroutine boundcond_domainfill(itime,loutend)
C                                       i      i
********************************************************************************
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subr. boundcond_domainfill.   *
*            The computational grid is the WRF x-y grid rather than lat-lon.   *
*                                                                              *
* Particles are created by this subroutine continuously throughout the         *
* simulation at the boundaries of the domain-filling box.                      *
* All particles carry the same amount of mass which alltogether comprises the  *
* mass of air within the box, which remains (more or less) constant.           *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     16 October 2002                                                          *
*                                                                              *
*    26 Oct 2005, R. Easter - changes to calc. of boundarea                    *
*                             associated with WRF horizontal grid.             *
*                             Also need to get true ylat for pv calcs.         *
*    11 Nov 2005, R. Easter - fixed error involving xy to latlong              *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
*                                                                              *
* nx_we(2)       grid indices for western and eastern boundary of domain-      *
*                filling trajectory calculations                               *
* ny_sn(2)       grid indices for southern and northern boundary of domain-    *
*                filling trajectory calculations                               *
*                                                                              *
********************************************************************************


      include 'includepar'
      include 'includecom'

      real dz,dz1,dz2,ran1,dt1,dt2,dtt,ylat,xm,cosfact,accmasst
      real dumx,dumy,xlon
      integer itime,in,indz,indzp,i,idummy,loutend
      integer j,k,ix,jy,m,indzh,indexh,minpart,ipart,mmass
      integer numactiveparticles

      real windl(2),rhol(2)
      real windhl(2),rhohl(2)
      real windx,rhox
      real deltaz,boundarea,fluxofmass

      integer ixm,ixp,jym,jyp,indzm,mm
      real pvpart,ddx,ddy,rddx,rddy,p1,p2,p3,p4,y1(2),yh1(2)

      save idummy
      data idummy/-11/


C If domain-filling is global, no boundary conditions are needed
****************************************************************

      if (gdomainfill) return
    
      accmasst=0.
      numactiveparticles=0

C Terminate trajectories that have left the domain, if domain-filling
C trajectory calculation domain is not global
*********************************************************************

      do 85 i=1,numpart
        if (itra1(i).eq.itime) then
          if ((ytra1(i).gt.float(ny_sn(2))).or.
     +    (ytra1(i).lt.float(ny_sn(1)))) itra1(i)=-999999999
          if (((.not.xglobal).or.(nx_we(2).ne.(nx-2))).and.
     +    ((xtra1(i).lt.float(nx_we(1))).or.
     +    (xtra1(i).gt.float(nx_we(2))))) itra1(i)=-999999999
        endif
        if (itra1(i).ne.-999999999) numactiveparticles=
     +  numactiveparticles+1
85      continue


C Determine auxiliary variables for time interpolation
******************************************************

      dt1=float(itime-memtime(1))
      dt2=float(memtime(2)-itime)
      dtt=1./(dt1+dt2)

C Initialize auxiliary variable used to search for vacant storage space
***********************************************************************

      minpart=1

****************************************
C Western and eastern boundary condition
****************************************

C Loop from south to north
**************************

      do 70 jy=ny_sn(1),ny_sn(2)

C Loop over western (index 1) and eastern (index 2) boundary
************************************************************

        do 70 k=1,2

c for FLEXPART_WRF, x & y coords are in meters.
c In the "do 70" loop, ylat is only needed for for pv calcs,
c     "if (ylat.lt.0.) pvpart=-1.*pvpart"
c Note: in the FLEXPART_ECMWF code, ylat was not defined 
c     in the "do 70" loop (a bug).
          dumx=float(nx_we(k))
          dumy=float(jy)
c Are these dumx,dumy correct ???
          call xyindex_to_ll_wrf( 0, dumx, dumy, xlon, ylat )

C Loop over all release locations in a column
*********************************************

          do 70 j=1,numcolumn_we(k,jy)

C Determine, for each release location, the area of the corresponding boundary
******************************************************************************

            if (j.eq.1) then
              deltaz=(zcolumn_we(k,jy,2)+zcolumn_we(k,jy,1))/2.
            else if (j.eq.numcolumn_we(k,jy)) then
c             deltaz=height(nz)-(zcolumn_we(k,jy,j-1)+
c    +        zcolumn_we(k,jy,j))/2.
C In order to avoid taking a very high column for very many particles,
C use the deltaz from one particle below instead
              deltaz=(zcolumn_we(k,jy,j)-zcolumn_we(k,jy,j-2))/2.
            else
              deltaz=(zcolumn_we(k,jy,j+1)-zcolumn_we(k,jy,j-1))/2.
            endif

c for FLEXPART_ECMWF, dy is in degrees-lat, and 111198.5 converts
c   from degrees-latitude to m 
c for FLEXPART_WRF, dy is in meters
c           if ((jy.eq.ny_sn(1)).or.(jy.eq.ny_sn(2))) then
c             boundarea=deltaz*111198.5/2.*dy
c           else
c             boundarea=deltaz*111198.5*dy
c           endif
            if ((jy.eq.ny_sn(1)).or.(jy.eq.ny_sn(2))) then
              boundarea=deltaz/2.*dy
            else
              boundarea=deltaz*dy
            endif


C Interpolate the wind velocity and density to the release location
*******************************************************************

C Determine the model level below the release position
******************************************************

            do 5 i=2,nz
              if (height(i).gt.zcolumn_we(k,jy,j)) then
                indz=i-1
                indzp=i
                goto 6
              endif
5             continue
6           continue

C Vertical distance to the level below and above current position
*****************************************************************

            dz1=zcolumn_we(k,jy,j)-height(indz)
            dz2=height(indzp)-zcolumn_we(k,jy,j)
            dz=1./(dz1+dz2)

C Vertical and temporal interpolation
*************************************

            do 31 m=1,2
              indexh=memind(m)
              do 32 in=1,2
                indzh=indz+in-1
                windl(in)=uu(nx_we(k),jy,indzh,indexh)
                rhol(in)=rho(nx_we(k),jy,indzh,indexh)
32              continue

              windhl(m)=(dz2*windl(1)+dz1*windl(2))*dz
              rhohl(m)=(dz2*rhol(1)+dz1*rhol(2))*dz
31            continue

            windx=(windhl(1)*dt2+windhl(2)*dt1)*dtt
            rhox=(rhohl(1)*dt2+rhohl(2)*dt1)*dtt

C Calculate mass flux
*********************

            fluxofmass=windx*rhox*boundarea*float(lsynctime)


C If the mass flux is directed into the domain, add it to previous mass fluxes;
C if it is out of the domain, set accumulated mass flux to zero
*******************************************************************************

            if (k.eq.1) then
              if (fluxofmass.ge.0.) then
                acc_mass_we(k,jy,j)=acc_mass_we(k,jy,j)+fluxofmass
              else
                acc_mass_we(k,jy,j)=0.
              endif
            else
              if (fluxofmass.le.0.) then
                acc_mass_we(k,jy,j)=acc_mass_we(k,jy,j)+abs(fluxofmass)
              else
                acc_mass_we(k,jy,j)=0.
              endif
            endif
            accmasst=accmasst+acc_mass_we(k,jy,j)

C If the accumulated mass exceeds half the mass that each particle shall carry,
C one (or more) particle(s) is (are) released and the accumulated mass is
C reduced by the mass of this (these) particle(s)
*******************************************************************************

            if (acc_mass_we(k,jy,j).ge.xmassperparticle/2.) then
              mmass=int((acc_mass_we(k,jy,j)+xmassperparticle/2.)/
     +        xmassperparticle)
              acc_mass_we(k,jy,j)=acc_mass_we(k,jy,j)-
     +        float(mmass)*xmassperparticle
            else
              mmass=0
            endif

            do 71 m=1,mmass
              do 72 ipart=minpart,maxpart

C If a vacant storage space is found, attribute everything to this array element
********************************************************************************

                if (itra1(ipart).ne.itime) then

C Assign particle positions
***************************

                  xtra1(ipart)=float(nx_we(k))
                  if (jy.eq.ny_sn(1)) then
                    ytra1(ipart)=float(jy)+0.5*ran1(idummy)
                  else if (jy.eq.ny_sn(2)) then
                    ytra1(ipart)=float(jy)-0.5*ran1(idummy)
                  else
                    ytra1(ipart)=float(jy)+(ran1(idummy)-.5)
                  endif
                  if (j.eq.1) then
                    ztra1(ipart)=zcolumn_we(k,jy,1)+(zcolumn_we(k,jy,2)-
     +              zcolumn_we(k,jy,1))/4.
                  else if (j.eq.numcolumn_we(k,jy)) then
                    ztra1(ipart)=(2.*zcolumn_we(k,jy,j)+
     +              zcolumn_we(k,jy,j-1)+height(nz))/4.
                  else
                    ztra1(ipart)=zcolumn_we(k,jy,j-1)+ran1(idummy)*
     +              (zcolumn_we(k,jy,j+1)-zcolumn_we(k,jy,j-1))
                  endif

C Interpolate PV to the particle position
*****************************************
                  ixm=int(xtra1(ipart))
                  jym=int(ytra1(ipart))
                  ixp=ixm+1
                  jyp=jym+1
                  ddx=xtra1(ipart)-float(ixm)
                  ddy=ytra1(ipart)-float(jym)
                  rddx=1.-ddx
                  rddy=1.-ddy
                  p1=rddx*rddy
                  p2=ddx*rddy
                  p3=rddx*ddy
                  p4=ddx*ddy
                  do 25 i=2,nz
                    if (height(i).gt.ztra1(ipart)) then
                      indzm=i-1
                      indzp=i
                      goto 26
                    endif
25                  continue
26                continue
                  dz1=ztra1(ipart)-height(indzm)
                  dz2=height(indzp)-ztra1(ipart)
                  dz=1./(dz1+dz2)
                  do 61 mm=1,2
                    indexh=memind(mm)
                    do 60 in=1,2
                      indzh=indzm+in-1
60                    y1(in)=p1*pv(ixm,jym,indzh,indexh)
     +                      +p2*pv(ixp,jym,indzh,indexh)
     +                      +p3*pv(ixm,jyp,indzh,indexh)
     +                      +p4*pv(ixp,jyp,indzh,indexh)
61                yh1(mm)=(dz2*y1(1)+dz1*y1(2))*dz
                  pvpart=(yh1(1)*dt2+yh1(2)*dt1)*dtt
                  if (ylat.lt.0.) pvpart=-1.*pvpart


C For domain-filling option 2 (stratospheric O3), do the rest only in the stratosphere
**************************************************************************************

                  if (((ztra1(ipart).gt.3000.).and.
     +            (pvpart.gt.pvcrit)).or.(mdomainfill.eq.1)) then
                    nclass(ipart)=min(int(ran1(idummy)*
     +              float(nclassunc))+1,nclassunc)
                    numactiveparticles=numactiveparticles+1
                    numparticlecount=numparticlecount+1
                    npoint(ipart)=numparticlecount
                    idt(ipart)=mintime
                    itra1(ipart)=itime
                    itramem(ipart)=itra1(ipart)
                    itrasplit(ipart)=itra1(ipart)+ldirect*itsplit
                    xmass1(ipart,1)=xmassperparticle
                    if (mdomainfill.eq.2) xmass1(ipart,1)=
     +              xmass1(ipart,1)*pvpart*48./29.*ozonescale/10.**9
                  else
                    goto 71
                  endif


C Increase numpart, if necessary
********************************

                  numpart=max(numpart,ipart)
                  goto 73      ! Storage space has been found, stop searching
                endif
72              continue
              if (ipart.gt.maxpart)
     +        stop 'boundcond_domainfill.f: too many particles required'
73            minpart=ipart+1
71            continue


70          continue


******************************************
C Southern and northern boundary condition
******************************************

C Loop from west to east
************************

      do 170 ix=nx_we(1),nx_we(2)

C Loop over southern (index 1) and northern (index 2) boundary
**************************************************************

        do 170 k=1,2

c for FLEXPART_WRF, x & y coords are in meters.
c         ylat=ylat0+float(ny_sn(k))*dy
c         cosfact=cos(ylat*pi180)
c In the "do 170" loop, ylat is only needed for for pv calcs,
c    "if (ylat.lt.0.) pvpart=-1.*pvpart"
          dumx=float(ix)
          dumy=float(ny_sn(k))
c Are these dumx,dumy correct ???
          call xyindex_to_ll_wrf( 0, dumx, dumy, xlon, ylat )

C Loop over all release locations in a column
*********************************************

          do 170 j=1,numcolumn_sn(k,ix)

C Determine, for each release location, the area of the corresponding boundary
******************************************************************************

            if (j.eq.1) then
              deltaz=(zcolumn_sn(k,ix,2)+zcolumn_sn(k,ix,1))/2.
            else if (j.eq.numcolumn_sn(k,ix)) then
c             deltaz=height(nz)-(zcolumn_sn(k,ix,j-1)+
c    +        zcolumn_sn(k,ix,j))/2.
C In order to avoid taking a very high column for very many particles,
C use the deltaz from one particle below instead
              deltaz=(zcolumn_sn(k,ix,j)-zcolumn_sn(k,ix,j-2))/2.
            else
              deltaz=(zcolumn_sn(k,ix,j+1)-zcolumn_sn(k,ix,j-1))/2.
            endif

c for FLEXPART_ECMWF, dx is in degrees-long, and 111198.5*cosfact converts
c   from degrees-longitude to m 
c for FLEXPART_WRF, dx is in meters
c           if ((ix.eq.nx_we(1)).or.(ix.eq.nx_we(2))) then
c             boundarea=deltaz*111198.5/2.*cosfact*dx
c           else
c             boundarea=deltaz*111198.5*cosfact*dx
c           endif
            if ((ix.eq.nx_we(1)).or.(ix.eq.nx_we(2))) then
              boundarea=deltaz/2.*dx
            else
              boundarea=deltaz*dx
            endif


C Interpolate the wind velocity and density to the release location
*******************************************************************

C Determine the model level below the release position
******************************************************

            do 15 i=2,nz
              if (height(i).gt.zcolumn_sn(k,ix,j)) then
                indz=i-1
                indzp=i
                goto 16
              endif
15            continue
16          continue

C Vertical distance to the level below and above current position
*****************************************************************

            dz1=zcolumn_sn(k,ix,j)-height(indz)
            dz2=height(indzp)-zcolumn_sn(k,ix,j)
            dz=1./(dz1+dz2)

C Vertical and temporal interpolation
*************************************

            do 131 m=1,2
              indexh=memind(m)
              do 132 in=1,2
                indzh=indz+in-1
                windl(in)=vv(ix,ny_sn(k),indzh,indexh)
                rhol(in)=rho(ix,ny_sn(k),indzh,indexh)
132             continue

              windhl(m)=(dz2*windl(1)+dz1*windl(2))*dz
              rhohl(m)=(dz2*rhol(1)+dz1*rhol(2))*dz
131           continue

            windx=(windhl(1)*dt2+windhl(2)*dt1)*dtt
            rhox=(rhohl(1)*dt2+rhohl(2)*dt1)*dtt

C Calculate mass flux
*********************

            fluxofmass=windx*rhox*boundarea*float(lsynctime)

C If the mass flux is directed into the domain, add it to previous mass fluxes;
C if it is out of the domain, set accumulated mass flux to zero
*******************************************************************************

            if (k.eq.1) then
              if (fluxofmass.ge.0.) then
                acc_mass_sn(k,ix,j)=acc_mass_sn(k,ix,j)+fluxofmass
              else
                acc_mass_sn(k,ix,j)=0.
              endif
            else
              if (fluxofmass.le.0.) then
                acc_mass_sn(k,ix,j)=acc_mass_sn(k,ix,j)+abs(fluxofmass)
              else
                acc_mass_sn(k,ix,j)=0.
              endif
            endif
            accmasst=accmasst+acc_mass_sn(k,ix,j)

C If the accumulated mass exceeds half the mass that each particle shall carry,
C one (or more) particle(s) is (are) released and the accumulated mass is
C reduced by the mass of this (these) particle(s)
*******************************************************************************

            if (acc_mass_sn(k,ix,j).ge.xmassperparticle/2.) then
              mmass=int((acc_mass_sn(k,ix,j)+xmassperparticle/2.)/
     +        xmassperparticle)
              acc_mass_sn(k,ix,j)=acc_mass_sn(k,ix,j)-
     +        float(mmass)*xmassperparticle
            else
              mmass=0
            endif

            do 171 m=1,mmass
              do 172 ipart=minpart,maxpart

C If a vacant storage space is found, attribute everything to this array element
********************************************************************************

                if (itra1(ipart).ne.itime) then

C Assign particle positions
***************************

                  ytra1(ipart)=float(ny_sn(k))
                  if (ix.eq.nx_we(1)) then
                    xtra1(ipart)=float(ix)+0.5*ran1(idummy)
                  else if (ix.eq.nx_we(2)) then
                    xtra1(ipart)=float(ix)-0.5*ran1(idummy)
                  else
                    xtra1(ipart)=float(ix)+(ran1(idummy)-.5)
                  endif
                  if (j.eq.1) then
                    ztra1(ipart)=zcolumn_sn(k,ix,1)+(zcolumn_sn(k,ix,2)-
     +              zcolumn_sn(k,ix,1))/4.
                  else if (j.eq.numcolumn_sn(k,ix)) then
                    ztra1(ipart)=(2.*zcolumn_sn(k,ix,j)+
     +              zcolumn_sn(k,ix,j-1)+height(nz))/4.
                  else
                    ztra1(ipart)=zcolumn_sn(k,ix,j-1)+ran1(idummy)*
     +              (zcolumn_sn(k,ix,j+1)-zcolumn_sn(k,ix,j-1))
                  endif


C Interpolate PV to the particle position
*****************************************
                  ixm=int(xtra1(ipart))
                  jym=int(ytra1(ipart))
                  ixp=ixm+1
                  jyp=jym+1
                  ddx=xtra1(ipart)-float(ixm)
                  ddy=ytra1(ipart)-float(jym)
                  rddx=1.-ddx
                  rddy=1.-ddy
                  p1=rddx*rddy
                  p2=ddx*rddy
                  p3=rddx*ddy
                  p4=ddx*ddy
                  do 125 i=2,nz
                    if (height(i).gt.ztra1(ipart)) then
                      indzm=i-1
                      indzp=i
                      goto 126
                    endif
125                 continue
126               continue
                  dz1=ztra1(ipart)-height(indzm)
                  dz2=height(indzp)-ztra1(ipart)
                  dz=1./(dz1+dz2)
                  do 161 mm=1,2
                    indexh=memind(mm)
                    do 160 in=1,2
                      indzh=indzm+in-1
160                   y1(in)=p1*pv(ixm,jym,indzh,indexh)
     +                      +p2*pv(ixp,jym,indzh,indexh)
     +                      +p3*pv(ixm,jyp,indzh,indexh)
     +                      +p4*pv(ixp,jyp,indzh,indexh)
161               yh1(mm)=(dz2*y1(1)+dz1*y1(2))*dz
                  pvpart=(yh1(1)*dt2+yh1(2)*dt1)*dtt
                  if (ylat.lt.0.) pvpart=-1.*pvpart


C For domain-filling option 2 (stratospheric O3), do the rest only in the stratosphere
**************************************************************************************

                  if (((ztra1(ipart).gt.3000.).and.
     +            (pvpart.gt.pvcrit)).or.(mdomainfill.eq.1)) then
                    nclass(ipart)=min(int(ran1(idummy)*
     +              float(nclassunc))+1,nclassunc)
                    numactiveparticles=numactiveparticles+1
                    numparticlecount=numparticlecount+1
                    npoint(ipart)=numparticlecount
                    idt(ipart)=mintime
                    itra1(ipart)=itime
                    itramem(ipart)=itra1(ipart)
                    itrasplit(ipart)=itra1(ipart)+ldirect*itsplit
                    xmass1(ipart,1)=xmassperparticle
                    if (mdomainfill.eq.2) xmass1(ipart,1)=
     +              xmass1(ipart,1)*pvpart*48./29.*ozonescale/10.**9
                  else
                    goto 171
                  endif


C Increase numpart, if necessary
********************************
                  numpart=max(numpart,ipart)
                  goto 173      ! Storage space has been found, stop searching
                endif
172             continue
              if (ipart.gt.maxpart)
     +        stop 'boundcond_domainfill.f: too many particles required'
173           minpart=ipart+1
171           continue


170         continue


      xm=0.
      do 55 i=1,numpart
        if (itra1(i).eq.itime) xm=xm+xmass1(i,1)
55      continue

c     write(*,*) itime,numactiveparticles,numparticlecount,numpart,
c    +xm,accmasst,xm+accmasst


C If particles shall be dumped, then accumulated masses at the domain boundaries
C must be dumped, too, to be used for later runs
********************************************************************************

      if ((ipout.gt.0).and.(itime.eq.loutend)) then
        open(unitboundcond,file=path(2)(1:len(2))//'boundcond.bin',
     +  form='unformatted')
        write(unitboundcond) numcolumn_we,numcolumn_sn,
     +  zcolumn_we,zcolumn_sn,acc_mass_we,acc_mass_sn
        close(unitboundcond)
      endif

      end
