      subroutine init_domainfill()
C
********************************************************************************
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subroutine init_domainfill.   *
*            The computational grid is the WRF x-y grid rather than lat-lon.   *
*                                                                              *
* Initializes particles equally distributed over the first release location    *
* specified in file RELEASES. This box is assumed to be the domain for doing   *
* domain-filling trajectory calculations.                                      *
* All particles carry the same amount of mass which alltogether comprises the  *
* mass of air within the box.                                                  *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     15 October 2002                                                          *
*                                                                              *
*    26 Oct 2005, R. Easter - changes for gridarea                             *
*                             associated with WRF horizontal grid.             *
*                             Also calc. true ylat for pv stuff.               *
*    11 Nov 2005, R. Easter - fixed error involving xy to latlon               *
*                                                                              *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
*                                                                              *
* numparticlecount    consecutively counts the number of particles released    *
* nx_we(2)       grid indices for western and eastern boundary of domain-      *
*                filling trajectory calculations                               *
* ny_sn(2)       grid indices for southern and northern boundary of domain-    *
*                filling trajectory calculations                               *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      integer j,ix,jy,kz,ncolumn,numparttot,idummy

c     real gridarea(0:nymax-1),pp(nzmax),ylat,ylatp,ylatm,hzone,ran1
      real gridarea(0:nymax-1),pp(nzmax),ylat,                  ran1
c     real cosfactm,cosfactp,pih,deltacol,dz1,dz2,dz,pnew,fractus
      real                   pih,deltacol,dz1,dz2,dz,pnew,fractus
      real xlon

      parameter(pih=pi/180.)
      real colmass(0:nxmax-1,0:nymax-1),colmasstotal,zposition

      integer ixm,ixp,jym,jyp,indzm,indzp,in,indzh,i,jj
      real pvpart,ddx,ddy,rddx,rddy,p1,p2,p3,p4,y1(2)

      data idummy/-11/


C Determine the release region (only full grid cells), over which particles
C shall be initialized
C Use 2 fields for west/east and south/north boundary
***************************************************************************

      nx_we(1)=max(int(xpoint1(1)),0)
      nx_we(2)=min((int(xpoint2(1))+1),nxmin1)
      ny_sn(1)=max(int(ypoint1(1)),0)
      ny_sn(2)=min((int(ypoint2(1))+1),nymin1)

C For global simulations (both global wind data and global domain-filling),
C set a switch, such that no boundary conditions are used
***************************************************************************
      if (xglobal.and.sglobal.and.nglobal) then
        if ((nx_we(1).eq.0).and.(nx_we(2).eq.nxmin1).and.
     +  (ny_sn(1).eq.0).and.(ny_sn(2).eq.nymin1)) then
          gdomainfill=.true.
        else
          gdomainfill=.false.
        endif
      endif

C Do not release particles twice (i.e., not at both in the leftmost and rightmost
C grid cell) for a global domain
*********************************************************************************
      if (xglobal) nx_we(2)=min(nx_we(2),nx-2)


C Calculate area of grid cell with formula M=2*pi*R*h*dx/360,
C see Netz, Formeln der Mathematik, 5. Auflage (1983), p.90
*************************************************************

      do 10 jy=ny_sn(1),ny_sn(2)      ! loop about latitudes

c        ylat=ylat0+float(jy)*dy
c        ylatp=ylat+0.5*dy
c        ylatm=ylat-0.5*dy
c        if ((ylatm.lt.0).and.(ylatp.gt.0.)) then
c          hzone=1./dyconst
c        else
c          cosfactp=cos(ylatp*pih)*r_earth
c          cosfactm=cos(ylatm*pih)*r_earth
c          if (cosfactp.lt.cosfactm) then
c            hzone=sqrt(r_earth**2-cosfactp**2)-
c     +      sqrt(r_earth**2-cosfactm**2)
c          else
c            hzone=sqrt(r_earth**2-cosfactm**2)-
c     +      sqrt(r_earth**2-cosfactp**2)
c          endif
c        endif
c10      gridarea(jy)=2.*pi*r_earth*hzone*dx/360.

c for FLEXPART_WRF, dx & dy are in meters, and no cos(lat) is needed
c ??? should maybe include map factor here ???
10      gridarea(jy)=dx*dy

C Do the same for the south pole

      if (sglobal) then
         write(*,*)
         write(*,*) '*** stopping in init_domainfill ***'
         write(*,*) '    the s-pole code section should not be active'
         write(*,*)
c        ylat=ylat0
c        ylatp=ylat+0.5*dy
c        ylatm=ylat
c        cosfactm=0.
c        cosfactp=cos(ylatp*pih)*r_earth
c        hzone=sqrt(r_earth**2-cosfactm**2)-
c     +  sqrt(r_earth**2-cosfactp**2)
c        gridarea(0)=2.*pi*r_earth*hzone*dx/360.
      endif

C Do the same for the north pole

      if (nglobal) then
         write(*,*)
         write(*,*) '*** stopping in init_domainfill ***'
         write(*,*) '    the s-pole code section should not be active'
         write(*,*)
c        ylat=ylat0+float(nymin1)*dy
c        ylatp=ylat
c        ylatm=ylat-0.5*dy
c        cosfactp=0.
c        cosfactm=cos(ylatm*pih)*r_earth
c        hzone=sqrt(r_earth**2-cosfactp**2)-
c     +  sqrt(r_earth**2-cosfactm**2)
c        gridarea(nymin1)=2.*pi*r_earth*hzone*dx/360.
      endif


C Calculate total mass of each grid column and of the whole atmosphere
**********************************************************************

      colmasstotal=0.
      do 20 jy=ny_sn(1),ny_sn(2)          ! loop about latitudes
        do 20 ix=nx_we(1),nx_we(2)      ! loop about longitudes
          pp(1)=rho(ix,jy,1,1)*r_air*tt(ix,jy,1,1)
          pp(nz)=rho(ix,jy,nz,1)*r_air*tt(ix,jy,nz,1)
          colmass(ix,jy)=(pp(1)-pp(nz))/ga*gridarea(jy)
20        colmasstotal=colmasstotal+colmass(ix,jy)

               write(*,*) 'Atm. mass: ',colmasstotal


      if (ipin.eq.0) numpart=0

C Determine the particle positions
**********************************

      numparttot=0
      numcolumn=0
      do 30 jy=ny_sn(1),ny_sn(2)      ! loop about latitudes
c       ylat=ylat0+float(jy)*dy
        do 30 ix=nx_we(1),nx_we(2)      ! loop about longitudes

c for FLEXPART_WRF, x & y coords are in meters.
c In the "do 30" loop, ylat is only needed for pv calcs.
          call xyindex_to_ll_wrf( 0, float(ix), float(jy), xlon, ylat )

          ncolumn=nint(0.999*float(npart(1))*colmass(ix,jy)/
     +    colmasstotal)
          if (ncolumn.eq.0) goto 30
          if (ncolumn.gt.numcolumn) numcolumn=ncolumn

C Calculate pressure at the altitudes of model surfaces, using the air density
C information, which is stored as a 3-d field
******************************************************************************

          do 31 kz=1,nz 
31          pp(kz)=rho(ix,jy,kz,1)*r_air*tt(ix,jy,kz,1)


          deltacol=(pp(1)-pp(nz))/float(ncolumn)
          pnew=pp(1)+deltacol/2.
          jj=0
          do 32 j=1,ncolumn
            jj=jj+1


C For columns with many particles (i.e. around the equator), distribute
C the particles equally, for columns with few particles (i.e. around the
C poles), distribute the particles randomly
************************************************************************


            if (ncolumn.gt.20) then
              pnew=pnew-deltacol
            else
              pnew=pp(1)-ran1(idummy)*(pp(1)-pp(nz))
            endif

            do 33 kz=1,nz-1
              if ((pp(kz).ge.pnew).and.(pp(kz+1).lt.pnew)) then
                dz1=pp(kz)-pnew
                dz2=pnew-pp(kz+1)
                dz=1./(dz1+dz2)

C Assign particle position
**************************
C Do the following steps only if particles are not read in from previous model run
**********************************************************************************
                if (ipin.eq.0) then
                  xtra1(numpart+jj)=float(ix)-0.5+ran1(idummy)
                  if (ix.eq.0) xtra1(numpart+jj)=ran1(idummy)
                  if (ix.eq.nxmin1) xtra1(numpart+jj)=
     +            float(nxmin1)-ran1(idummy)
                  ytra1(numpart+jj)=float(jy)-0.5+ran1(idummy)
                  ztra1(numpart+jj)=(height(kz)*dz2+height(kz+1)*dz1)*dz
                  if (ztra1(numpart+jj).gt.height(nz)-0.5)
     +            ztra1(numpart+jj)=height(nz)-0.5


C Interpolate PV to the particle position
*****************************************
                  ixm=int(xtra1(numpart+jj))
                  jym=int(ytra1(numpart+jj))
                  ixp=ixm+1
                  jyp=jym+1
                  ddx=xtra1(numpart+jj)-float(ixm)
                  ddy=ytra1(numpart+jj)-float(jym)
                  rddx=1.-ddx
                  rddy=1.-ddy
                  p1=rddx*rddy
                  p2=ddx*rddy
                  p3=rddx*ddy
                  p4=ddx*ddy
                  do 5 i=2,nz
                    if (height(i).gt.ztra1(numpart+jj)) then
                      indzm=i-1
                      indzp=i
                      goto 6
                    endif
5                   continue
6                 continue
                  dz1=ztra1(numpart+jj)-height(indzm)
                  dz2=height(indzp)-ztra1(numpart+jj)
                  dz=1./(dz1+dz2)
                  do 60 in=1,2     
                    indzh=indzm+in-1
60                  y1(in)=p1*pv(ixm,jym,indzh,1)
     +                    +p2*pv(ixp,jym,indzh,1)
     +                    +p3*pv(ixm,jyp,indzh,1)
     +                    +p4*pv(ixp,jyp,indzh,1)
                  pvpart=(dz2*y1(1)+dz1*y1(2))*dz
                  if (ylat.lt.0.) pvpart=-1.*pvpart


C For domain-filling option 2 (stratospheric O3), do the rest only in the stratosphere
**************************************************************************************

                  if (((ztra1(numpart+jj).gt.3000.).and.
     +            (pvpart.gt.pvcrit)).or.(mdomainfill.eq.1)) then

C Assign certain properties to the particle
*******************************************
                    nclass(numpart+jj)=min(int(ran1(idummy)*
     +              float(nclassunc))+1,nclassunc)
                    numparticlecount=numparticlecount+1
                    npoint(numpart+jj)=numparticlecount
                    idt(numpart+jj)=mintime
                    itra1(numpart+jj)=0
                    itramem(numpart+jj)=0
                    itrasplit(numpart+jj)=itra1(numpart+jj)+ldirect*
     +              itsplit
                    xmass1(numpart+jj,1)=colmass(ix,jy)/float(ncolumn)
                    if (mdomainfill.eq.2) xmass1(numpart+jj,1)=
     +             xmass1(numpart+jj,1)*pvpart*48./29.*ozonescale/10.**9
                  else
                    jj=jj-1
                  endif
                endif
              endif
33            continue
32          continue
          numparttot=numparttot+ncolumn
          if (ipin.eq.0) numpart=numpart+jj
30        continue

               
C Check whether numpart is really smaller than maxpart
******************************************************

      if (numpart.gt.maxpart) then
        write(*,*) 'numpart too large: change source in init_atm_mass.f'
        write(*,*) 'numpart: ',numpart,' maxpart: ',maxpart
      endif


      xmassperparticle=colmasstotal/float(numparttot)


C Make sure that all particles are within domain
************************************************

      do 40 j=1,numpart
        if ((xtra1(j).lt.0.).or.(xtra1(j).ge.float(nxmin1)).or.
     +  (ytra1(j).lt.0.).or.(ytra1(j).ge.float(nymin1))) then
          itra1(j)=-999999999
        endif
40      continue




C For boundary conditions, we need fewer particle release heights per column,
C because otherwise it takes too long until enough mass has accumulated to
C release a particle at the boundary (would take dx/u seconds), leading to
C relatively large position errors of the order of one grid distance.
C It's better to release fewer particles per column, but to do so more often.
C Thus, use on the order of nz starting heights per column.
C We thus repeat the above to determine fewer starting heights, that are
C used furtheron in subroutine boundcond_domainfill.f.
*****************************************************************************

      fractus=float(numcolumn)/float(nz)
      write(*,*) 'Total number of particles at model start: ',numpart
      write(*,*) 'Maximum number of particles per column: ',numcolumn
      write(*,*) 'If ',fractus,' <1, better use more particles'
      fractus=sqrt(max(fractus,1.))/2.

      do 80 jy=ny_sn(1),ny_sn(2)      ! loop about latitudes
        do 80 ix=nx_we(1),nx_we(2)      ! loop about longitudes
          ncolumn=nint(0.999/fractus*float(npart(1))*colmass(ix,jy)
     +    /colmasstotal)
          if (ncolumn.gt.maxcolumn) stop 'maxcolumn too small'
          if (ncolumn.eq.0) goto 80


C Memorize how many particles per column shall be used for all boundaries
C This is further used in subroutine boundcond_domainfill.f
C Use 2 fields for west/east and south/north boundary
*************************************************************************

          if (ix.eq.nx_we(1)) numcolumn_we(1,jy)=ncolumn
          if (ix.eq.nx_we(2)) numcolumn_we(2,jy)=ncolumn
          if (jy.eq.ny_sn(1)) numcolumn_sn(1,ix)=ncolumn
          if (jy.eq.ny_sn(2)) numcolumn_sn(2,ix)=ncolumn

C Calculate pressure at the altitudes of model surfaces, using the air density
C information, which is stored as a 3-d field
******************************************************************************

          do 81 kz=1,nz 
81          pp(kz)=rho(ix,jy,kz,1)*r_air*tt(ix,jy,kz,1)

C Determine the reference starting altitudes
********************************************

          deltacol=(pp(1)-pp(nz))/float(ncolumn)
          pnew=pp(1)+deltacol/2.
          do 82 j=1,ncolumn
            pnew=pnew-deltacol
            do 83 kz=1,nz-1
              if ((pp(kz).ge.pnew).and.(pp(kz+1).lt.pnew)) then
                dz1=pp(kz)-pnew
                dz2=pnew-pp(kz+1)
                dz=1./(dz1+dz2)
                zposition=(height(kz)*dz2+height(kz+1)*dz1)*dz
               if (zposition.gt.height(nz)-0.5) zposition=height(nz)-0.5

C Memorize vertical positions where particles are introduced
C This is further used in subroutine boundcond_domainfill.f
************************************************************

                if (ix.eq.nx_we(1)) zcolumn_we(1,jy,j)=zposition
                if (ix.eq.nx_we(2)) zcolumn_we(2,jy,j)=zposition
                if (jy.eq.ny_sn(1)) zcolumn_sn(1,ix,j)=zposition
                if (jy.eq.ny_sn(2)) zcolumn_sn(2,ix,j)=zposition

C Initialize mass that has accumulated at boundary to zero
**********************************************************

                acc_mass_we(1,jy,j)=0.
                acc_mass_we(2,jy,j)=0.
                acc_mass_sn(1,jy,j)=0.
                acc_mass_sn(2,jy,j)=0.
              endif
83            continue
82          continue
80        continue

C If particles shall be read in to continue an existing run,
C then the accumulated masses at the domain boundaries must be read in, too.
C This overrides any previous calculations.
****************************************************************************

      if (ipin.eq.1) then
        open(unitboundcond,file=path(2)(1:len(2))//'boundcond.bin',
     +  form='unformatted')
        read(unitboundcond) numcolumn_we,numcolumn_sn,
     +  zcolumn_we,zcolumn_sn,acc_mass_we,acc_mass_sn
        close(unitboundcond)
      endif




      end
