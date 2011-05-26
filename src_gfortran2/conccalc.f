      subroutine conccalc(itime,weight)
C                           i     i
********************************************************************************
*                                                                              *
*     Calculation of the concentrations on a regular grid using volume sampling*
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     24 May 1996                                                              *
*                                                                              *
*     April 2000: Update to calculate age spectra                              *
*                 Bug fix to avoid negative conc. at the domain boundaries,    *
*                 as suggested by Petra Seibert                                *
*                                                                              *
*     2 July 2002: re-order if-statements in order to optimize CPU time        *
*                                                                              *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* nspeciesdim     = nspec for forward runs, 1 for backward runs                *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      integer itime,itage,i,ix,jy,ixp,jyp,kz,k,n,nage,nspeciesdim
      integer il,ind,indz,indzp
      real rddx,rddy,p1,p2,p3,p4,dz1,dz2,dz
      real weight,hx,hy,hz,h,xd,yd,zd,factor,xkern,r2,c(maxspec),ddx,ddy
      real rhoprof(2),rhoi
      real hxmax,hymax,hzmax,xl,yl,wx,wy,w
      integer nspecpointer(maxspec)
      parameter(factor=.596831,hxmax=6.0,hymax=4.0,hzmax=150.)
c flexpart_wrf 06-nov-2005
      integer iconccalc_diagaa


c flexpart_wrf 06-nov-2005
      write(*,*) 'conccalc -- itime, weight =', itime, weight
      iconccalc_diagaa = 0


C For forward simulations, make a loop over the number of species;
C for backward simulations, there is only one species, but the output
C grid is dimensioned with numpoint (number of "backward release"
C (i.e., receptor) points),and each particle makes it's contribution only
C once to this output field
****************************************************************************

      if (ldirect.eq.1) then
        nspeciesdim=nspec
        do 27 k=1,nspec
27        nspecpointer(k)=k
      else
        nspeciesdim=1
      endif


      do 20 i=1,numpart
        if (itra1(i).ne.itime) goto 20

C Determine age class of the particle
        itage=abs(itra1(i)-itramem(i))
        do 32 nage=1,nageclass
          if (itage.lt.lage(nage)) goto 33
32        continue
33      continue


C For special runs, interpolate the air density to the particle position
*************************************************************************
************************************************************************
cAF IND_SOURCE switches between different units for concentrations at the source
cAf    NOTE that in backward simulations the release of particles takes place 
cAf    at the receptor and the sampling at the source.
cAf          1="mass" 
cAf          2="mass mixing ratio" 
cAf IND_RECEPTOR switches between different units for concentrations at the receptor
cAf          1="mass" 
cAf          2="mass mixing ratio" 

cAf switches for the conccalcfile:
cAF IND_SAMP =  0 : xmass * 1
cAf IND_SAMP = -1 : xmass / rho

cAf ind_samp is defined in readcommand.f

        if ( ind_samp .eq. -1 ) then

          ix=int(xtra1(i))
          jy=int(ytra1(i))
          ixp=ix+1
          jyp=jy+1
          ddx=xtra1(i)-float(ix)
          ddy=ytra1(i)-float(jy)
          rddx=1.-ddx
          rddy=1.-ddy
          p1=rddx*rddy
          p2=ddx*rddy
          p3=rddx*ddy
          p4=ddx*ddy

          do 5 il=2,nz
            if (height(il).gt.ztra1(i)) then
              indz=il-1
              indzp=il
              goto 6
            endif
5           continue
6         continue

          dz1=ztra1(i)-height(indz)
          dz2=height(indzp)-ztra1(i)
          dz=1./(dz1+dz2)

C Take density from 2nd wind field in memory (accurate enough, no time interpolation needed)
********************************************************************************************
          do 70 ind=indz,indzp
70          rhoprof(ind-indz+1)=p1*rho(ix ,jy ,ind,2)
     +                         +p2*rho(ixp,jy ,ind,2)
     +                         +p3*rho(ix ,jyp,ind,2)
     +                         +p4*rho(ixp,jyp,ind,2)
          rhoi=(dz1*rhoprof(2)+dz2*rhoprof(1))*dz
       elseif (ind_samp.eq.0) then
          rhoi = 1.
       endif


*****************************************************************************
C 1. Evaluate grid concentrations using a uniform kernel of bandwidths dx, dy
*****************************************************************************


C For backward simulations, look from which release point the particle comes from
C For domain-filling trajectory option, npoint contains a consecutive particle
C number, not the release point information. Therefore, nspecpointer is set to 1
C for the domain-filling option.
*********************************************************************************

        if (ldirect.ne.1) then 
          if (mdomainfill.eq.0) then
            nspecpointer(1)=npoint(i)
          else
            nspecpointer(1)=1
          endif
        endif
        
        do 30 kz=1,numzgrid                ! determine height of cell
          if (outheight(kz).gt.ztra1(i)) goto 21
30        continue
21      continue
        if (kz.le.numzgrid) then           ! inside output domain


*********************************
C Do everything for mother domain
*********************************

          xl=(xtra1(i)*dx+xoutshift)/dxout
          yl=(ytra1(i)*dy+youtshift)/dyout
          ix=int(xl)
          if (xl.lt.0.) ix=ix-1
          jy=int(yl)
          if (yl.lt.0.) jy=jy-1

c      if (i.eq.10000) write(*,*) itime,xtra1(i),ytra1(i),ztra1(i),xl,yl


C For particles aged less than 3 hours, attribute particle mass to grid cell
C it resides in rather than use the kernel, in order to avoid its smoothing effect.
C For older particles, use the uniform kernel.
C If a particle is close to the domain boundary, do not use the kernel either.
***********************************************************************************

          if ((itage.lt.10800).or.(xl.lt.0.5).or.(yl.lt.0.5).or.
     +    (xl.gt.float(numxgrid-1)-0.5).or.
     +    (yl.gt.float(numygrid-1)-0.5)) then             ! no kernel, direct attribution to grid cell
            if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgrid-1).and.
     +      (jy.le.numygrid-1)) then
              do 42 k=1,nspeciesdim
42              gridunc(ix,jy,kz,nspecpointer(k),nclass(i),nage)=
     +          gridunc(ix,jy,kz,nspecpointer(k),nclass(i),nage)+
     +          xmass1(i,k)/rhoi*weight
c flexpart_wrf 06-nov-2005
              if (iconccalc_diagaa .gt. 0)
     &        write(*,'(a,4i4,1p,e10.2,0p,3f8.2)') 
     &        'conccalc aa', i, ix, jy, kz,
     &        gridunc(ix,jy,kz,nspecpointer(1),nclass(i),nage),
     &        xtra1(i)*dx*1.0e-3, ytra1(i)*dx*1.0e-3, ztra1(i)*1.0e-3
            endif

          else                                 ! attribution via uniform kernel

            ddx=xl-float(ix)                   ! distance to left cell border
            ddy=yl-float(jy)                   ! distance to lower cell border
            if (ddx.gt.0.5) then
              ixp=ix+1
              wx=1.5-ddx
            else
              ixp=ix-1
              wx=0.5+ddx
            endif

            if (ddy.gt.0.5) then
              jyp=jy+1
              wy=1.5-ddy
            else
              jyp=jy-1
              wy=0.5+ddy
            endif


C Determine mass fractions for four grid points
***********************************************

            if ((ix.ge.0).and.(ix.le.numxgrid-1)) then
              if ((jy.ge.0).and.(jy.le.numygrid-1)) then
                w=wx*wy
                do 22 k=1,nspeciesdim
22                gridunc(ix,jy,kz,nspecpointer(k),nclass(i),nage)=
     +            gridunc(ix,jy,kz,nspecpointer(k),nclass(i),nage)+
     +            xmass1(i,k)/rhoi*weight*w
c flexpart_wrf 06-nov-2005
              if (iconccalc_diagaa .gt. 0)
     &        write(*,'(a,4i4,1p,e10.2,0p,3f8.2)') 
     &        'conccalc bb', i, ix, jy, kz,
     &        gridunc(ix,jy,kz,nspecpointer(1),nclass(i),nage),
     &        xtra1(i)*dx*1.0e-3, ytra1(i)*dx*1.0e-3, ztra1(i)*1.0e-3
              endif

              if ((jyp.ge.0).and.(jyp.le.numygrid-1)) then
                w=wx*(1.-wy)
                do 25 k=1,nspeciesdim
25                gridunc(ix,jyp,kz,nspecpointer(k),nclass(i),nage)=
     +            gridunc(ix,jyp,kz,nspecpointer(k),nclass(i),nage)+
     +            xmass1(i,k)/rhoi*weight*w
              endif
c flexpart_wrf 06-nov-2005
              if (iconccalc_diagaa .gt. 0)
     &        write(*,'(a,4i4,1p,e10.2,0p,3f8.2)') 
     &        'conccalc cc', i, ix, jy, kz,
     &        gridunc(ix,jy,kz,nspecpointer(1),nclass(i),nage),
     &        xtra1(i)*dx*1.0e-3, ytra1(i)*dx*1.0e-3, ztra1(i)*1.0e-3
            endif


            if ((ixp.ge.0).and.(ixp.le.numxgrid-1)) then
              if ((jyp.ge.0).and.(jyp.le.numygrid-1)) then
                w=(1.-wx)*(1.-wy)
                do 23 k=1,nspeciesdim
23                gridunc(ixp,jyp,kz,nspecpointer(k),nclass(i),nage)=
     +            gridunc(ixp,jyp,kz,nspecpointer(k),nclass(i),nage)+
     +            xmass1(i,k)/rhoi*weight*w
c flexpart_wrf 06-nov-2005
              if (iconccalc_diagaa .gt. 0)
     &        write(*,'(a,4i4,1p,e10.2,0p,3f8.2)') 
     &        'conccalc dd', i, ix, jy, kz,
     &        gridunc(ix,jy,kz,nspecpointer(1),nclass(i),nage),
     &        xtra1(i)*dx*1.0e-3, ytra1(i)*dx*1.0e-3, ztra1(i)*1.0e-3
              endif

              if ((jy.ge.0).and.(jy.le.numygrid-1)) then
                w=(1.-wx)*wy
                do 24 k=1,nspeciesdim
24                gridunc(ixp,jy,kz,nspecpointer(k),nclass(i),nage)=
     +            gridunc(ixp,jy,kz,nspecpointer(k),nclass(i),nage)+
     +            xmass1(i,k)/rhoi*weight*w
c flexpart_wrf 06-nov-2005
              if (iconccalc_diagaa .gt. 0)
     &        write(*,'(a,4i4,1p,e10.2,0p,3f8.2)') 
     &        'conccalc ee', i, ix, jy, kz,
     &        gridunc(ix,jy,kz,nspecpointer(1),nclass(i),nage),
     &        xtra1(i)*dx*1.0e-3, ytra1(i)*dx*1.0e-3, ztra1(i)*1.0e-3
              endif
            endif
          endif



*************************************
C Do everything for the nested domain
*************************************

          if (nested_output.eq.1) then
            xl=(xtra1(i)*dx+xoutshiftn)/dxoutn
            yl=(ytra1(i)*dy+youtshiftn)/dyoutn
            ix=int(xl)
            if (xl.lt.0.) ix=ix-1
            jy=int(yl)
            if (yl.lt.0.) jy=jy-1


C For particles aged less than 3 hours, attribute particle mass to grid cell
C it resides in rather than use the kernel, in order to avoid its smoothing effect.
C For older particles, use the uniform kernel.
C If a particle is close to the domain boundary, do not use the kernel either.
***********************************************************************************

            if ((itage.lt.10800).or.(xl.lt.0.5).or.(yl.lt.0.5).or.
     +      (xl.gt.float(numxgridn-1)-0.5).or.
     +      (yl.gt.float(numygridn-1)-0.5)) then             ! no kernel, direct attribution to grid cell
              if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgridn-1).and.
     +        (jy.le.numygridn-1)) then
                do 142 k=1,nspeciesdim
142               griduncn(ix,jy,kz,nspecpointer(k),nclass(i),nage)=
     +            griduncn(ix,jy,kz,nspecpointer(k),nclass(i),nage)+
     +            xmass1(i,k)/rhoi*weight
              endif

            else                                 ! attribution via uniform kernel

              ddx=xl-float(ix)                   ! distance to left cell border
              ddy=yl-float(jy)                   ! distance to lower cell border
              if (ddx.gt.0.5) then
                ixp=ix+1
                wx=1.5-ddx
              else
                ixp=ix-1
                wx=0.5+ddx
              endif

              if (ddy.gt.0.5) then
                jyp=jy+1
                wy=1.5-ddy
              else
                jyp=jy-1
                wy=0.5+ddy
              endif


C Determine mass fractions for four grid points
***********************************************

              if ((ix.ge.0).and.(ix.le.numxgridn-1)) then
                if ((jy.ge.0).and.(jy.le.numygridn-1)) then
                  w=wx*wy
                  do 122 k=1,nspeciesdim
122                 griduncn(ix,jy,kz,nspecpointer(k),nclass(i),nage)=
     +              griduncn(ix,jy,kz,nspecpointer(k),nclass(i),nage)+
     +              xmass1(i,k)/rhoi*weight*w
                endif

                if ((jyp.ge.0).and.(jyp.le.numygridn-1)) then
                  w=wx*(1.-wy)
                  do 125 k=1,nspeciesdim
125                 griduncn(ix,jyp,kz,nspecpointer(k),nclass(i),nage)=
     +              griduncn(ix,jyp,kz,nspecpointer(k),nclass(i),nage)+
     +              xmass1(i,k)/rhoi*weight*w
                endif
              endif


              if ((ixp.ge.0).and.(ixp.le.numxgridn-1)) then
                if ((jyp.ge.0).and.(jyp.le.numygridn-1)) then
                  w=(1.-wx)*(1.-wy)
                  do 123 k=1,nspeciesdim
123                 griduncn(ixp,jyp,kz,nspecpointer(k),nclass(i),nage)=
     +              griduncn(ixp,jyp,kz,nspecpointer(k),nclass(i),nage)+
     +              xmass1(i,k)/rhoi*weight*w
                endif

                if ((jy.ge.0).and.(jy.le.numygridn-1)) then
                  w=(1.-wx)*wy
                  do 124 k=1,nspeciesdim
124                 griduncn(ixp,jy,kz,nspecpointer(k),nclass(i),nage)=
     +              griduncn(ixp,jy,kz,nspecpointer(k),nclass(i),nage)+
     +              xmass1(i,k)/rhoi*weight*w
                endif
              endif
            endif

          endif
        endif
20      continue

************************************************************************
C 2. Evaluate concentrations at receptor points, using the kernel method
************************************************************************

      do 50 n=1,numreceptor


C Reset concentrations
**********************

        do 35 k=1,nspec
35        c(k)=0.


C Estimate concentration at receptor
************************************

        do 40 i=1,numpart

          if (itra1(i).ne.itime) goto 40
          itage=abs(itra1(i)-itramem(i))

          hz=min(50.+0.3*sqrt(float(itage)),hzmax)
          zd=ztra1(i)/hz
          if (zd.gt.1.) goto 40          ! save computing time, leave loop

          hx=min((0.29+2.222e-3*sqrt(float(itage)))*dx+
     +    float(itage)*1.2e-5,hxmax)                     ! 80 km/day
          xd=(xtra1(i)-xreceptor(n))/hx
          if (xd*xd.gt.1.) goto 40       ! save computing time, leave loop

          hy=min((0.18+1.389e-3*sqrt(float(itage)))*dy+
     +    float(itage)*7.5e-6,hymax)                     ! 80 km/day
          yd=(ytra1(i)-yreceptor(n))/hy
          if (yd*yd.gt.1.) goto 40       ! save computing time, leave loop
          h=hx*hy*hz

          r2=xd*xd+yd*yd+zd*zd
          if (r2.lt.1.) then
            xkern=factor*(1.-r2)
            do 45 k=1,nspec
45            c(k)=c(k)+xmass1(i,k)*xkern/h
          endif
40        continue

        do 50 k=1,nspec
50        creceptor(n,k)=creceptor(n,k)+2.*weight*c(k)/receptorarea(n)

      end
