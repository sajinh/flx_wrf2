      subroutine convmix(itime)
c                          i
c**************************************************************
C     handles all the calculations related to convective mixing
C     Petra Seibert, Bernd C. Krueger, Feb 2001
c     nested grids included, Bernd C. Krueger, May 2001
c
c     Changes by Caroline Forster, April 2004 - February 2005:
c       convmix called every lsynctime seconds
C     CHANGES by A. Stohl:
C       various run-time optimizations - February 2005
C**************************************************************

      include 'includepar'
      include 'includecom'
      include 'includeconv'

      integer igr,igrold, ipart, itime, ix, j, inest
      integer ipconv
      integer jy, kpart, ktop, kz, kzp, ngrid
      integer igrid(maxpart), ipoint(maxpart), igridn(maxpart,maxnests)
C itime [s]                 current time
C igrid(maxpart)            horizontal grid position of each particle
C igridn(maxpart,maxnests)  dto. for nested grids
C ipoint(maxpart)           pointer to access particles according to grid position

      logical lconv
      real x, y, xtn,ytn, ztold, delt
      real dt1,dt2,dtt
      real duma, dumz(nuvzmax+1)
      integer mind1,mind2
C dt1,dt2,dtt,mind1,mind2       variables used for time interpolation
      integer itage,nage

c     monitoring variables
c     real sumconv,sumall


      write(*,'(//a,a//)')
     &    '*** Stopping in subr. convmix ***',
     &    '    This is not implemented for FLEXPART_WRF'
      stop


C Calculate auxiliary variables for time interpolation
******************************************************

      dt1=float(itime-memtime(1))
      dt2=float(memtime(2)-itime)
      dtt=1./(dt1+dt2)
      mind1=memind(1)
      mind2=memind(2)
      delt=float(abs(lsynctime))


      lconv = .false. 

c if no particles are present return after initialization
*********************************************************

      if (numpart.le.0) return

C Assign igrid and igridn, which are pseudo grid numbers indicating particles
C that are outside the part of the grid under consideration
C (e.g. particles near the poles or particles in other nests).
C Do this for all nests but use only the innermost nest; for all others
C igrid shall be -1
C Also, initialize index vector ipoint
*************************************************************************

      do 20 ipart=1,numpart
        igrid(ipart)=-1
        do 21 j=numbnests,1,-1
21        igridn(ipart,j)=-1
        ipoint(ipart)=ipart
C do not consider particles that are (yet) not part of simulation
        if (itra1(ipart).ne.itime) goto 20
        x = xtra1(ipart)
        y = ytra1(ipart)
        
C Determine which nesting level to be used
***********************************************************

        ngrid=0
        do 22 j=numbnests,1,-1
          if ( x.gt.xln(j) .and. x.lt.xrn(j) .and.
     +         y.gt.yln(j) .and. y.lt.yrn(j) ) then
            ngrid=j
            goto 23
          endif
 22       continue
 23     continue
      
C Determine nested grid coordinates
***********************************

        if (ngrid.gt.0) then
C nested grids
          xtn=(x-xln(ngrid))*xresoln(ngrid)
          ytn=(y-yln(ngrid))*yresoln(ngrid)
          ix=nint(xtn)
          jy=nint(ytn)
          igridn(ipart,ngrid) = 1 + jy*nxn(ngrid) + ix
        else if(ngrid.eq.0) then
c mother grid
          ix=nint(x)
          jy=nint(y)
          igrid(ipart) = 1 + jy*nx + ix
        endif

 20     continue

c     sumall = 0. 
c     sumconv = 0.

***************************************************************************************
C 1. Now, do everything for the mother domain and, later, for all of the nested domains
C While all particles have to be considered for redistribution, the Emanuel convection
C scheme only needs to be called once for every grid column where particles are present.
C Therefore, particles are sorted according to their grid position. Whenever a new grid
C cell is encountered by looping through the sorted particles, the convection scheme is called.
***************************************************************************************

C sort particles according to horizontal position and calculate index vector IPOINT

      call sort2(numpart,igrid,ipoint)

C Now visit all grid columns where particles are present
C by going through the sorted particles

      igrold = -1
      do 50 kpart=1,numpart
        igr = igrid(kpart)
        if (igr .eq. -1) goto 50
        ipart = ipoint(kpart)
c       sumall = sumall + 1
        if (igr. ne. igrold) then
C we are in a new grid column
          jy = (igr-1)/nx
          ix = igr - jy*nx - 1

C Interpolate all meteorological data needed for the convection scheme
C Note that tconv & qconv are shifted downward
          psconv=(ps(ix,jy,1,mind1)*dt2+ps(ix,jy,1,mind2)*dt1)*dtt
          tt2conv=(tt2(ix,jy,1,mind1)*dt2+tt2(ix,jy,1,mind2)*dt1)*dtt
          td2conv=(td2(ix,jy,1,mind1)*dt2+td2(ix,jy,1,mind2)*dt1)*dtt
          do 54 kz=1,nconvlev+1
c FLEXPART_WRF - used add_sfc_level for the shifting
c           tconv(kz)=(tth(ix,jy,kz+1,mind1)*dt2+
c    +      tth(ix,jy,kz+1,mind2)*dt1)*dtt
c54         qconv(kz)=(qvh(ix,jy,kz+1,mind1)*dt2+
c    +      qvh(ix,jy,kz+1,mind2)*dt1)*dtt
            kzp = kz + add_sfc_level
            tconv(kz)=(tth(ix,jy,kzp, mind1)*dt2+
     +      tth(ix,jy,kzp, mind2)*dt1)*dtt
54          qconv(kz)=(qvh(ix,jy,kzp, mind1)*dt2+
     +      qvh(ix,jy,kzp, mind2)*dt1)*dtt

C pconv(kz)  = pressure at center of layer k
C phconv(kz) = pressure at bottom boundary of layer k
          do kz = 1, nuvz-add_sfc_level
            kzp = kz + add_sfc_level
            pconv(kz)=(pph(ix,jy,kzp, mind1)*dt2+
     +                 pph(ix,jy,kzp, mind2)*dt1)*dtt
          end do
          do kz = 1, nuvz-add_sfc_level+1
            kzp = kz + add_sfc_level
            dumz(kz) =(zzh(ix,jy,kzp, mind1)*dt2+
     +                 zzh(ix,jy,kzp, mind2)*dt1)*dtt
          end do
          phconv(1) = psconv
          do kz = 2, nuvz-add_sfc_level
            duma = (dumz(kz)-dumz(kz-1))/(dumz(kz+1)-dumz(kz-1))
            phconv(kz)=phconv(kz-1)*(1.0-duma) + pconv(kz)*duma
          end do

C Calculate translocation matrix
          call calcmatrix(lconv,delt,cbaseflux(ix,jy))
          igrold = igr
          ktop = 0
        endif
        
C treat particle only if column has convection
        if (lconv .eqv. .true.) then                
C assign new vertical position to particle

          ztold=ztra1(ipart)
          call redist(ipart,ktop,ipconv)
c         if (ipconv.le.0) sumconv = sumconv+1

C Calculate the gross fluxes across layer interfaces
****************************************************

          if (iflux.eq.1) then
            itage=abs(itra1(ipart)-itramem(ipart))
            do 36 nage=1,nageclass
              if (itage.lt.lage(nage)) goto 37
 36           continue
 37         continue

            if (nage.le.nageclass)
     $      call calcfluxes(nage,ipart,sngl(xtra1(ipart)),
     $      sngl(ytra1(ipart)),ztold)
          endif

        endif   !(lconv .eqv. .true)
 50     continue


************************************************************************************
C 2. Nested domains
************************************************************************************

C sort particles according to horizontal position and calculate index vector IPOINT

      do 70 inest=1,numbnests
        do ipart=1,numpart
          ipoint(ipart)=ipart
          igrid(ipart) = igridn(ipart,inest)
        enddo
        call sort2(numpart,igrid,ipoint)

C Now visit all grid columns where particles are present
C by going through the sorted particles

        igrold = -1
        do 60 kpart=1,numpart
          igr = igrid(kpart)
          if (igr .eq. -1) goto 60
          ipart = ipoint(kpart)
c         sumall = sumall + 1
          if (igr. ne. igrold) then
C we are in a new grid column
            jy = (igr-1)/nxn(inest)
            ix = igr - jy*nxn(inest) - 1

C Interpolate all meteorological data needed for the convection scheme
C Note that tconv & qconv are shifted downward
            psconv=(psn(ix,jy,1,mind1,inest)*dt2+
     +              psn(ix,jy,1,mind2,inest)*dt1)*dtt
            tt2conv=(tt2n(ix,jy,1,mind1,inest)*dt2+
     +               tt2n(ix,jy,1,mind2,inest)*dt1)*dtt
            td2conv=(td2n(ix,jy,1,mind1,inest)*dt2+
     +               td2n(ix,jy,1,mind2,inest)*dt1)*dtt
            do 55 kz=1,nconvlev+1
c FLEXPART_WRF - used add_sfc_level for the shifting
c             tconv(kz)=(tthn(ix,jy,kz+1,mind1,inest)*dt2+
c    +        tthn(ix,jy,kz+1,mind2,inest)*dt1)*dtt
c55           qconv(kz)=(qvhn(ix,jy,kz+1,mind1,inest)*dt2+
c    +        qvhn(ix,jy,kz+1,mind2,inest)*dt1)*dtt
              kzp = kz + add_sfc_level
              tconv(kz)=(tthn(ix,jy,kzp, mind1,inest)*dt2+
     +        tthn(ix,jy,kzp, mind2,inest)*dt1)*dtt
55            qconv(kz)=(qvhn(ix,jy,kzp, mind1,inest)*dt2+
     +        qvhn(ix,jy,kzp, mind2,inest)*dt1)*dtt

C pconv(kz)  = pressure at center of layer k
C phconv(kz) = pressure at bottom boundary of layer k
            do kz = 1, nuvz-add_sfc_level
              kzp = kz + add_sfc_level
              pconv(kz)=(pphn(ix,jy,kzp, mind1,inest)*dt2+
     +                   pphn(ix,jy,kzp, mind2,inest)*dt1)*dtt
            end do
            do kz = 1, nuvz-add_sfc_level+1
              kzp = kz + add_sfc_level
              dumz(kz) =(zzhn(ix,jy,kzp, mind1,inest)*dt2+
     +                   zzhn(ix,jy,kzp, mind2,inest)*dt1)*dtt
            end do
            phconv(1) = psconv
            do kz = 2, nuvz-add_sfc_level
              duma = (dumz(kz)-dumz(kz-1))/(dumz(kz+1)-dumz(kz-1))
              phconv(kz)=phconv(kz-1)*(1.0-duma) + pconv(kz)*duma
            end do

C calculate translocation matrix
********************************
            call calcmatrix(lconv,delt,cbasefluxn(ix,jy,inest))
            igrold = igr
            ktop = 0
          endif
        
C treat particle only if column has convection
          if (lconv .eqv. .true.) then                
C assign new vertical position to particle
            ztold=ztra1(ipart)
            call redist(ipart,ktop,ipconv)
c           if (ipconv.le.0) sumconv = sumconv+1

C Calculate the gross fluxes across layer interfaces
****************************************************

            if (iflux.eq.1) then
              itage=abs(itra1(ipart)-itramem(ipart))
              do 46 nage=1,nageclass
                if (itage.lt.lage(nage)) goto 47
 46             continue
 47           continue

              if (nage.le.nageclass)
     $       call calcfluxes(nage,ipart,sngl(xtra1(ipart)),
     $       sngl(ytra1(ipart)),ztold)
            endif

          endif !(lconv .eqv. .true.)


 60       continue
 70     continue    !inest - loop
c--------------------------------------------------------------------------
c     write(*,*)'############################################'
c     write(*,*)'TIME=',
c    &  itime
c     write(*,*)'fraction of particles under convection',
c    &  sumconv/(sumall+0.001) 
c     write(*,*)'total number of particles',
c    &  sumall 
c     write(*,*)'number of particles under convection',
c    &  sumconv
c     write(*,*)'############################################'

      return
      end
