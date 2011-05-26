      subroutine timemanager()                                  

********************************************************************************
*                                                                              *
* Handles the computation of trajectories, i.e. determines which               *
* trajectories have to be computed at what time.                               *
* Manages dry+wet deposition routines, radioactive decay and the computation   *
* of concentrations.                                                           *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     20 May 1996                                                              *
*                                                                              *
*     Dec 2005, J. Fast - Only call conccalc & concoutput when (iout.ge.1)     *
*                                                                              *
********************************************************************************
*  Changes, Bernd C. Krueger, Feb. 2001:                                       *
*        Call of convmix when new windfield is read                            *
*------------------------------------                                          *
*  Changes Petra Seibert, Sept 2002                                            *
*     fix wet scavenging problem                                               *
*     Code may not be correct for decay of deposition!                         *
*  Changes Petra Seibert, Nov 2002                                             *
*     call convection BEFORE new fields are read in BWD mode                   *
*  Changes Caroline Forster, Feb 2005
*     new interface between flexpart and convection scheme
*     Emanuel's latest subroutine convect43c.f is used
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* DEP                .true. if either wet or dry deposition is switched on     *
* decay(maxspec) [1/s] decay constant for radioactive decay                    *
* DRYDEP             .true. if dry deposition is switched on                   *
* ideltas [s]        modelling period                                          *
* itime [s]          actual temporal position of calculation                   *
* ldeltat [s]        time since computation of radioact. decay of depositions  *
* loutaver [s]       averaging period for concentration calculations           *
* loutend [s]        end of averaging for concentration calculations           *
* loutnext [s]       next time at which output fields shall be centered        *
* loutsample [s]     sampling interval for averaging of concentrations         *
* loutstart [s]      start of averaging for concentration calculations         *
* loutstep [s]       time interval for which concentrations shall be calculated*
* npoint(maxpart)    index, which starting point the trajectory has            *
*                    starting positions of trajectories                        *
* nstop              serves as indicator for fate of particles                 *
*                    in the particle loop                                      *
* nstop1             serves as indicator for wind fields (see getfields)       *
* outnum             number of samples for each concentration calculation      *
* outnum             number of samples for each concentration calculation      *
* prob               probability of absorption at ground due to dry deposition *
* WETDEP             .true. if wet deposition is switched on                   *
* weight             weight for each concentration sample (1/2 or 1)           *
* uap(maxpart),ucp(maxpart),uzp(maxpart) = random velocities due to turbulence *
* us(maxpart),vs(maxpart),ws(maxpart) = random velocities due to interpolation *
* xtra1(maxpart), ytra1(maxpart), ztra1(maxpart) =                             *
*                    spatial positions of trajectories                         *
*                                                                              *
* Constants:                                                                   *
* maxpart            maximum number of trajectories                            *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      integer j,k,l,n,itime,nstop,nstop1
      integer loutnext,loutstart,loutend
      integer ix,jy,ldeltat,itage,nage
      real outnum,weight,prob(maxspec)
      real uap(maxpart),ucp(maxpart),uzp(maxpart),decfact
      real us(maxpart),vs(maxpart),ws(maxpart),cbt(maxpart)
      real drydeposit(maxspec),gridtotalunc,wetgridtotalunc
      real drygridtotalunc,xold,yold,zold
c     real xm,xm1


C First output for time 0
*************************

      loutnext=loutstep/2
      outnum=0.
      loutstart=loutnext-loutaver/2
      loutend=loutnext+loutaver/2


***********************************************************************
C Loop over the whole modelling period in time steps of mintime seconds
***********************************************************************

      do 10 itime=0,ideltas,lsynctime

C Computation of wet deposition every lsynctime seconds
C maybe wet depo frequency can be relaxed later but better be on safe side
C wetdepo must be called BEFORE new fields are read in but should not
C be called in the very beginning before any fields are loaded, or
C before particles are in the system
C Code may not be correct for decay of deposition
C changed by Petra Seibert 9/02
*********************************************************************

        if (WETDEP .and. itime .ne. 0 .and. numpart .gt. 0)
     &    call wetdepo(itime,lsynctime,loutnext)

c compute convection for backward runs
**************************************

          if ((ldirect.eq.-1).and.(lconvection.eq.1).and.(itime.lt.0))
     &    call convmix(itime)

C Get necessary wind fields if not available
********************************************

        call getfields(itime,nstop1)
        if (nstop1.gt.1) stop 'NO METEO FIELDS AVAILABLE'

C Release particles
*******************

        if (mdomainfill.ge.1) then
          if (itime.eq.0) then
            call init_domainfill()
          else
            call boundcond_domainfill(itime,loutend)
          endif
        else
          call releaseparticles(itime)
        endif


C Compute convective mixing for forward runs
c for backward runs it is done before next windfield is read in
***************************************************************

          if ((ldirect.eq.1).and.(lconvection.eq.1))
     &    call convmix(itime)


C If middle of averaging period of output fields is reached, accumulated
C deposited mass radioactively decays 
************************************************************************

        if (DEP.and.(itime.eq.loutnext)) then 
          do 41 k=1,nspec
            if (decay(k).gt.0.) then
              do 42 nage=1,nageclass
                do 42 l=1,nclassunc
C Mother output grid
                  do 43 jy=0,numygrid-1
                    do 43 ix=0,numxgrid-1
                      wetgridunc(ix,jy,k,l,nage)=
     +                wetgridunc(ix,jy,k,l,nage)*
     +                exp(-1.*outstep*decay(k))
43                    drygridunc(ix,jy,k,l,nage)=
     +                drygridunc(ix,jy,k,l,nage)*
     +                exp(-1.*outstep*decay(k))
C Nested output grid
                  if (nested_output.eq.1) then
                    do 44 jy=0,numygridn-1
                      do 44 ix=0,numxgridn-1
                        wetgriduncn(ix,jy,k,l,nage)=
     +                  wetgriduncn(ix,jy,k,l,nage)*
     +                  exp(-1.*outstep*decay(k))
44                      drygriduncn(ix,jy,k,l,nage)=
     +                  drygriduncn(ix,jy,k,l,nage)*
     +                  exp(-1.*outstep*decay(k))
                  endif
42              continue
            endif
41          continue
        endif

!!! CHANGE: These lines may be switched on to check the conservation
!!! of mass within FLEXPART

c       if (mod(itime,loutsample).eq.0) then 
c          xm=0.
c          xm1=0.
c          do 247 j=1,numpart
c47          if (itra1(j).eq.itime) xm1=xm1+xmass1(j,1)
c          xm=xm1
c          do 248 nage=1,nageclass
c            do 248 ix=0,numxgrid-1
c              do 248 jy=0,numygrid-1
c                do 248 l=1,nclassunc
c48        xm=xm+wetgridunc(ix,jy,1,l,nage)+drygridunc(ix,jy,1,l,nage)
c          write(*,'(i6,4f10.3)') itime,xm,xm1
c       endif
!!! CHANGE

          
C Check whether concentrations are to be calculated
***************************************************

        if ((ldirect*itime.ge.ldirect*loutstart).and.
     +  (ldirect*itime.le.ldirect*loutend)) then ! add to grid
          if (mod(itime-loutstart,loutsample).eq.0) then

C If we are exactly at the start or end of the concentration averaging interval,
C give only half the weight to this sample
********************************************************************************

            if ((itime.eq.loutstart).or.(itime.eq.loutend)) then
              weight=0.5
            else
              weight=1.0
            endif
            outnum=outnum+weight
            if(iout.ge.1) call conccalc(itime,weight)
          endif


          if ((mquasilag.eq.1).and.(itime.eq.(loutstart+loutend)/2))
     +    call partoutput_short(itime)    ! dump particle positions in extremely compressed format


C Output and reinitialization of grid
C If necessary, first sample of new grid is also taken
******************************************************

          if ((itime.eq.loutend).and.(outnum.gt.0.)) then
            if ((iout.le.3.).or.(iout.eq.5)) then 
              if(iout.ge.1) call concoutput(itime,outnum,gridtotalunc,
     +        wetgridtotalunc,drygridtotalunc)
              if (nested_output.eq.1.and.iout.ge.1)
     +           call concoutput_nest(itime,outnum)
              outnum=0.
            endif
            if ((iout.eq.4).or.(iout.eq.5)) call plumetraj(itime)
            if (iflux.eq.1) call fluxoutput(itime)
            write(*,45) itime,numpart,gridtotalunc,wetgridtotalunc,
     +      drygridtotalunc
45          format(i9,' SECONDS SIMULATED: ',i8,
     +      ' PARTICLES:    Uncertainty: ',3f7.3)
            if (ipout.ge.1) call partoutput(itime)    ! dump particle positions
            loutnext=loutnext+loutstep
            loutstart=loutnext-loutaver/2
            loutend=loutnext+loutaver/2
            if (itime.eq.loutstart) then
              weight=0.5
              outnum=outnum+weight
              if(iout.ge.1) call conccalc(itime,weight)
            endif


C Check, whether particles are to be split:
C If so, create new particles and attribute all information from the old
C particles also to the new ones; old and new particles both get half the
C mass of the old ones
*************************************************************************

            if (ldirect*itime.ge.ldirect*itsplit) then
              n=numpart
              do 30 j=1,numpart
                if (ldirect*itime.ge.ldirect*itrasplit(j)) then
                  if (n.lt.maxpart) then
                    n=n+1
                    itrasplit(j)=2*(itrasplit(j)-itramem(j))+itramem(j)
                    itrasplit(n)=itrasplit(j)
                    itramem(n)=itramem(j)
                    itra1(n)=itra1(j)
                    idt(n)=idt(j)
                    npoint(n)=npoint(j)
                    nclass(n)=nclass(j)
                    xtra1(n)=xtra1(j)
                    ytra1(n)=ytra1(j)
                    ztra1(n)=ztra1(j)
                    uap(n)=uap(j)
                    ucp(n)=ucp(j)
                    uzp(n)=uzp(j)
                    us(n)=us(j)
                    vs(n)=vs(j)
                    ws(n)=ws(j)
                    cbt(n)=cbt(j)
                    do 35 k=1,nspec
                      xmass1(j,k)=xmass1(j,k)/2.
35                    xmass1(n,k)=xmass1(j,k)
                  endif
                endif
30              continue
              numpart=n
            endif
          endif
        endif
        

        if (itime.eq.ideltas) goto 99

C Compute interval since radioactive decay of deposited mass was computed
*************************************************************************

        if (itime.lt.loutnext) then
          ldeltat=itime-(loutnext-loutstep)
        else                                  ! first half of next interval
          ldeltat=itime-loutnext
        endif


C Loop over all particles
*************************

        do 20 j=1,numpart


C If integration step is due, do it
***********************************

          if (itra1(j).eq.itime) then

C Determine age class of the particle
            itage=abs(itra1(j)-itramem(j))
            do 36 nage=1,nageclass
              if (itage.lt.lage(nage)) goto 37
36            continue
37          continue

C Initialize newly released particle
************************************

            if ((itramem(j).eq.itime).or.(itime.eq.0))
     +      call initialize(itime,idt(j),uap(j),ucp(j),uzp(j),
     +      us(j),vs(j),ws(j),xtra1(j),ytra1(j),ztra1(j),cbt(j))

C Memorize particle positions
*****************************

            xold=xtra1(j)
            yold=ytra1(j)
            zold=ztra1(j)

C Integrate Lagevin equation for lsynctime seconds
**************************************************

            call advance(itime,idt(j),uap(j),ucp(j),uzp(j),us(j),
     +      vs(j),ws(j),nstop,xtra1(j),ytra1(j),ztra1(j),prob,cbt(j))


C Calculate the gross fluxes across layer interfaces
****************************************************

            if (iflux.eq.1) call calcfluxes(nage,j,xold,yold,zold)


C Determine, when next time step is due
C If trajectory is terminated, mark it
***************************************

            if (nstop.gt.1) then
              itra1(j)=-999999999
            else
              itra1(j)=itime+lsynctime


C Dry deposition and radioactive decay for each species
*******************************************************

              do 23 k=1,nspec
                if (decay(k).gt.0.) then             ! radioactive decay
                  decfact=exp(-float(abs(lsynctime))*decay(k))
                else
                  decfact=1.
                endif

                if (DRYDEPSPEC(k)) then        ! dry deposition
                  drydeposit(k)=xmass1(j,k)*prob(k)*decfact
                  xmass1(j,k)=xmass1(j,k)*(1.-prob(k))*decfact
                  if (decay(k).gt.0.) then   ! correct for decay (see wetdepo)
                    drydeposit(k)=drydeposit(k)*
     +              exp(float(abs(ldeltat))*decay(k))
                  endif
                else                           ! no dry deposition
                  xmass1(j,k)=xmass1(j,k)*decfact
                endif
23              continue

              if (DRYDEP) then
                call drydepokernel(nclass(j),drydeposit,sngl(xtra1(j)),
     +          sngl(ytra1(j)),nage)
                if (nested_output.eq.1) call drydepokernel_nest(
     +          nclass(j),drydeposit,sngl(xtra1(j)),sngl(ytra1(j)),nage)
              endif

C Terminate trajectories that are older than maximum allowed age
****************************************************************

              if (abs(itra1(j)-itramem(j)).ge.lage(nageclass))
     +        itra1(j)=-999999999
            endif

          endif

20        continue                         ! end of loop over particles
          

10      continue                           ! end of loop over simulation time

99    continue

      if (ipout.eq.2) call partoutput(itime)     ! dump particle positions

      end
