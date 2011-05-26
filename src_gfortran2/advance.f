      subroutine advance(itime,ldt,up,vp,wp,
     +usigold,vsigold,wsigold,nstop,xt,yt,zt,prob,cbt)
C                          i    i  i/oi/oi/o
C       i/o     i/o     i/o     o  i/oi/oi/o i/o  i/o
********************************************************************************
*                                                                              *
*  Note:  This is the FLEXPART_WRF version of subroutine gridcheck.            *
*    The computational grid is the WRF x-y grid rather than lat-lon.           *
*                                                                              *
*  Calculation of turbulent particle trajectories utilizing a                  *
*  zero-acceleration scheme, which is corrected by a numerically more accurate *
*  Petterssen scheme whenever possible.                                        *
*                                                                              *
*  Particle positions are read in, incremented, and returned to the calling    *
*  program.                                                                    *
*                                                                              *
*  In different regions of the atmosphere (PBL vs. free troposphere),          *
*  different parameters are needed for advection, parameterizing turbulent     *
*  velocities, etc. For efficiency, different interpolation routines have      *
*  been written for these different cases, with the disadvantage that there    *
*  exist several routines doing almost the same. They all share the            *
*  included file 'includeinterpol'. The following                              *
*  interpolation routines are used:                                            *
*                                                                              *
*  interpol_all(_nests)     interpolates everything (called inside the PBL)    *
*  interpol_misslev(_nests) if a particle moves vertically in the PBL,         *
*                           additional parameters are interpolated if it       *
*                           crosses a model level                              *
*  interpol_wind(_nests)    interpolates the wind and determines the           *
*                           standard deviation of the wind (called outside PBL)*
*                           also interpolates potential vorticity              *
*  interpol_wind_short(_nests) only interpolates the wind (needed for the      *
*                           Petterssen scheme)                                 *
*  interpol_vdep(_nests)    interpolates deposition velocities                 *
*                                                                              *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     16 December 1997                                                         *
*                                                                              *
*  Changes:                                                                    *
*                                                                              *
*  8 April 2000: Deep convection parameterization                              *
*                                                                              *
*  May 2002: Petterssen scheme introduced                                      *
*                                                                              *
*  26 Oct 2005, R. Easter - changes for horizontal grid in m instead of lat,lon*
*  10 Nov 2005, R. Easter - zero turbulent wind components is                  *
*                           turbulence is turned off                           *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* cbt                1 if particle not transferred to forbidden state, else -1 *
* dawsave            accumulated displacement in along-wind direction          *
* dcwsave            accumulated displacement in cross-wind direction          *
* dxsave             accumulated displacement in longitude                     *
* dysave             accumulated displacement in latitude                      *
* h [m]              Mixing height                                             *
* lwindinterv [s]    time interval between two wind fields                     *
* itime [s]          time at which this subroutine is entered                  *
* itimec [s]         actual time, which is incremented in this subroutine      *
* href [m]           height for which dry deposition velocity is calculated    *
* ladvance [s]       Total integration time period                             *
* ldirect            1 forward, -1 backward                                    *
* ldt [s]            Time step for the next integration                        *
* lsynctime [s]      Synchronisation interval of FLEXPART                      *
* ngrid              index which grid is to be used                            *
* nrand              index for a variable to be picked from rannumb            *
* nstop              if > 1 particle has left domain and must be stopped       *
* prob               probability of absorption due to dry deposition           *
* rannumb(maxrand)   normally distributed random variables                     *
* rhoa               air density                                               *
* rhograd            vertical gradient of the air density                      *
* up,vp,wp           random velocities due to turbulence (along wind, cross    *
*                    wind, vertical wind                                       *
* usig,vsig,wsig     mesoscale wind fluctuations                               *
* usigold,vsigold,wsigold  like usig, etc., but for the last time step         *
* vdepo              Deposition velocities for all species                     *
* xt,yt,zt           Particle position                                         *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'
      include 'includeinterpol'
      include 'includehanna'

      double precision xt,yt
      real zt,xts,yts
      integer itime,itimec,nstop,idummy,ldt,i,j,k,nrand,loop,memindnext
      integer ngr
      real dz,dz1,dz2,cbt,xlon,ylat,xpol,ypol,gridsize,cgszll
      real ru,rv,rw,dt,ux,vy,cosfact,xtn,ytn
      real prob(maxspec),up,vp,wp,dxsave,dysave,dawsave,dcwsave
      real usigold,vsigold,wsigold,r,rs
      real uold,vold,wold,vdepo(maxspec)
c     real uprof(nzmax),vprof(nzmax),wprof(nzmax)
c     real usigprof(nzmax),vsigprof(nzmax),wsigprof(nzmax)
c     real rhoprof(nzmax),rhogradprof(nzmax)
      real rhoa,rhograd,ran3,delz,dtf,rhoaux,dtftlw,eps,uxscale,wpscale
      parameter(eps=nxmax/3.e5)
      save idummy


!!! CHANGE: TEST OF THE WELL-MIXED CRITERION
c       integer iclass
c       parameter(iclass=10)
c       double precision zacc,tacc,t(iclass),th(0:iclass),hsave
c       logical dump
c       save zacc,tacc,t,th,hsave,dump
!!! CHANGE

      data idummy/-7/


!!! CHANGE: TEST OF THE WELL-MIXED CRITERION
c     if (idummy.eq.-7) then
c     open(550,file='WELLMIXEDTEST')
c     do 17 i=0,iclass
c7      th(i)=float(i)/float(iclass)
c     endif
!!! CHANGE


      nstop=0
      do 70 i=1,nmixz
70      indzindicator(i)=.true.


      if (DRYDEP) then    ! reset probability for deposition
        do 11 i=1,nspec
          depoindicator(i)=.true.
11        prob(i)=0.
      endif

      dxsave=0.           ! reset position displacements
      dysave=0.           ! due to mean wind
      dawsave=0.          ! and turbulent wind
      dcwsave=0.

      itimec=itime

      nrand=int(ran3(idummy)*float(maxrand-1))+1


C Determine whether lat/long grid or polarstereographic projection
C is to be used
C Furthermore, determine which nesting level to be used
******************************************************************

      if (nglobal.and.(yt.gt.switchnorthg)) then
        write(*,*)
        write(*,*) '*** stopping in advance ***'
        write(*,*) '    the n-pole code section should not be active'
        write(*,*)
        ngrid=-1
      else if (sglobal.and.(yt.lt.switchsouthg)) then
        write(*,*)
        write(*,*) '*** stopping in advance ***'
        write(*,*) '    the s-pole code section should not be active'
        write(*,*)
        ngrid=-2
      else
        ngrid=0
        do 22 j=numbnests,1,-1
          if ((xt.gt.xln(j)).and.(xt.lt.xrn(j)).and.
     +    (yt.gt.yln(j)).and.(yt.lt.yrn(j))) then
            ngrid=j
            goto 23
          endif
22        continue
23      continue
      endif


****************************
C Interpolate necessary data
****************************

      if (abs(itime-memtime(1)).lt.abs(itime-memtime(2))) then
        memindnext=1
      else
        memindnext=2
      endif

C Determine nested grid coordinates
***********************************

      if (ngrid.gt.0) then
        xtn=(xt-xln(ngrid))*xresoln(ngrid)
        ytn=(yt-yln(ngrid))*yresoln(ngrid)
        ix=int(xtn)
        jy=int(ytn)
      else
        ix=int(xt)
        jy=int(yt)
      endif
      ixp=ix+1
      jyp=jy+1


C Compute maximum mixing height around particle position
********************************************************

      h=0.
      if (ngrid.le.0) then
        do 27 k=1,2
          do 27 j=jy,jyp
            do 27 i=ix,ixp
              if (hmix(i,j,1,k).gt.h) h=hmix(i,j,1,k)
27            continue
      else
        do 28 k=1,2
          do 28 j=jy,jyp
            do 28 i=ix,ixp
              if (hmixn(i,j,1,k,ngrid).gt.h) h=hmixn(i,j,1,k,ngrid)
28            continue
      endif

      zeta=zt/h
 


**************************************************************
C If particle is in the PBL, interpolate once and then make a
C time loop until end of interval is reached
**************************************************************

      if (zeta.le.1.) then

C BEGIN TIME LOOP
C================

        loop=0
100       loop=loop+1
          if (method.eq.1) then
            ldt=min(ldt,abs(lsynctime-itimec+itime))
            itimec=itimec+ldt*ldirect
          else
            ldt=abs(lsynctime)
            itimec=itime+lsynctime
          endif
          dt=float(ldt)

          zeta=zt/h


          if (loop.eq.1) then
            if (ngrid.le.0) then
              xts=sngl(xt)
              yts=sngl(yt)
              call interpol_all(itime,xts,yts,zt)
            else
              call interpol_all_nests(itime,xtn,ytn,zt)
            endif

          else


C Determine the level below the current position for u,v,rho
************************************************************

            do 5 i=2,nz
              if (height(i).gt.zt) then
                indz=i-1
                indzp=i
                goto 6
              endif
5             continue
6           continue

C If one of the levels necessary is not yet available,
C calculate it
******************************************************

            do 12 i=indz,indzp
              if (indzindicator(i)) then
                if (ngrid.le.0) then
                  call interpol_misslev(i)
                else
                  call interpol_misslev_nests(i)
                endif
              endif
12            continue
          endif


C Vertical interpolation of u,v,w,rho and drhodz
************************************************

C Vertical distance to the level below and above current position
C both in terms of (u,v) and (w) fields
*****************************************************************

          dz=1./(height(indzp)-height(indz))
          dz1=(zt-height(indz))*dz
          dz2=(height(indzp)-zt)*dz

          u=dz1*uprof(indzp)+dz2*uprof(indz)
          v=dz1*vprof(indzp)+dz2*vprof(indz)
          w=dz1*wprof(indzp)+dz2*wprof(indz)
          rhoa=dz1*rhoprof(indzp)+dz2*rhoprof(indz)
          rhograd=dz1*rhogradprof(indzp)+dz2*rhogradprof(indz)


C Compute the turbulent disturbances
C Determine the sigmas and the timescales
*****************************************

          if (turbswitch) then
            call hanna(zt)
          else
            call hanna1(zt)
          endif
            

******************************************
C Determine the new diffusivity velocities
******************************************

C Horizontal components
***********************

          if (nrand+1.gt.maxrand) nrand=1
          if (dt/tlu.lt..5) then
            up=(1.-dt/tlu)*up+rannumb(nrand)*sigu*sqrt(2.*dt/tlu)
          else
            ru=exp(-dt/tlu)
            up=ru*up+rannumb(nrand)*sigu*sqrt(1.-ru**2)
          endif
          if (dt/tlv.lt..5) then
            vp=(1.-dt/tlv)*vp+rannumb(nrand+1)*sigv*sqrt(2.*dt/tlv)
          else
            rv=exp(-dt/tlv)
            vp=rv*vp+rannumb(nrand+1)*sigv*sqrt(1.-rv**2)
          endif
          nrand=nrand+2


          if (nrand+ifine.gt.maxrand) nrand=1
          rhoaux=rhograd/rhoa
          dtf=dt*fine

          dtftlw=dtf/tlw

C Loop over ifine short time steps for vertical component
*********************************************************

          do 20 i=1,ifine

C Determine the drift velocity and density correction velocity
**************************************************************

            if (turbswitch) then
              if (dtftlw.lt..5) then
                wp=((1.-dtftlw)*wp+rannumb(nrand+i)*sqrt(2.*dtftlw)
     +          +dtf*(dsigwdz+rhoaux*sigw))*cbt
              else
                rw=exp(-dtftlw)
                wp=(rw*wp+rannumb(nrand+i)*sqrt(1.-rw**2)
     +          +tlw*(1.-rw)*(dsigwdz+rhoaux*sigw))*cbt
              endif
              delz=wp*sigw*dtf
            else
              rw=exp(-dtftlw)
              wp=(rw*wp+rannumb(nrand+i)*sqrt(1.-rw**2)*sigw
     +        +tlw*(1.-rw)*(dsigw2dz+rhoaux*sigw**2))*cbt
              delz=wp*dtf
            endif

c FLEXPART_WRF - zero up,vp,wp if turbulence is turned off
            if (turb_option .eq. turb_option_none) then
              up=0.0
              vp=0.0
              wp=0.0
            end if

*****************************************************
C Compute turbulent vertical displacement of particle
*****************************************************

            if (abs(delz).gt.h) delz=mod(delz,h)

C Determine if particle transfers to a "forbidden state" below the ground
C or above the mixing height
*************************************************************************

            if (delz.lt.-zt) then         ! reflection at ground
              cbt=-1.
              zt=-zt-delz
            else if (delz.gt.(h-zt)) then ! reflection at h
              cbt=-1.
              zt=-zt-delz+2.*h
            else                         ! no reflection
              cbt=1.
              zt=zt+delz
            endif

            if (i.ne.ifine) then
              zeta=zt/h
              call hanna_short(zt)
            endif

20          continue
          nrand=nrand+i

C Determine time step for next integration
******************************************

          if (turbswitch) then
            ldt=int(min(tlw,h/max(2.*abs(wp*sigw),1.e-5),
     +      0.5/abs(dsigwdz))*ctl)
          else
            ldt=int(min(tlw,h/max(2.*abs(wp),1.e-5))*ctl)
          endif
          ldt=max(ldt,mintime)


C If particle represents only a single species, add gravitational settling
C velocity. The settling velocity is zero for gases, or if particle
C represents more than one species
**************************************************************************

          if (nspec.eq.1) w=w+vsetaver(1)

C Horizontal displacements during time step dt are small real values compared
C to the position; adding the two, would result in large numerical errors.
C Thus, displacements are accumulated during lsynctime and are added to the
C position at the end
*****************************************************************************

          dxsave=dxsave+u*dt
          dysave=dysave+v*dt
          dawsave=dawsave+up*dt
          dcwsave=dcwsave+vp*dt
          zt=zt+w*dt*float(ldirect)

          if (zt.gt.h) then
            if (itimec.eq.itime+lsynctime) goto 99
            goto 700    ! complete the current interval above PBL
          endif
          if (zt.lt.0.) zt=-1.*zt    ! if particle below ground -> refletion


!!! CHANGE: TEST OF THE WELL-MIXED CRITERION
!!! These lines may be switched on to test the well-mixed criterion
c     if (zt.le.h) then
c       zacc=zacc+zt/h*dt
c       hsave=hsave+h*dt
c       tacc=tacc+dt
c       do 67 i=1,iclass
c         if ((zt/h.gt.th(i-1)).and.(zt/h.le.th(i)))
c    +    t(i)=t(i)+dt
c7        continue
c     endif
c     if ((mod(itime,10800).eq.0).and.dump) then
c      dump=.false.
c      write(550,'(i6,12f10.3)') itime,hsave/tacc,zacc/tacc,
c    + (t(i)/tacc*float(iclass),i=1,iclass)
c       zacc=0.
c       tacc=0.
c       do 68 i=1,iclass
c8        t(i)=0.
c       hsave=0.
c     endif
c     if (mod(itime,10800).ne.0) dump=.true.
!!! CHANGE


C Determine probability of deposition
*************************************

          if ((DRYDEP).and.(zt.lt.2.*href)) then
            do 10 i=1,nspec
              if (DRYDEPSPEC(i)) then
                if (depoindicator(i)) then
                  if (ngrid.le.0) then
                    call interpol_vdep(i,vdepo(i))
                  else
                    call interpol_vdep_nests(i,vdepo(i))
                  endif
                endif
c correction by Petra Seibert, 10 April 2001
c   this formulation means that prob(n) = 1 - f(0)*...*f(n)
c   where f(n) is the exponential term
                prob(i)=1.+(prob(i)-1.)*exp(-vdepo(i)*abs(dt)/(2.*href))
              endif
10          continue
          endif

          if (itimec.eq.(itime+lsynctime)) then
            usig=0.5*(usigprof(indzp)+usigprof(indz))
            vsig=0.5*(vsigprof(indzp)+vsigprof(indz))
            wsig=0.5*(wsigprof(indzp)+wsigprof(indz))
            goto 99  ! finished
          endif
          goto 100

C END TIME LOOP
C==============


      endif



***********************************************************
C For all particles that are outside the PBL, make a single
C time step. Only horizontal turbulent disturbances are
C calculated. Vertical disturbances are reset.
***********************************************************


C Interpolate the wind
**********************

700   continue
      if (ngrid.le.0) then
        xts=sngl(xt)
        yts=sngl(yt)
        call interpol_wind(itime,xts,yts,zt)
      else
        call interpol_wind_nests(itime,xtn,ytn,zt)
      endif


C Compute everything for above the PBL

C Assume constant, uncorrelated, turbulent perturbations
C In the stratosphere, use a small vertical diffusivity d_strat,
C in the troposphere, use a larger horizontal diffusivity d_trop.
C Turbulent velocity scales are determined based on sqrt(d_trop/dt)
*******************************************************************

      ldt=abs(lsynctime-itimec+itime)
      dt=float(ldt)
      
      if ((zt.lt.5000.).or.((zt.lt.18000.).and.(abs(pvi).lt.2.))) then  ! in the troposphere
        uxscale=sqrt(d_trop/dt)
        if (nrand+1.gt.maxrand) nrand=1
        ux=rannumb(nrand)*uxscale
        vy=rannumb(nrand+1)*uxscale
        nrand=nrand+2
        wp=0.
      else                 ! in the stratosphere
        ux=0.
        vy=0.
        wpscale=sqrt(d_strat/dt)
        wp=rannumb(nrand)*wpscale
      endif

c FLEXPART_WRF - zero ux,vy,wp if turbulence is turned off
      if (turb_option .eq. turb_option_none) then
        ux=0.0
        vy=0.0
        wp=0.0
      end if


C If particle represents only a single species, add gravitational settling
C velocity. The settling velocity is zero for gases
**************************************************************************

      if (nspec.eq.1) w=w+vsetaver(1)

C Calculate position at time step itime+lsynctime
*************************************************

      dxsave=dxsave+(u+ux)*dt
      dysave=dysave+(v+vy)*dt
      zt=zt+(w+wp)*dt*float(ldirect)
      if (zt.lt.0.) zt=-1.*zt    ! if particle below ground -> refletion

99    continue



*****************************************************************
C Add mesoscale random disturbances
C This is done only once for the whole lsynctime interval to save
C computation time
*****************************************************************


C Mesoscale wind velocity fluctuations are obtained by scaling
C with the standard deviation of the grid-scale winds surrounding
C the particle location, multiplied by a factor turbmesoscale.
C The autocorrelation time constant is taken as half the
C time interval between wind fields
*****************************************************************

      r=exp(-2.*float(abs(lsynctime))/float(lwindinterv))
      rs=sqrt(1.-r**2)
      if (nrand+2.gt.maxrand) nrand=1
      usigold=r*usigold+rs*rannumb(nrand)*usig*turbmesoscale
      vsigold=r*vsigold+rs*rannumb(nrand+1)*vsig*turbmesoscale
      wsigold=r*wsigold+rs*rannumb(nrand+2)*wsig*turbmesoscale

c FLEXPART_WRF - zero u/v/wsigold if turbulence is turned off
c Note: for mesoscale model applications this component should be ignored!
c     if (turb_option .eq. turb_option_none) then
        usigold=0.0
        vsigold=0.0
        wsigold=0.0
c     end if

      dxsave=dxsave+usigold*float(lsynctime)
      dysave=dysave+vsigold*float(lsynctime)

      zt=zt+wsigold*float(lsynctime)
      if (zt.lt.0.) zt=-1.*zt    ! if particle below ground -> refletion

**************************************************************
C Transform along and cross wind components to xy coordinates,
C add them to u and v, transform u,v to grid units/second
C and calculate new position
**************************************************************

      call windalign(dxsave,dysave,dawsave,dcwsave,ux,vy)
      dxsave=dxsave+ux
      dysave=dysave+vy
      if (ngrid.ge.0) then
c for FLEXPART_WRF, dx & dy are in meters,
c dxconst=1/dx, dyconst=1/dy, and no cos(lat) is needed
c       cosfact=dxconst/cos((yt*dy+ylat0)*pi180)
c       xt=xt+dble(dxsave*cosfact*float(ldirect))
c       yt=yt+dble(dysave*dyconst*float(ldirect))
        xt=xt+dble(dxsave*dxconst*float(ldirect))
        yt=yt+dble(dysave*dyconst*float(ldirect))
c     else if (ngrid.eq.-1) then      ! around north pole
c       xlon=xlon0+xt*dx
c       ylat=ylat0+yt*dy
c       call cll2xy(northpolemap,ylat,xlon,xpol,ypol)
c       gridsize=1000.*cgszll(northpolemap,ylat,xlon)
c       dxsave=dxsave/gridsize
c       dysave=dysave/gridsize
c       xpol=xpol+dxsave*float(ldirect)
c       ypol=ypol+dysave*float(ldirect)
c       call cxy2ll(northpolemap,xpol,ypol,ylat,xlon)
c       xt=(xlon-xlon0)/dx
c       yt=(ylat-ylat0)/dy
c     else if (ngrid.eq.-2) then    ! around south pole
c       xlon=xlon0+xt*dx
c       ylat=ylat0+yt*dy
c       call cll2xy(southpolemap,ylat,xlon,xpol,ypol)
c       gridsize=1000.*cgszll(southpolemap,ylat,xlon)
c       dxsave=dxsave/gridsize
c       dysave=dysave/gridsize
c       xpol=xpol+dxsave*float(ldirect)
c       ypol=ypol+dysave*float(ldirect)
c       call cxy2ll(southpolemap,xpol,ypol,ylat,xlon)
c       xt=(xlon-xlon0)/dx
c       yt=(ylat-ylat0)/dy
      else
        write(*,*) 'advance -- bad ngrid = ', ngrid
        stop
      endif


C If global data are available, use cyclic boundary condition
*************************************************************

      if (xglobal) then
        if (xt.ge.float(nxmin1)) xt=xt-float(nxmin1)
        if (xt.lt.0.) xt=xt+float(nxmin1)
        if (xt.le.eps) xt=eps
        if (abs(xt-float(nxmin1)).le.eps) xt=float(nxmin1)-eps
      endif


C Check position: If trajectory outside model domain, terminate it
******************************************************************

      if ((xt.lt.0.).or.(xt.ge.float(nxmin1)).or.(yt.lt.0.).or.
     +(yt.ge.float(nymin1))) then
        nstop=3
        return
      endif

C If particle above highest model level, set it back into the domain
********************************************************************

      if (zt.ge.height(nz)) zt=height(nz)-100.*eps


*************************************************************************
C Now we could finish, as this was done in FLEXPART versions up to 4.0.
C However, truncation errors of the advection can be significantly
C reduced by doing one iteration of the Petterssen scheme, if this is
C possible.
C Note that this is applied only to the grid-scale winds, not to
C the turbulent winds.
*************************************************************************

C The Petterssen scheme can only applied with long time steps (only then u
C is the "old" wind as required by the scheme); otherwise do nothing
**************************************************************************

      if (ldt.ne.abs(lsynctime)) return

C The Petterssen scheme can only be applied if the ending time of the time step
C (itime+ldt*ldirect) is still between the two wind fields held in memory;
C otherwise do nothing
*******************************************************************************

      if (abs(itime+ldt*ldirect).gt.abs(memtime(2))) return  

C Apply it also only if starting and ending point of current time step are on
C the same grid; otherwise do nothing
******************************************************************************
      if (nglobal.and.(yt.gt.switchnorthg)) then
        ngr=-1
      else if (sglobal.and.(yt.lt.switchsouthg)) then
        ngr=-2
      else
        ngr=0
        do 42 j=numbnests,1,-1
          if ((xt.gt.xln(j)).and.(xt.lt.xrn(j)).and.
     +    (yt.gt.yln(j)).and.(yt.lt.yrn(j))) then
            ngr=j
            goto 43
          endif
42        continue
43      continue
      endif

      if (ngr.ne.ngrid) return

C Determine nested grid coordinates
***********************************

      if (ngrid.gt.0) then
        xtn=(xt-xln(ngrid))*xresoln(ngrid)
        ytn=(yt-yln(ngrid))*yresoln(ngrid)
        ix=int(xtn)
        jy=int(ytn)
      else
        ix=int(xt)
        jy=int(yt)
      endif 
      ixp=ix+1
      jyp=jy+1


C Memorize the old wind
***********************

      uold=u
      vold=v
      wold=w

C Interpolate wind at new position and time
*******************************************

      if (ngrid.le.0) then
        xts=sngl(xt)
        yts=sngl(yt)
        call interpol_wind_short(itime+ldt*ldirect,xts,yts,zt)
      else
        call interpol_wind_short_nests(itime+ldt*ldirect,xtn,ytn,zt)
      endif

C Determine the difference vector between new and old wind
C (use half of it to correct position according to Petterssen)
**************************************************************

      u=(u-uold)/2.
      v=(v-vold)/2.
      w=(w-wold)/2.


C Finally, correct the old position
***********************************

      zt=zt+w*float(ldt*ldirect)
      if (zt.lt.0.) zt=-1.*zt    ! if particle below ground -> refletion
      if (ngrid.ge.0) then
c for FLEXPART_WRF, dx & dy are in meters,
c dxconst=1/dx, dyconst=1/dy, and no cos(lat) is needed
c       cosfact=dxconst/cos((yt*dy+ylat0)*pi180)
c       xt=xt+dble(u*cosfact*float(ldt*ldirect))
c       yt=yt+dble(v*dyconst*float(ldt*ldirect))
        xt=xt+dble(u*dxconst*float(ldt*ldirect))
        yt=yt+dble(v*dyconst*float(ldt*ldirect))
c     else if (ngrid.eq.-1) then      ! around north pole
c       xlon=xlon0+xt*dx
c       ylat=ylat0+yt*dy
c       call cll2xy(northpolemap,ylat,xlon,xpol,ypol)
c       gridsize=1000.*cgszll(northpolemap,ylat,xlon)
c       u=u/gridsize
c       v=v/gridsize
c       xpol=xpol+u*float(ldt*ldirect)
c       ypol=ypol+v*float(ldt*ldirect)
c       call cxy2ll(northpolemap,xpol,ypol,ylat,xlon)
c       xt=(xlon-xlon0)/dx
c       yt=(ylat-ylat0)/dy
c     else if (ngrid.eq.-2) then    ! around south pole
c       xlon=xlon0+xt*dx
c       ylat=ylat0+yt*dy
c       call cll2xy(southpolemap,ylat,xlon,xpol,ypol)
c       gridsize=1000.*cgszll(southpolemap,ylat,xlon)
c       u=u/gridsize
c       v=v/gridsize
c       xpol=xpol+u*float(ldt*ldirect)
c       ypol=ypol+v*float(ldt*ldirect)
c       call cxy2ll(southpolemap,xpol,ypol,ylat,xlon)
c       xt=(xlon-xlon0)/dx
c       yt=(ylat-ylat0)/dy
      else
        write(*,*) 'advance -- bad ngrid = ', ngrid
        stop
      endif

C If global data are available, use cyclic boundary condition
*************************************************************

      if (xglobal) then
        if (xt.ge.float(nxmin1)) xt=xt-float(nxmin1)
        if (xt.lt.0.) xt=xt+float(nxmin1)
        if (xt.le.eps) xt=eps
        if (abs(xt-float(nxmin1)).le.eps) xt=float(nxmin1)-eps
      endif

C Check position: If trajectory outside model domain, terminate it
******************************************************************

      if ((xt.lt.0.).or.(xt.ge.float(nxmin1)).or.(yt.lt.0.).or.
     +(yt.ge.float(nymin1))) then
        nstop=3
        return
      endif

C If particle above highest model level, set it back into the domain
********************************************************************

      if (zt.ge.height(nz)) zt=height(nz)-100.*eps


      end
