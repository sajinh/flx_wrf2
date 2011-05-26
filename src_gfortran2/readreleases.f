      subroutine readreleases()
********************************************************************************
*                                                                              *
*     This routine reads the release point specifications for the current      *
*     model run. Several release points can be used at the same time.          *
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subroutine readreleases.      *
*            The computational grid is the WRF x-y grid rather than lat-lon.   *
*                                                                              *
*     Author: A. Stohl                                                         *
*     18 May 1996                                                              *
*                                                                              *
*     Update: 29 January 2001                                                  *
*     Release altitude can be either in magl or masl                           *
*                                                                              *
*     Nov 2005, R. Easter - Do not adjust xpoint1 & 2 by +/-360 degrees        *
*     Dec 2005, R. Easter - x/ypoint1/2 values may be input either as          *
*                           degrees-latlon or grid-meters                      *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* decay               decay constant of species                                *
* dquer [um]          mean particle diameters                                  *
* dsigma              e.g. dsigma=10 or dsigma=0.1 means that 68% of the mass  *
*                     are between 0.1*dquer and 10*dquer                       *
* ireleasestart, ireleaseend [s] starting time and ending time of each release *
* kindz               1: zpoint is in m agl, 2: zpoint is in m asl, 3: zpoint  *
*                     is in hPa                                                *
* npart               number of particles to be released                       *
* nspec               number of species to be released                         *
* density [kg/m3]     density of the particles                                 *
* rm [s/m]            Mesophyll resistance                                     *
* species             name of species                                          *
* xmass               total mass of each species                               *
* xpoint1,ypoint1     geograf. coordinates of lower left corner of release area*
* xpoint2,ypoint2     geograf. coordinates of upper right corner of release are*
* weta, wetb          parameters to determine the wet scavenging coefficient   *
* zpoint1,zpoint2     height range, over which release takes place             *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      integer numpartmax,i,j,id1,it1,id2,it2,idow,ihour
      real vsh(ni),fracth(ni),schmih(ni),releaserate
      double precision jul1,jul2,juldate
      character*50 line
      character*3 aspecnumb
      logical old

      DEP=.false.
      DRYDEP=.false.
      WETDEP=.false.
      do 5 i=1,maxspec
5       DRYDEPSPEC(i)=.false.

C Open the releases file and read user options
**********************************************

      open(unitreleases,file=path(1)(1:len(1))//'RELEASES',status='old',
     +err=999)


C Check the format of the RELEASES file (either in free format,
C or using a formatted mask)
C Use of formatted mask is assumed if line 10 contains the word 'DIRECTION'
***************************************************************************

      call skplin(12,unitreleases)
      read (unitreleases,901) line
901   format (a)
      if (index(line,'Total') .eq. 0) then
        old = .false.
      else
        old = .true.
      endif
      rewind(unitreleases)


C Skip first 11 lines (file header)
***********************************

      call skplin(11,unitreleases)


C Read the number of species and the link to the species information table
C Assign species-specific parameters needed for physical processes
**************************************************************************

      read(unitreleases,*,err=998) nspec
      if (nspec.gt.maxspec) goto 994
      if (old) call skplin(2,unitreleases)
      do 35 i=1,nspec
        read(unitreleases,*,err=998) link(i)
        if (old) call skplin(2,unitreleases)
        species(i)=specname(link(i))

C For backward runs, only 1 species is allowed
**********************************************

      if ((ldirect.lt.0).and.(nspec.gt.1)) then
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
      write(*,*) '#### FOR BACKWARD RUNS, ONLY 1 SPECIES IS ALLOWED####'
      write(*,*) '#####################################################'
        stop
      endif 

C Molecular weight
******************

        weightmolar(i)=weightmol(link(i))
        if (((iout.eq.2).or.(iout.eq.3)).and.
     +  (weightmolar(i).lt.0.)) then
          write(*,*) 'For mixing ratio output, valid molar weight'
          write(*,*) 'must be specified for all simulated species.'
          write(*,*) 'Check table SPECIES or choose concentration'
          write(*,*) 'output instead if molar weight is not known.'
          stop
        endif


C Radioactive decay
*******************

        decay(i)=0.693147/decaytime(link(i)) !conversion half life to decay constant

C Wet deposition
****************

        weta(i)=wetscava(link(i))
        wetb(i)=wetscavb(link(i))

C Dry deposition of gases
*************************

        reldiff(i)=drydiff(link(i))             ! Diffusivity rel. to H20
        henry(i)=dryhenry(link(i))              ! Henry constant
        f0(i)=dryactiv(link(i))                 ! activity
        if (reldiff(i).gt.0.)
     +  rm(i)=1./(henry(i)/3000.+100.*f0(i))    ! mesophyll resistance

C Dry deposition of particles
*****************************

        vsetaver(i)=0.
        density(i)=partrho(link(i))                 ! Particle density
        dquer(i)=partmean(link(i))*1000000.         ! Conversion m to um
        dsigma(i)=partsig(link(i))
        if (density(i).gt.0.) then                  ! Additional parameters
          call part0(dquer(i),dsigma(i),density(i),fracth,schmih,vsh)
          do 80 j=1,ni
            fract(i,j)=fracth(j)
            schmi(i,j)=schmih(j)
            vset(i,j)=vsh(j)
80          vsetaver(i)=vsetaver(i)-vset(i,j)*fract(i,j)
        endif

C Dry deposition for constant deposition velocity
*************************************************

        dryvel(i)=dryvelo(link(i))*0.01         ! conversion to m/s

        if (weta(i).gt.0.) WETDEP=.true.
        if ((reldiff(i).gt.0.).or.(density(i).gt.0.).or.
     +  (dryvel(i).gt.0.)) then
          DRYDEP=.true.
          DRYDEPSPEC(i)=.true.
        endif


C Read in daily and day-of-week variation of emissions, if available
********************************************************************

        do 92 j=1,24           ! initialize everything to no variation
          area_hour(i,j)=1.
92        point_hour(i,j)=1.
        do 93 j=1,7
          area_dow(i,j)=1.
93        point_dow(i,j)=1.

        write(aspecnumb,'(i3.3)') link(i)
        open(unitemissvar,file=path(1)(1:len(1))//'EMISSION_VARIATION_'
     +  //aspecnumb//'.dat',status='old',err=35)
        read(unitemissvar,*)
        do 90 j=1,24     ! 24 hours, starting with 0-1 local time
90        read(unitemissvar,*) ihour,area_hour(i,j),point_hour(i,j)
        read(unitemissvar,*)
        do 91 j=1,7      ! 7 days of the week, starting with Monday
91        read(unitemissvar,*) idow,area_dow(i,j),point_dow(i,j)
        close(unitemissvar)

35      continue

        if (WETDEP.or.DRYDEP) DEP=.true.


C Read specifications for each release point
********************************************

      numpoint=0
      numpartmax=0
      releaserate=0.
100   numpoint=numpoint+1
      read(unitreleases,*,end=25)

      read(unitreleases,*,err=998,end=25) id1,it1
      if (numpoint.gt.maxpoint) goto 997
      if (old) call skplin(2,unitreleases)
      read(unitreleases,*,err=998) id2,it2
      if (old) call skplin(2,unitreleases)
      read(unitreleases,*,err=998) xpoint1(numpoint)
      if (old) call skplin(2,unitreleases)
      read(unitreleases,*,err=998) ypoint1(numpoint)
      if (old) call skplin(2,unitreleases)
      read(unitreleases,*,err=998) xpoint2(numpoint)
      if (old) call skplin(2,unitreleases)
      read(unitreleases,*,err=998) ypoint2(numpoint)
      if (old) call skplin(2,unitreleases)
      read(unitreleases,*,err=998) kindz(numpoint)
      if (old) call skplin(2,unitreleases)
      read(unitreleases,*,err=998) zpoint1(numpoint)
      if (old) call skplin(2,unitreleases)
      read(unitreleases,*,err=998) zpoint2(numpoint)
      if (old) call skplin(2,unitreleases)
      read(unitreleases,*,err=998) npart(numpoint)
      if (old) call skplin(2,unitreleases)
      do 20 i=1,nspec
        read(unitreleases,*,err=998) xmass(numpoint,i)
20      if (old) call skplin(2,unitreleases)
      compoint(numpoint) = ' '
      read(unitreleases,'(a40)',err=998) compoint(numpoint)(1:40)
      if (old) call skplin(1,unitreleases)
      if((xpoint1(numpoint).eq.0.).and.(ypoint1(numpoint).eq.0.).and.
     +(xpoint2(numpoint).eq.0.).and.(ypoint2(numpoint).eq.0.).and.
     +(compoint(numpoint)(1:8).eq.'        ')) goto 25


      j = numpoint
      write(*,'(/a,i7)') 'readreleases diagnostics - numpoint = ', j
      write(*,'(a,1p,2e18.10)') 'x, ypoint1 (in) ', 
     &   xpoint1(j), ypoint1(j)
      write(*,'(a,1p,2e18.10)') 'x, ypoint2 (in) ', 
     &   xpoint2(j), ypoint2(j)
      if (iomode_xycoord .eq. iomode_xycoord_latlon) then
C In this case, the above inputs are the actual geographical lat/lon
C   of the southwest & northeast corners of the release area
C Need to convert from lat/lon to grid-meters
         releases_swlon(j) = xpoint1(j)
         releases_swlat(j) = ypoint1(j)
         releases_nelon(j) = xpoint2(j)
         releases_nelat(j) = ypoint2(j)
         call ll_to_xymeter_wrf( releases_swlon(j), releases_swlat(j), 
     &      xpoint1(j), ypoint1(j) )
         call ll_to_xymeter_wrf( releases_nelon(j), releases_nelat(j), 
     &      xpoint2(j), ypoint2(j) )
         write(*,'(a,1p,2e18.10)') 'x, ypoint1      ', 
     &      xpoint1(j), ypoint1(j)
         write(*,'(a,1p,2e18.10)') 'x, ypoint2      ', 
     &      xpoint2(j), ypoint2(j)
      else
C In this case, the above inputs are in grid-meters 
C Need to convert from grid-meters to lat/lon
         call xymeter_to_ll_wrf( xpoint1(j), ypoint1(j),
     &      releases_swlon(j), releases_swlat(j) )
         call xymeter_to_ll_wrf( xpoint2(j), ypoint2(j),
     &      releases_nelon(j), releases_nelat(j) )
         write(*,'(f15.10,5x,a)') releases_swlon(j), 'releases_swlon'
         write(*,'(f15.10,5x,a)') releases_swlat(j), 'releases_swlat'
         write(*,'(f15.10,5x,a)') releases_nelon(j), 'releases_nelon'
         write(*,'(f15.10,5x,a)') releases_nelat(j), 'releases_nelat'
      end if

C If a release point contains no particles, stop and issue error message
************************************************************************

      if (npart(numpoint).eq.0) then
        write(*,*) 'FLEXPART MODEL ERROR'
        write(*,*) 'RELEASES file is corrupt.'
        write(*,*) 'At least for one release point, there are zero'
        write(*,*) 'particles released. Make changes to RELEASES.'
        stop
      endif

C Check whether x coordinates of release point are within model domain
**********************************************************************

c FLEXPART_WRF - x & y coords are in meters, so the following lines 
c   (which adjust longitude by +/-360 degrees) are not needed
c
c      if (xpoint1(numpoint).lt.xlon0) 
c    +       xpoint1(numpoint)=xpoint1(numpoint)+360.
c      if (xpoint1(numpoint).gt.xlon0+(nxmin1)*dx)
c    +       xpoint1(numpoint)=xpoint1(numpoint)-360.
c      if (xpoint2(numpoint).lt.xlon0) 
c    +       xpoint2(numpoint)=xpoint2(numpoint)+360.
c      if (xpoint2(numpoint).gt.xlon0+(nxmin1)*dx)
c    +       xpoint2(numpoint)=xpoint2(numpoint)-360.

C Determine relative beginning and ending times of particle release
*******************************************************************

      jul1=juldate(id1,it1)
      jul2=juldate(id2,it2)
      if (jul1.gt.jul2) then
        write(*,*) 'FLEXPART MODEL ERROR'
        write(*,*) 'Release stops before it begins.'
        write(*,*) 'Make changes to file RELEASES.'
        stop
      endif
      if (mdomainfill.eq.0) then   ! no domain filling
        if (ldirect.eq.1) then
          if ((jul1.lt.bdate).or.(jul2.gt.edate)) then
            write(*,*) 'FLEXPART MODEL ERROR'
            write(*,*) 'Release starts before simulation begins or ends'
            write(*,*) 'after simulation stops.'
            write(*,*) 'Make files COMMAND and RELEASES consistent.'
            stop
          endif
          ireleasestart(numpoint)=int((jul1-bdate)*86400.)
          ireleaseend(numpoint)=int((jul2-bdate)*86400.)
        else if (ldirect.eq.-1) then
          if ((jul1.lt.edate).or.(jul2.gt.bdate)) then
            write(*,*) 'FLEXPART MODEL ERROR'
            write(*,*) 'Release starts before simulation begins or ends'
            write(*,*) 'after simulation stops.'
            write(*,*) 'Make files COMMAND and RELEASES consistent.'
            stop
          endif
          ireleasestart(numpoint)=int((jul1-bdate)*86400.)
          ireleaseend(numpoint)=int((jul2-bdate)*86400.)
        endif
      endif


C Check, whether the total number of particles may exceed totally allowed
C number of particles at some time during the simulation
*************************************************************************

C Determine the release rate (particles per second) and total number
C of particles released during the simulation
********************************************************************

      if (ireleasestart(numpoint).ne.ireleaseend(numpoint)) then
        releaserate=releaserate+float(npart(numpoint))/
     +  float(ireleaseend(numpoint)-ireleasestart(numpoint))
      else
        releaserate=99999999
      endif
      numpartmax=numpartmax+npart(numpoint)
      goto 100

25    close(unitreleases)
      numpoint=numpoint-1

      if (releaserate.gt.
     +0.99*float(maxpart)/float(lage(nageclass))) then
        if (numpartmax.gt.maxpart) then
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
      write(*,*) '####                                             ####'
      write(*,*) '####WARNING - TOTAL NUMBER OF PARTICLES SPECIFIED####'
      write(*,*) '#### IN FILE "RELEASES" MAY AT SOME POINT DURING ####'
      write(*,*) '#### THE SIMULATION EXCEED THE MAXIMUM ALLOWED   ####'
      write(*,*) '#### NUMBER (MAXPART).IF RELEASES DO NOT OVERLAP,####'
      write(*,*) '#### FLEXPART CAN POSSIBLY COMPLETE SUCCESSFULLY.####'
      write(*,*) '#### HOWEVER, FLEXPART MAY HAVE TO STOP          ####'
      write(*,*) '#### AT SOME TIME DURING THE SIMULATION. PLEASE  ####'
      write(*,*) '#### MAKE SURE THAT YOUR SETTINGS ARE CORRECT.   ####'
      write(*,*) '#####################################################'
          write(*,*) 'Maximum release rate may be: ',releaserate,
     +    ' particles per second'
          write(*,*) 'Maximum allowed release rate is: ',
     +    float(maxpart)/float(lage(nageclass)),' particles per second'
          write(*,*)
     +'Total number of particles released during the simulation is: ',
     +    numpartmax
          write(*,*) 'Maximum allowed number of particles is: ',maxpart
        endif
      endif


C Make a consistency check, whether the forward/backward switch is correctly set
********************************************************************************

      if (ldirect.eq.1) then
        if (maxpointspec.lt.nspec) then
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - PARAMETER MAXPOINTSPEC IS NOT       ####'
      write(*,*) '#### CORRECTLY SET FOR A FORWARD SIMULATION.     ####'
      write(*,*) '#### CHANGE APPROPRIATELY IN FILE INCLUDEPAR.    ####'
      write(*,*) '#####################################################'
        endif
      else
        if (maxpointspec.lt.numpoint) then
      write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - PARAMETER MAXPOINTSPEC IS NOT       ####'
      write(*,*) '#### CORRECTLY SET FOR A BACKWARD SIMULATION.    ####'
      write(*,*) '#### CHANGE APPROPRIATELY IN FILE INCLUDEPAR.    ####'
      write(*,*) '#####################################################'
        endif
      endif

      return


994   write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - MAXIMUM NUMBER OF EMITTED SPECIES IS####'
      write(*,*) '#### TOO LARGE. PLEASE REDUCE NUMBER OF SPECIES. ####'
      write(*,*) '#####################################################'
      stop


997   write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### ERROR - NUMBER OF RELEASE POINTS SPECIFIED  ####'
      write(*,*) '#### IN FILE "RELEASES" EXCEEDS THE MAXIMUM      ####'
      write(*,*) '#### ALLOWED NUMBER.                             ####'
      write(*,*) '#####################################################'
      stop


998   write(*,*) '#####################################################'
      write(*,*) '#### FLEXPART MODEL SUBROUTINE READRELEASES:     ####'
      write(*,*) '####                                             ####'
      write(*,*) '#### FATAL ERROR - FILE "RELEASES" IS            ####'
      write(*,*) '#### CORRUPT. PLEASE CHECK YOUR INPUTS FOR       ####'
      write(*,*) '#### MISTAKES OR GET A NEW "RELEASES"-           ####'
      write(*,*) '#### FILE ...                                    ####'
      write(*,*) '#####################################################'
      stop


999   write(*,*) '#####################################################'
      write(*,*) '   FLEXPART MODEL SUBROUTINE READRELEASES: '
      write(*,*)
      write(*,*) 'FATAL ERROR - FILE CONTAINING PARTICLE RELEASE POINTS'
      write(*,*) 'POINTS IS NOT AVAILABLE OR YOU ARE NOT'
      write(*,*) 'PERMITTED FOR ANY ACCESS'
      write(*,*) '#####################################################'
      stop

      end
