      subroutine readcommand()
********************************************************************************
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subroutine readcommand.       *
*                                                                              *
*     This routine reads the user specifications for the current model run.    *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     18 May 1996                                                              *
*                                                                              *
*     Nov-Dec-2005, R. Easter - input turb_option, add_sfc_level, iouttype     *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* bdate                beginning date as Julian date                           *
* ctl                  factor by which time step must be smaller than          *
*                      Lagrangian time scale                                   *
* edate                ending date as Julian date                              *
* hhh                  hour                                                    *
* ibdate,ibtime        beginnning date and time (YYYYMMDD, HHMISS)             *
* ideltas [s]          modelling period                                        *
* iedate,ietime        ending date and time (YYYYMMDD, HHMISS)                 *
* ifine                reduction factor for vertical wind time step            *
* iflux                switch to turn on (1)/off (0) flux calculations         *
* iout                 1 for conc. (residence time for backward runs) output,  *
*                      2 for mixing ratio output, 3 both, 4 for plume          *
*                      trajectory output, 5 = options 1 and 4                  *
* ipin                 1 continue simulation with dumped particle data, 0 no   *
* ipout                0 no particle dump, 1 every output time, 3 only at end  *
* itsplit [s]          time constant for particle splitting                    *
* loutaver [s]         concentration output is an average over loutaver seconds*
* loutsample [s]       average is computed from samples taken every [s] seconds*
* loutstep [s]         time interval of concentration output                   *
* lsynctime [s]        synchronisation time interval for all particles         *
* lagespectra          switch to turn on (1)/off (0) calculation of age spectra*
* lconvection          value of either 0 and 1 indicating mixing by convection *
*                      = 0 .. no convection                                    *
*                      + 1 .. parameterisation of mixing by subgrid-scale      *
*                              convection = on                                 *
* lsubgrid             switch to turn on (1)/off (0) subgrid topography        *
*                      parameterization                                        *
* method               method used to compute the particle pseudovelocities    *
* mdomainfill          1 use domain-filling option, 0 not, 2 use strat. O3     *
* mi                   minute                                                  *
* ss                   second                                                  *
*                                                                              *
* Constants:                                                                   *
* unitcommand          unit connected to file COMMAND                          *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      integer hhh,mi,ss
      double precision juldate
      character*50 line
      logical old


C Open the command file and read user options
*********************************************


      open(unitcommand,file=path(1)(1:len(1))//'COMMAND',status='old',
     +err=999)

C Check the format of the COMMAND file (either in free format,
C or using formatted mask)
C Use of formatted mask is assumed if line 10 contains the word 'DIRECTION'
***************************************************************************

      call skplin(9,unitcommand)
      read (unitcommand,901) line
901   format (a)
      if (index(line,'LDIRECT') .eq. 0) then
        old = .false.
      else
        old = .true.
      endif
      rewind(unitcommand)


C Read parameters
*****************

      call skplin(7,unitcommand)
      if (old) call skplin(1,unitcommand)

      read(unitcommand,*) ldirect
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) ibdate,ibtime
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) iedate,ietime
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) loutstep
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) loutaver
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) loutsample
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) itsplit
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) lsynctime
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) ctl
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) ifine
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) iout
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) ipout
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) lsubgrid
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) lconvection
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) lagespectra
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) ipin
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) iflux
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) mdomainfill
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) ind_source
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) ind_receptor
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) mquasilag
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) nested_output
c FLEXPART_WRF - read turb_option, add_sfc_level
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) turb_option
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) add_sfc_level
c FLEXPART_WRF - read sfc_option
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) sfc_option
c FLEXPART_WRF - read iouttype
      if (old) call skplin(3,unitcommand)
      read(unitcommand,*) iouttype
      close(unitcommand)

      ifine=max(ifine,1)


C Determine how Markov chain is formulated (for w or for w/sigw)
****************************************************************

      if (ctl.ge.0.1) then
        turbswitch=.true.
      else
        turbswitch=.false.
        ifine=1
      endif
      fine=1./float(ifine)
      ctl=1./ctl

C Set the switches required for the various options for input/output units
**************************************************************************
cAF Set the switches IND_REL and IND_SAMP for the release and sampling
cAf switches for the releasefile:
cAf IND_REL =  1 : xmass * rho
cAf IND_REL =  0 : xmass * 1

cAf switches for the conccalcfile:
cAF IND_SAMP =  0 : xmass * 1
cAf IND_SAMP = -1 : xmass / rho

cAF IND_SOURCE switches between different units for concentrations at the source
cAf   NOTE that in backward simulations the release of computational particles 
cAf   takes place at the "receptor" and the sampling of p[articles at the "source".
cAf          1 = mass units 
cAf          2 = mass mixing ratio units 
cAf IND_RECEPTOR switches between different units for concentrations at the receptor
cAf          1 = mass units 
cAf          2 = mass mixing ratio units 

      if ( ldirect .eq. 1 ) then  ! FWD-Run
cAf set release-switch
         if ( IND_SOURCE .eq. 1 ) then !mass
            ind_rel = 0
         else ! mass mix
            ind_rel = 1
         endif
cAf set sampling switch
         if ( IND_RECEPTOR .eq. 1) then !mass
            ind_samp = 0
         else ! mass mix
            ind_samp = -1
         endif
      elseif (ldirect .eq. -1 ) then !BWD-Run
cAf set sampling switch
         if ( IND_SOURCE .eq. 1 ) then !mass
            ind_samp = -1
         else ! mass mix
            ind_samp = 0
         endif
cAf set release-switch
         if ( IND_RECEPTOR .eq. 1) then !mass
            ind_rel = 1
         else ! mass mix
            ind_rel = 0
         endif
      endif



**************************************************************
C Check whether valid options have been chosen in file COMMAND
**************************************************************

C Check input dates
*******************

      if (iedate.lt.ibdate) then
        write(*,*) ' #### FLEXPART MODEL ERROR! BEGINNING DATE    #### ' 
        write(*,*) ' #### IS LARGER THAN ENDING DATE. CHANGE      #### '
        write(*,*) ' #### EITHER POINT 2 OR POINT 3 IN FILE       #### '
        write(*,*) ' #### "COMMAND".                              #### '
        stop
      else if (iedate.eq.ibdate) then
        if (ietime.lt.ibtime) then
        write(*,*) ' #### FLEXPART MODEL ERROR! BEGINNING TIME    #### ' 
        write(*,*) ' #### IS LARGER THAN ENDING TIME. CHANGE      #### '
        write(*,*) ' #### EITHER POINT 2 OR POINT 3 IN FILE       #### '
        write(*,*) ' #### "COMMAND".                              #### '
        stop
        endif
      endif


C Determine kind of dispersion method
*************************************

      if (ctl.gt.0.) then    
        method=1
        mintime=minstep
      else
        method=0
        mintime=lsynctime
      endif

C Check whether a valid option for gridded model output has been chosen
***********************************************************************

      if ((iout.lt.0).or.(iout.gt.5)) then
        write(*,*) ' #### FLEXPART MODEL ERROR! FILE COMMAND:     #### ' 
        write(*,*) ' #### IOUT MUST BE 0, 1, 2, 3, 4, OR 5!       #### '
        stop
      endif

cAF check consistency between units and volume mixing ratio
      if ( iout.ne.1.and.(ind_source.gt.1 .or. ind_receptor.gt.1) ) then
        write(*,*) ' #### FLEXPART MODEL ERROR! FILE COMMAND:     #### ' 
        write(*,*) ' #### VOLUME MIXING RATIO ONLY SUPPORTED      #### '
        write(*,*) ' #### FOR MASS UNITS (at the moment)          #### '
        stop
      endif
    

C For backward runs, only residence time output (iout=1) or plume trajectories (iout=4),
C or both (iout=5) makes sense; other output options are "forbidden"
****************************************************************************************

      if (ldirect.lt.0) then
        if ((iout.eq.2).or.(iout.eq.3)) then
          write(*,*) '#### FLEXPART MODEL ERROR! FILE COMMAND:     ####' 
          write(*,*) '#### FOR BACKWARD RUNS, IOUT MUST BE 1,4,OR 5####' 
          stop
        endif
      endif


C For domain-filling trajectories, a plume centroid trajectory makes no sense,
C and is "forbidden"
******************************************************************************

      if (mdomainfill.ge.1) then
        if ((iout.eq.4).or.(iout.eq.5)) then
          write(*,*) '#### FLEXPART MODEL ERROR! FILE COMMAND:     ####' 
          write(*,*) '#### FOR DOMAIN-FILLING TRAJECTORY OPTION,   ####' 
          write(*,*) '#### IOUT MUST NOT BE SET TO 4 OR 5.         ####' 
          stop
        endif
      endif
       


C Check whether a valid options for particle dump has been chosen
*****************************************************************

      if ((ipout.ne.0).and.(ipout.ne.1).and.(ipout.ne.2)) then
        write(*,*) ' #### FLEXPART MODEL ERROR! FILE COMMAND:     #### ' 
        write(*,*) ' #### IPOUT MUST BE 1, 2 OR 3!                #### '
        stop
      endif

      if(lsubgrid.ne.1) then
        write(*,*) '             ----------------               '
        write(*,*) ' INFORMATION: SUBGRIDSCALE TERRAIN EFFECT IS'
        write(*,*) ' NOT PARAMETERIZED DURING THIS SIMULATION.  '
        write(*,*) '             ----------------               '
      endif
   

C Check whether convection scheme is either turned on or off
************************************************************

      if ((lconvection.ne.0).and.(lconvection.ne.1)) then
        write(*,*) ' #### FLEXPART MODEL ERROR! FILE COMMAND:     #### ' 
        write(*,*) ' #### LCONVECTION MUST BE SET TO EITHER 1 OR 0#### ' 
        stop
      endif


C Check whether synchronisation interval is sufficiently short
**************************************************************

      if (lsynctime.gt.(idiffnorm/2)) then
        write(*,*) ' #### FLEXPART MODEL ERROR! SYNCHRONISATION   #### ' 
        write(*,*) ' #### TIME IS TOO LONG. MAKE IT SHORTER.      #### '
        stop
      endif


C Check consistency of the intervals, sampling periods, etc., for model output
******************************************************************************

      if (loutaver.eq.0) then
        write(*,*) ' #### FLEXPART MODEL ERROR! TIME AVERAGE OF   #### ' 
        write(*,*) ' #### CONCENTRATION FIELD OUTPUT MUST NOT BE  #### '
        write(*,*) ' #### ZERO.                                   #### '
        write(*,*) ' #### CHANGE INPUT IN FILE COMMAND.           #### '
        stop
      endif

      if (loutaver.gt.loutstep) then
        write(*,*) ' #### FLEXPART MODEL ERROR! TIME AVERAGE OF   #### ' 
        write(*,*) ' #### CONCENTRATION FIELD OUTPUT MUST NOT BE  #### '
        write(*,*) ' #### GREATER THAN INTERVAL OF OUTPUT.        #### '
        write(*,*) ' #### CHANGE INPUT IN FILE COMMAND.           #### '
        stop
      endif

      if (loutsample.gt.loutaver) then
        write(*,*) ' #### FLEXPART MODEL ERROR! SAMPLING TIME OF  #### ' 
        write(*,*) ' #### CONCENTRATION FIELD OUTPUT MUST NOT BE  #### '
        write(*,*) ' #### GREATER THAN TIME AVERAGE OF OUTPUT.    #### '
        write(*,*) ' #### CHANGE INPUT IN FILE COMMAND.           #### '
        stop
      endif

      if (mod(loutaver,lsynctime).ne.0) then
        write(*,*) ' #### FLEXPART MODEL ERROR! AVERAGING TIME OF #### ' 
        write(*,*) ' #### CONCENTRATION FIELD MUST BE A MULTIPLE  #### '
        write(*,*) ' #### OF THE SYNCHRONISATION INTERVAL         #### '
        stop
      endif

      if ((loutaver/lsynctime).lt.2) then
        write(*,*) ' #### FLEXPART MODEL ERROR! AVERAGING TIME OF #### ' 
        write(*,*) ' #### CONCENTRATION FIELD MUST BE AT LEAST    #### '
        write(*,*) ' #### TWICE THE SYNCHRONISATION INTERVAL      #### '
        stop
      endif

      if (mod(loutstep,lsynctime).ne.0) then
        write(*,*) ' #### FLEXPART MODEL ERROR! INTERVAL BETWEEN  #### ' 
        write(*,*) ' #### CONCENTRATION FIELDS MUST BE A MULTIPLE #### '
        write(*,*) ' #### OF THE SYNCHRONISATION INTERVAL         #### '
        stop
      endif

      if ((loutstep/lsynctime).lt.2) then
        write(*,*) ' #### FLEXPART MODEL ERROR! INTERVAL BETWEEN  #### ' 
        write(*,*) ' #### CONCENTRATION FIELDS MUST BE AT LEAST   #### '
        write(*,*) ' #### TWICE THE SYNCHRONISATION INTERVAL      #### '
        stop
      endif

      if (mod(loutsample,lsynctime).ne.0) then
        write(*,*) ' #### FLEXPART MODEL ERROR! SAMPLING TIME OF  #### ' 
        write(*,*) ' #### CONCENTRATION FIELD MUST BE A MULTIPLE  #### '
        write(*,*) ' #### OF THE SYNCHRONISATION INTERVAL         #### '
        stop
      endif

      if (itsplit.lt.loutaver) then
        write(*,*) ' #### FLEXPART MODEL ERROR! SPLITTING TIME FOR#### ' 
        write(*,*) ' #### PARTICLES IS TOO SHORT. PLEASE INCREASE #### '
        write(*,*) ' #### SPLITTING TIME CONSTANT.                #### '
        stop
      endif

      if ((mquasilag.eq.1).and.(iout.ge.4)) then
        write(*,*) ' #### FLEXPART MODEL ERROR! CONFLICTING       #### ' 
        write(*,*) ' #### OPTIONS: IF MQUASILAG=1, PLUME          #### ' 
        write(*,*) ' #### TRAJECTORY OUTPUT IS IMPOSSIBLE.        #### ' 
        stop
      endif

c FLEXPART_WRF - check turb_option, add_sfc_level
      if ((turb_option.ne.turb_option_none     ) .and.
     &    (turb_option.ne.turb_option_diagnosed) .and.
     &    (turb_option.ne.turb_option_tke      )) then
        write(*,*) ' #### FLEXPART MODEL ERROR!                   #### ' 
        write(*,*) ' #### TURB_OPTION MUST BE ONE OF:             #### ' 
        write(*,'(5x,5i5)') turb_option_none, turb_option_diagnosed,
     &             turb_option_tke
        write(*,*) ' #### ---------------------------------       #### ' 
        stop
      endif

      if ((add_sfc_level.ne.0) .and. (add_sfc_level.ne.1)) then
        write(*,*) ' #### FLEXPART MODEL ERROR!                   #### ' 
        write(*,*) ' #### ADD_SFC_LAYER MUST BE 0 or 1            #### ' 
        stop
      endif

      if ((sfc_option.ne.sfc_option_diagnosed) .and.
     &    (sfc_option.ne.sfc_option_wrf      )) then
        write(*,*) ' #### FLEXPART MODEL ERROR!                   #### ' 
        write(*,*) ' #### SFC_OPTION MUST BE ONE OF:              #### ' 
        write(*,'(5x,5i5)') sfc_option_diagnosed, sfc_option_wrf
        write(*,*) ' #### ---------------------------------       #### ' 
        stop
      endif

C iouttype -- convert negative values to 0; positive values to 1
      if (iouttype .lt. 0) iouttype = 0
      if (iouttype .gt. 1) iouttype = 1

C Conversion of format HHHMISS to seconds
*****************************************

      hhh=ideltas/10000
      mi=(ideltas-10000*hhh)/100
      ss=ideltas-10000*hhh-100*mi
      ideltas=hhh*3600+60*mi+ss


C Compute modeling time in seconds and beginning date in Julian date
********************************************************************
 
      outstep=float(abs(loutstep))
      if (ldirect.eq.1) then
        bdate=juldate(ibdate,ibtime)
        edate=juldate(iedate,ietime)
        ideltas=nint((edate-bdate)*86400.)
      else if (ldirect.eq.-1) then
        loutaver=-1*loutaver
        loutstep=-1*loutstep
        loutsample=-1*loutsample
        lsynctime=-1*lsynctime
        bdate=juldate(iedate,ietime)
        edate=juldate(ibdate,ibtime)
        ideltas=nint((edate-bdate)*86400.)
      else
        write(*,*) ' #### FLEXPART MODEL ERROR! DIRECTION IN      #### '
        write(*,*) ' #### FILE "COMMAND" MUST BE EITHER -1 OR 1.  #### '
        stop
      endif

      return    

999   write(*,*) ' #### FLEXPART MODEL ERROR! FILE "COMMAND"    #### ' 
      write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
      write(*,'(a)') path(1)(1:len(1))
      stop

      end
