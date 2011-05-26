      program flexpart_wrf
********************************************************************************
*                                                                              *
*     This is the Lagrangian Particle Dispersion Model FLEXPART_WRF.           *
*                                                                              *
*         FLEXPART uses met. files from the ECMWF model (in grib format),      *
*             and its internal computational grid is latitude-longitude.       *
*                                                                              *
*         FLEXPART_WRF uses met. files from the WRF model (in NetCDF format),  *
*             and its internal computational grid is the WRF x-y grid.         *
*                                                                              *
*     The main program manages the reading of model run specifications, etc.   *
*     All actual computing is done within subroutine timemanager.              *
*                                                                              *
*     Author: A. Stohl                                                         *
*     18 May 1996                                                              *
*                                                                              *
*     Nov 2005, R. Easter - Added the above comments and changed               *
*                           the program name to "flexpart_wrf"                 *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
*                                                                              *
* Constants:                                                                   *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'
      include 'includeconv'

      integer idummy,i,j,ix,jy,inest

      data idummy/-320/


C Generate a large number of random numbers
      write(*,'(/a)') 'FLEXPART_WRF iomode_xycoord setting:'
      if (iomode_xycoord .eq. iomode_xycoord_latlon) then
         write(*,9100) 'degrees longitude,latitude', 
     +                 'degrees longitude,latitude'
      else
         write(*,9100) 'grid-meters', 'grid-meters'
      endif
9100  format( '    x,y coordinates in input  files must be in ', a /
     +        '    x,y coordinates in output files will be in ', a / )


C Generate a large number of random numbers
*******************************************

      do 10 i=1,maxrand-1,2
10      call gasdev1(idummy,rannumb(i),rannumb(i+1))
      call gasdev1(idummy,rannumb(maxrand),rannumb(maxrand-1))


C Read the pathnames where input/output files are stored
********************************************************

      call readpaths()


C Read the user specifications for the current model run
********************************************************

      call readcommand()
       

C Read the age classes to be used
*********************************

      call readageclasses()
       

C Read, which wind fields are available within the modelling period
*******************************************************************

      call readavailable()


C Read the model grid specifications,
C both for the mother domain and eventual nests
***********************************************

      call gridcheck()
      call gridcheck_nests()


C Read the output grid specifications
*************************************

      call readoutgrid()
c     if (nested_output.eq.1) call readoutgrid_nest()
c Nested grid output is not fully implemented in FLEXPART or FLEXPART_WRF
      if (nested_output.ge.1) then
          write(*,'(/a/a/)') 
     &        '*** Nested grid output is not fully implemented  ***',
     &        '*** Set NESTED_OUTPUT=0 in COMMAND file          ***'
          stop
      end if


C Read the receptor points for which extra concentrations are to be calculated
******************************************************************************

      call readreceptors()


C Read the physico-chemical species property table
**************************************************

      call readspecies()


C Read the landuse inventory
****************************

      call readlanduse()


C Assign fractional cover of landuse classes to each ECMWF grid point
*********************************************************************

      call assignland()


C Read the coordinates of the release locations
***********************************************

      call readreleases()


C Read and compute surface resistances to dry deposition of gases
*****************************************************************

      call readdepo()


C Convert the release point coordinates from geografical to grid coordinates
****************************************************************************

      call coordtrafo()


C Initialize all particles to non-existent
******************************************

      do 40 j=1,maxpart
40      itra1(j)=-999999999

C For continuation of previous run, read in particle positions
**************************************************************

      if (ipin.eq.1) then
        call readpartpositions()
      else
        numpart=0
      endif


C Calculate volume, surface area, etc., of all output grid cells
****************************************************************

      call outgrid_init()
      if (nested_output.eq.1) call outgrid_init_nest()


C Write basic information on the simulation to a file "header"
C and open files that are to be kept open throughout the simulation
*******************************************************************

      call writeheader()
      if (nested_output.eq.1) call writeheader_nest()
      open(unitdates,file=path(2)(1:len(2))//'dates')
      call openreceptors()
      if ((iout.eq.4).or.(iout.eq.5)) call openouttraj()

C Releases can only start and end at discrete times (multiples of lsynctime)
****************************************************************************
      
      do 20 i=1,numpoint
        ireleasestart(i)=nint(float(ireleasestart(i))/
     +  float(lsynctime))*lsynctime
20      ireleaseend(i)=nint(float(ireleaseend(i))/
     +  float(lsynctime))*lsynctime


C Initialize cloud-base mass fluxes for the convection scheme
*************************************************************

      do 30 jy=0,nymin1
        do 30 ix=0,nxmin1
30        cbaseflux(ix,jy)=0.
      do 31 inest=1,numbnests
        do 31 jy=0,nyn(inest)-1
          do 31 ix=0,nxn(inest)-1
31          cbasefluxn(ix,jy,inest)=0.


C Calculate particle trajectories
*********************************
      
      call timemanager()


      write(*,'(/a/)') 'CONGRATULATIONS: YOU HAVE SUCCESSFULLY ' // 
     +  'COMPLETED A FLEXPART_WRF MODEL RUN!'

      end
