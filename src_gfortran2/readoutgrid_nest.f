      subroutine readoutgrid_nest()
********************************************************************************
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subroutine readoutgrid_nest.  *
*                                                                              *
*     This routine reads the user specifications for the output nest.          *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     4 June 1996                                                              *
*                                                                              *
*     Dec 2005, R. Easter -                                                    *
*         The nested output grid is defined by specifying its southwest and    *
*         northeast corners, either in degrees-latlon or grid-meters           *
*         Changes to some read formats (wider fields).                         *
*         Changed names of "*lon0*" & "*lat0*" variables                       *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* dxoutn,dyoutn        grid distances of output nest                           *
* numxgridn,numygridn,numzgrid    nest dimensions                              *
* out_xm0n,out_ym0n    lower left corner of nest                               *
* outheight(maxzgrid)  height levels of output grid [m]                        *
*                                                                              *
* Constants:                                                                   *
* unitoutgrid          unit connected to file OUTGRID                          *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      real xr,xr1,yr,yr1,eps
      real xtmp, xtmp1, xtmp2, ytmp, ytmp1, ytmp2
      parameter(eps=1.e-4)



C Open the OUTGRID file and read output grid specifications
***********************************************************

      open(unitoutgrid,file=path(1)(1:len(1))//'OUTGRID_NEST',
     +status='old',err=999)


      call skplin(5,unitoutgrid)


C 1.  Read horizontal grid specifications
*****************************************

      call skplin(3,unitoutgrid)
      read(unitoutgrid,'(f15.8)') xtmp1
      call skplin(3,unitoutgrid)
      read(unitoutgrid,'(f15.8)') ytmp1
      call skplin(3,unitoutgrid)
      read(unitoutgrid,'(4x,i5)') numxgridn
      call skplin(3,unitoutgrid)
      read(unitoutgrid,'(4x,i5)') numygridn
      call skplin(3,unitoutgrid)
      read(unitoutgrid,'(f15.8)') xtmp2
      call skplin(3,unitoutgrid)
      read(unitoutgrid,'(f15.8)') ytmp2

      write(*,'(/a)') 'readoutgrid_nest diagnostics'
      write(*,'(a,1p,2e18.10)') 'xtmp1, ytmp1  (in)', xtmp1, ytmp1
      write(*,'(a,1p,2e18.10)') 'xtmp2, ytmp2  (in)', xtmp2, ytmp2
      if (iomode_xycoord .eq. iomode_xycoord_latlon) then
C In this case, the above inputs are the actual geographical lat/lon
C   of the southwest & northeast corners of the output grid
C Need to convert from lat/lon to grid-meters
         outgridn_swlon = xtmp1
         outgridn_swlat = ytmp1
         outgridn_nelon = xtmp2
         outgridn_nelat = ytmp2
         call ll_to_xymeter_wrf( outgridn_swlon, outgridn_swlat, 
     &      out_xm0n, out_ym0n )
         call ll_to_xymeter_wrf( outgridn_nelon, outgridn_nelat, 
     &      xtmp, ytmp )
         dxoutn = (xtmp - out_xm0n)/numxgridn
         dyoutn = (ytmp - out_ym0n)/numygridn
      else
C In this case, the above inputs are in grid-meters 
C Need to convert from grid-meters to lat/lon
         out_xm0n = xtmp1
         out_ym0n = ytmp1
         dxoutn = (xtmp2 - xtmp1)/numxgridn
         dyoutn = (ytmp2 - ytmp1)/numygridn
         call xymeter_to_ll_wrf( xtmp1, ytmp1, 
     &      outgridn_swlon, outgridn_swlat )
         call xymeter_to_ll_wrf( xtmp2, ytmp2, 
     &      outgridn_nelon, outgridn_nelat )
         write(*,'(f15.10,5x,a)') outgridn_swlon, 'outgridn_swlon'
         write(*,'(f15.10,5x,a)') outgridn_swlat, 'outgridn_swlat'
         write(*,'(f15.10,5x,a)') outgridn_nelon, 'outgridn_nelon'
         write(*,'(f15.10,5x,a)') outgridn_nelat, 'outgridn_nelat'
      end if
      write(*,'(a,1p,2e18.10)') 'out_xm0n, out_ym0n', out_xm0n, out_ym0n
      write(*,'(a,1p,2e18.10)') 'dxoutn,   dyoutn  ', dxoutn, dyoutn


C Check validity of output grid (shall be within model domain)
**************************************************************

      xr=out_xm0n+float(numxgridn)*dxoutn
      yr=out_ym0n+float(numygridn)*dyoutn
      xr1=xmet0+float(nxmin1)*dx
      yr1=ymet0+float(nymin1)*dy
      if ((out_xm0n+eps.lt.xmet0).or.(out_ym0n+eps.lt.ymet0)
     +.or.(xr.gt.xr1+eps).or.(yr.gt.yr1+eps)) then
        write(*,*) ' #### FLEXPART MODEL ERROR! PART OF OUTPUT    ####'
        write(*,*) ' #### NEST IS OUTSIDE MODEL DOMAIN. CHANGE    ####'
        write(*,*) ' #### FILE OUTGRID IN DIRECTORY               ####'
        write(*,'(a)') path(1)(1:len(1))
        stop
      endif
      if ((numxgridn.gt.maxxgridn).or.(numygridn.gt.maxygridn)) then
        write(*,*) ' #### FLEXPART MODEL ERROR! DIMENSIONS OF     ####'
        write(*,*) ' #### OUTPUT NEST EXCEED MAXIMUM VALUES.      ####'
        write(*,*) ' #### CHANGE FILE $FLEXPART/options/OUTGRID.  ####'
        stop
      endif

      xoutshiftn=xmet0-out_xm0n
      youtshiftn=ymet0-out_ym0n
      close(unitoutgrid)
      return


999   write(*,*) ' #### FLEXPART MODEL ERROR! FILE OUTGRID_NEST #### '
      write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
      write(*,*) ' #### xxx/flexpart/options                    #### '
      stop

      end
