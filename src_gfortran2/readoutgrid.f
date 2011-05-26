      subroutine readoutgrid()
********************************************************************************
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subroutine readoutgrid.       *
*                                                                              *
*     This routine reads the user specifications for the output grid.          *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     4 June 1996                                                              *
*                                                                              *
*     Dec 2005, R. Easter -                                                    *
*             The output grid is defined by specifying its southwest and       *
*             northease corners, either in degrees-latlon or grid-meters       *
*             Changes to some read formats (wider fields).                     *
*             Changed names of "*lon0*" & "*lat0*" variables                   *
*     10 Mar 2006, R. Easter -                                                 *
*             Change eps from 1.0e-4 (degrees value) to 10.0 (meters value)    *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* dxout,dyout          grid distance                                           *
* numxgrid,numygrid,numzgrid    grid dimensions                                *
* out_xm0,out_ym0      lower left corner of grid                               *
* outheight(maxzgrid)  height levels of output grid [m]                        *
*                                                                              *
* Constants:                                                                   *
* unitoutgrid          unit connected to file OUTGRID                          *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      integer i,j
      real outhelp,xr,xr1,yr,yr1,eps
      real xtmp, xtmp1, xtmp2, ytmp, ytmp1, ytmp2
C 10-mar-2006 rce - flexpart_wrf - eps should be a small dx/dy value in meters
C     parameter(eps=1.e-4)
      parameter(eps=10.0)



C Open the OUTGRID file and read output grid specifications
***********************************************************

      open(unitoutgrid,file=path(1)(1:len(1))//'OUTGRID',status='old',
     +err=999)


      call skplin(5,unitoutgrid)


C 1.  Read horizontal grid specifications
C
C *** NOTE ***
C [xtmp1, ytmp1] are the coordinates of the southwest corner
C    of the first (i.e., southwest or lower-left) output grid cell
C [xtmp2, ytmp2] are the coordinates of the northeast corner 
C    of the last (i.e,, northeast or upper-right) output grid cell
*****************************************

      call skplin(3,unitoutgrid)
      read(unitoutgrid,'(f15.8)') xtmp1
      call skplin(3,unitoutgrid)
      read(unitoutgrid,'(f15.8)') ytmp1
      call skplin(3,unitoutgrid)
      read(unitoutgrid,'(4x,i5)') numxgrid
      call skplin(3,unitoutgrid)
      read(unitoutgrid,'(4x,i5)') numygrid
      call skplin(3,unitoutgrid)
      read(unitoutgrid,'(f15.8)') xtmp2
      call skplin(3,unitoutgrid)
      read(unitoutgrid,'(f15.8)') ytmp2

      write(*,'(/a)') 'readoutgrid diagnostics'
      write(*,'(a,1p,2e18.10)') 'xtmp1, ytmp1 (in)', xtmp1, ytmp1
      write(*,'(a,1p,2e18.10)') 'xtmp2, ytmp2 (in)', xtmp2, ytmp2
      if (iomode_xycoord .eq. iomode_xycoord_latlon) then
C In this case, the above inputs are the actual geographical lat/lon
C   of the southwest & northeast corners of the output grid
C Need to convert from lat/lon to grid-meters
         outgrid_swlon = xtmp1
         outgrid_swlat = ytmp1
         outgrid_nelon = xtmp2
         outgrid_nelat = ytmp2
         call ll_to_xymeter_wrf( outgrid_swlon, outgrid_swlat, 
     &      out_xm0, out_ym0 )
         call ll_to_xymeter_wrf( outgrid_nelon, outgrid_nelat, 
     &      xtmp, ytmp )
C 10-mar-2006 rce
C If out_xm0 is very close to mother grid sw corner (=xmet0), set it to that
C If xtmp    is very close to mother grid ne corner (=xmet0+nxmin1*dx), set it to that
C Do similar for out_ym0 & ytmp
         if (abs(out_xm0-xmet0) .le. eps) out_xm0 = xmet0
         if (abs(out_ym0-ymet0) .le. eps) out_ym0 = ymet0
         xr1 = xmet0 + float(nxmin1)*dx
         yr1 = ymet0 + float(nymin1)*dy
         if (abs(xtmp-xr1) .le. eps) xtmp = xr1
         if (abs(ytmp-yr1) .le. eps) ytmp = yr1
         dxout = (xtmp - out_xm0)/numxgrid
         dyout = (ytmp - out_ym0)/numygrid
      else
C In this case, the above inputs are in grid-meters 
C Need to convert from grid-meters to lat/lon
         out_xm0 = xtmp1
         out_ym0 = ytmp1
         dxout = (xtmp2 - xtmp1)/numxgrid
         dyout = (ytmp2 - ytmp1)/numygrid
         call xymeter_to_ll_wrf( xtmp1, ytmp1, 
     &      outgrid_swlon, outgrid_swlat )
         call xymeter_to_ll_wrf( xtmp2, ytmp2, 
     &      outgrid_nelon, outgrid_nelat )
         write(*,'(f15.10,5x,a)') outgrid_swlon, 'outgrid_swlon'
         write(*,'(f15.10,5x,a)') outgrid_swlat, 'outgrid_swlat'
         write(*,'(f15.10,5x,a)') outgrid_nelon, 'outgrid_nelon'
         write(*,'(f15.10,5x,a)') outgrid_nelat, 'outgrid_nelat'
      end if
      write(*,'(a,1p,2e18.10)') 'out_xm0, out_ym0 ', out_xm0, out_ym0
      write(*,'(a,1p,2e18.10)') 'dxout,   dyout   ', dxout, dyout


C Check validity of output grid (shall be within model domain)
**************************************************************

      xr=out_xm0+float(numxgrid)*dxout
      yr=out_ym0+float(numygrid)*dyout
      xr1=xmet0+float(nxmin1)*dx
      yr1=ymet0+float(nymin1)*dy
      if ((out_xm0+eps.lt.xmet0).or.(out_ym0+eps.lt.ymet0)
     +.or.(xr.gt.xr1+eps).or.(yr.gt.yr1+eps)) then
        write(*,*) ' #### FLEXPART MODEL ERROR! PART OF OUTPUT    ####'
        write(*,*) ' #### GRID IS OUTSIDE MODEL DOMAIN. CHANGE    ####'
        write(*,*) ' #### FILE OUTGRID IN DIRECTORY               ####'
        write(*,'(a)') path(1)(1:len(1))
        stop
      endif
      if ((numxgrid.gt.maxxgrid).or.(numygrid.gt.maxygrid)) then
        write(*,*) ' #### FLEXPART MODEL ERROR! DIMENSIONS OF     ####'
        write(*,*) ' #### OUTPUT GRID EXCEED MAXIMUM VALUES.      ####'
        write(*,*) ' #### CHANGE FILE $FLEXPART/options/OUTGRID.  ####'
        stop
      endif

C 2. Vertical levels of output grid
***********************************

      j=0
100     j=j+1
        do 20 i=1,3
20        read(unitoutgrid,*,end=99)
        read(unitoutgrid,'(4x,f7.1)',end=99) outhelp
        if (outhelp.eq.0.) goto 99
        if (j.gt.maxzgrid) then
        write(*,*) ' #### FLEXPART MODEL ERROR! TOO MANY HEIGHT   #### ' 
        write(*,*) ' #### LEVELS ARE GIVEN FOR OUTPUT GRID.       #### ' 
        write(*,*) ' #### MAXIMUM NUMBER IS ',maxzgrid,'          #### ' 
        write(*,*) ' #### PLEASE MAKE CHANGES IN FILE OUTGRID.    #### ' 
        endif
        outheight(j)=outhelp
        goto 100
99    numzgrid=j-1


C Check whether vertical levels are specified in ascending order
****************************************************************

      do 40 j=2,numzgrid
        if (outheight(j).le.outheight(j-1)) then
        write(*,*) ' #### FLEXPART MODEL ERROR! YOUR SPECIFICATION#### ' 
        write(*,*) ' #### OF OUTPUT LEVELS IS CORRUPT AT LEVEL    #### ' 
        write(*,*) ' #### ',j,'                              #### ' 
        write(*,*) ' #### PLEASE MAKE CHANGES IN FILE OUTGRID.    #### ' 
        endif
40      continue

C Determine the half levels, i.e. middle levels of the output grid
******************************************************************

      outheighthalf(1)=outheight(1)/2.
      do 50 j=2,numzgrid
50      outheighthalf(j)=(outheight(j-1)+outheight(j))/2.


      xoutshift=xmet0-out_xm0
      youtshift=ymet0-out_ym0
      close(unitoutgrid)
      return


999   write(*,*) ' #### FLEXPART MODEL ERROR! FILE "OUTGRID"    #### '
      write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
      write(*,*) ' #### xxx/flexpart/options                    #### '
      stop

      end
