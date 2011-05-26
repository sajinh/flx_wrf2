      subroutine readreceptors()
********************************************************************************
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subroutine assignland.        *
*            The computational grid is the WRF x-y grid rather than lat-lon.   *
*                                                                              *
*     This routine reads the user specifications for the receptor points.      *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     1 August 1996                                                            *
*                                                                              *
*     Oct 2005, R. Easter - change calc of receptorarea()                      *
*     Dec 2005, R. Easter - x/yrecptor values may be input either as           *
*                           degrees-latlon or grid-meters                      *
*                           Changes to some read formats (wider fields).       *
*                           Changed names of "*lon0*" & "*lat0*" variables     *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* receptorarea(maxreceptor)  area of dx*dy at location of receptor             *
* receptorname(maxreceptor)  names of receptors                                *
* xreceptor,yreceptor  coordinates of receptor points                          *
*                                                                              *
* Constants:                                                                   *
* unitreceptor         unit connected to file RECEPTORS                        *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      integer j
      real x,y,xm,ym
      character*16 receptor


C For backward runs, do not allow receptor output. Thus, set number of receptors to zero
****************************************************************************************

      if (ldirect.lt.0) then
        numreceptor=0
        return
      endif
        

C Open the RECEPTORS file and read output grid specifications
*************************************************************

      open(unitreceptor,file=path(1)(1:len(1))//'RECEPTORS',
     +status='old',err=999)

      call skplin(5,unitreceptor)


C Read the names and coordinates of the receptors
*************************************************

      j=0
100     j=j+1
        read(unitreceptor,*,end=99)
        read(unitreceptor,*,end=99)
        read(unitreceptor,*,end=99)
        read(unitreceptor,'(4x,a16)',end=99) receptor
        call skplin(3,unitreceptor)
        read(unitreceptor,'(f15.8)',end=99) x
        call skplin(3,unitreceptor)
        read(unitreceptor,'(f15.8)',end=99) y
        if ((x.eq.0.).and.(y.eq.0.).and.
     +  (receptor.eq.'                ')) then
          j=j-1
          goto 100
        endif
        if (j.gt.maxreceptor) then
        write(*,*) ' #### FLEXPART MODEL ERROR! TOO MANY RECEPTOR #### ' 
        write(*,*) ' #### POINTS ARE GIVEN.                       #### ' 
        write(*,*) ' #### MAXIMUM NUMBER IS ',maxreceptor,'       #### ' 
        write(*,*) ' #### PLEASE MAKE CHANGES IN FILE RECEPTORS   #### ' 
        endif
        receptorname(j)=receptor

        write(*,'(/a,i5)') 'readreceptor diagnostics, j =', j
        write(*,'(a,1p,2e18.10)') 'x, y (in)   ', x, y
        if (iomode_xycoord .eq. iomode_xycoord_latlon) then
C In this case, the above inputs are the actual geographical lat/lon
C Need to convert from lat/lon to grid-index coordinates
           receptor_lon(j) = x
           receptor_lat(j) = y
           call ll_to_xyindex_wrf( receptor_lon(j), receptor_lat(j),
     &        xreceptor(j), yreceptor(j) )
        else
C In this case, the above inputs are in grid-meters 
C Need to convert from grid-meters to grid-index coordinates, then to lat/lon
           xreceptor(j)=(x-xmet0)/dx
           yreceptor(j)=(y-ymet0)/dy
           call xyindex_to_ll_wrf( 0, xreceptor(j), yreceptor(j),
     &        receptor_lon(j), receptor_lat(j) )
           write(*,'(f15.10,5x,a)') receptor_lon(j), 'receptor_lon'
           write(*,'(f15.10,5x,a)') receptor_lat(j), 'receptor_lat'
        end if
        write(*,'(a,1p,2e18.10)') 'x, yreceptor', 
     &     xreceptor(j), yreceptor(j)

c for FLEXPART_WRF, dx & dy are in meters,
c receptorarea() appears to be area in m**2 of a mother grid cell
c centers on x,y
c       xm=r_earth*cos(y*pi/180.)*dx/180.*pi
c       ym=r_earth*dy/180.*pi
        xm=dx
        ym=dy
        receptorarea(j)=xm*ym
        goto 100

99    numreceptor=j-1
      close(unitreceptor)
      return


999   write(*,*) ' #### FLEXPART MODEL ERROR! FILE "RECEPTORS"  #### '
      write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
      write(*,'(a)') path(1)(1:len(1))
      stop

      end
