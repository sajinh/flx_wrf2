      subroutine centerofmass(xl,yl,n,xcenter,ycenter)
C                             i  i  i    o       o
********************************************************************************
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subroutine assignland.        *
*            The computational grid is the WRF x-y grid rather than lat-lon.   *
*                                                                              *
*   This routine calculates the center of mass of n points on the Earth.       *
*   Input are the longitudes (xl) and latitudes (yl) of the individual         *
*   points, output is the longitude and latitude of the centre of mass.        *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     24 January 2002                                                          *
*                                                                              *
*    26 Oct 2005, R. Easter - changes associated with WRF horizontal grid.     *
*                             Since x & y coords are in meters,                *
*                             so just sum/average the xl & yl.                 *
*                                                                              *
********************************************************************************

      include 'includepar'

      integer n,l
      real xl(n),yl(n),xll,yll,xav,yav,zav,x,y,z,xcenter,ycenter


      xav=0.
      yav=0.
      zav=0.

      do 10 l=1,n

C Convert longitude and latitude from degrees to radians
********************************************************

c for FLEXPART_WRF, x & y coords are in meters, 
c so just sum/average the xl & yl
c       xll=xl(l)*pi180
c       yll=yl(l)*pi180

C Calculate 3D coordinates from longitude and latitude
******************************************************

c for FLEXPART_WRF, this isn't necessary
c       x = cos(yll)*sin(xll)
c       y = -1.*cos(yll)*cos(xll)
c       z = sin(yll)


C Find the mean location in Cartesian coordinates
*************************************************

c       xav=xav+x
c       yav=yav+y
c       zav=zav+z
        xav=xav+xl(l)
        yav=yav+yl(l)

10      continue

      xav=xav/float(n)
      yav=yav/float(n)
c     zav=zav/float(n)


C Project the point back onto Earth's surface
*********************************************

c for FLEXPART_WRF, this isn't necessary
c     xcenter=atan2(xav,-1.*yav)
c     ycenter=atan2(zav,sqrt(xav*xav+yav*yav))

C Convert back to degrees
*************************

c for FLEXPART_WRF, this isn't necessary
c     xcenter=xcenter/pi180
c     ycenter=ycenter/pi180
      xcenter=xav
      ycenter=yav

      end
