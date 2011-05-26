      real function qak_distance2(y1,x1,y2,x2)
C                             i  i  i  i
********************************************************************************
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of file distance2.f              *
*                                                                              *
*     Calculation of the distance (in km) between two points on the Earth's    *
*     surface. Input are latitudes and longitudes of the two points: (y1,x1)   *
*     and (y2,x2), which in this routine are required in radians, NOT in degree*
*     as is the case in distance.f                                             *
*     Output is the distance in km.                                            *
*                                                                              *
*     Author: G. Wotawa, adapted and documented by A. Stohl                    *
*                                                                              *
*     24 January 2002                                                          *
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of file distance.f               *
*     26 Oct 2005, R. Easter -                                                 *
*          The function name was changed from "distance2" to "qak_distance2"   *
*          to eliminate all usages of the "distance2" function in FLEXPART_WRF *
*                                                                              *
********************************************************************************

      include 'includepar'

      real y1,x1,y2,x2,ra
      double precision y10,x10,y20,x20,dl,dfl8
      parameter(ra=r_earth/1000.)                ! Earth's radius in km

      y10=y1
      y20=y2
      x10=x1
      x20=x2
      dl=(x10-x20)
      dfl8=1./(dcos(y10)*dcos(y20)*dcos(dl)+dsin(y10)*dsin(y20))**2-1.
      if((dfl8.lt.0.).and.(dfl8.gt.-0.0001)) dfl8=0.
      qak_distance2=ra*datan(dsqrt(dfl8))

      end

