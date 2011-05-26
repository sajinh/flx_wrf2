      real function qak_distance(y1,x1,y2,x2)
C                            i  i  i  i
********************************************************************************
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of file distance.f               *
*                                                                              *
*     Calculation of the distance (in km) between two points on the Earth's    *
*     surface. Input are latitudes and longitudes of the two points: (y1,x1)   *
*     and (y2,x2), output is the distance in km.                               *
*                                                                              *
*     Author: G. Wotawa, adapted and documented by A. Stohl                    *
*                                                                              *
*     24 January 2002                                                          *
*                                                                              *
*     26 Oct 2005, R. Easter -                                                 *
*          The function name was changed from "distance" to "qak_distance"     *
*          to eliminate all usages of the "distance" function in FLEXPART_WRF  *
*                                                                              *
********************************************************************************

      include 'includepar'

      real y1,x1,y2,x2,ra
      double precision y10,x10,y20,x20,dl,dfl8
      parameter(ra=r_earth/1000.)                ! Earth's radius in km

      y10=y1*pi180
      y20=y2*pi180
      x10=x1
      x20=x2
      dl=(x10-x20)*pi180
      dfl8=1./(dcos(y10)*dcos(y20)*dcos(dl)+dsin(y10)*dsin(y20))**2-1.
      if((dfl8.lt.0.).and.(dfl8.gt.-0.0001)) dfl8=0.
      qak_distance=ra*datan(dsqrt(dfl8))

      end

