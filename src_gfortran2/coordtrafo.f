      subroutine coordtrafo()
***********************************************************************
*                                                                     * 
* Note:  This is the FLEXPART_WRF version of subroutine coordtrafo.   *
*                                                                     * 
*             FLEXPART MODEL SUBROUTINE COORDTRAFO                    *
*                                                                     *
***********************************************************************
*                                                                     * 
* AUTHOR:      G. WOTAWA                                              *
* DATE:        1994-02-07                                             *
* LAST UPDATE: 1996-05-18   A. STOHL                                  *
*                                                                     * 
* Dec 2005, R. Easter - changed names of "*lon0*" & "*lat0*" variables*
*                                                                     * 
***********************************************************************
*                                                                     *
* DESCRIPTION: This subroutine transforms x and y coordinates of      *
* particle release points to grid coordinates.                        *
*                                                                     *
***********************************************************************

      include 'includepar'
      include 'includecom'

      integer i,j

      if (numpoint.eq.0) goto 30

* TRANSFORM X- AND Y- COORDINATES OF STARTING POINTS TO GRID COORDINATES
************************************************************************

      do 10 i=1,numpoint
        xpoint1(i)=(xpoint1(i)-xmet0)/dx
        xpoint2(i)=(xpoint2(i)-xmet0)/dx
        ypoint1(i)=(ypoint1(i)-ymet0)/dy
10      ypoint2(i)=(ypoint2(i)-ymet0)/dy

15    continue


* CHECK IF RELEASE POINTS ARE WITHIN DOMAIN
*******************************************

      do 25 i=1,numpoint
        if (sglobal.and.(ypoint1(i).lt.1.e-6)) ypoint1(i)=1.e-6
        if (nglobal.and.(ypoint2(i).gt.float(nymin1)-1.e-5))
     +  ypoint2(i)=float(nymin1)-1.e-5
      if ((ypoint1(i).lt.1.e-6).or.(ypoint1(i).ge.float(nymin1)-1.e-6)
     +.or.(ypoint2(i).lt.1.e-6).or.(ypoint2(i).ge.float(nymin1)-1.e-6)
     +.or.((.not.xglobal).and.((xpoint1(i).lt.1.e-6).or.
     +(xpoint1(i).ge.float(nxmin1)-1.e-6).or.(xpoint2(i).lt.1.e-6).or.
     +(xpoint2(i).ge.float(nxmin1)-1.e-6)))) then
          write(*,*) ' NOTICE: RELEASE POINT OUT OF DOMAIN DETECTED.'
          write(*,*) ' IT IS REMOVED NOW ... '
          write(*,*) ' COMMENT: ',compoint(i)

          if (i.lt.numpoint) then
            do 20 j=i+1,numpoint
              xpoint1(j-1)=xpoint1(j)
              ypoint1(j-1)=ypoint1(j)
              xpoint2(j-1)=xpoint2(j)
              ypoint2(j-1)=ypoint2(j)
              zpoint1(j-1)=zpoint1(j)
              zpoint2(j-1)=zpoint2(j)
              npart(j-1)=npart(j)
20            compoint(j-1)=compoint(j)         
          endif

          numpoint=numpoint-1
          if (numpoint.gt.0) goto 15
        endif
25      continue

30    if (numpoint.eq.0) then
        write(*,*) ' FLEXPART MODEL SUBROUTINE COORDTRAFO: ERROR ! '
        write(*,*) ' NO PARTICLE RELEASES ARE DEFINED!'
        write(*,*) ' CHECK FILE RELEASES...'
        stop
      endif

      end 
