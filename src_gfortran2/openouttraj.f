      subroutine openouttraj()
********************************************************************************
*                                                                              *
*   Note:  This is the FLEXPART_WRF version of subroutine openouttraj.         *
*                                                                              *
*   This routine opens the output file for the plume trajectory output         *
*   produced by the cluster analysis.                                          *
*                                                                              *
*     Author: A. Stohl                                                         *
*     27 January 2001                                                          *
*                                                                              *
*     Dec 2005, R. Easter - changed names of "*lon0*" & "*lat0*" variables     *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'
      
      integer i
      real xp1,yp1,xp2,yp2


C Open output file for trajectory output
****************************************

      open(unitouttraj,file=path(2)(1:len(2))//'trajectories.txt',
     +form='formatted',err=998)

      if (ldirect.eq.1) then
      write(unitouttraj,'(i8,1x,i6,1x,a)') ibdate,ibtime,'FLEXPART V5.1'
      else
      write(unitouttraj,'(i8,1x,i6,1x,a)') iedate,ietime,'FLEXPART V5.1'
      endif
      write(unitouttraj,*) method,lsubgrid,lconvection
      write(unitouttraj,*) numpoint
      do 13 i=1,numpoint
        xp1=xpoint1(i)*dx+xmet0
        yp1=ypoint1(i)*dy+ymet0
        xp2=xpoint2(i)*dx+xmet0
        yp2=ypoint2(i)*dy+ymet0
        write(unitouttraj,*) ireleasestart(i),ireleaseend(i),
     +  xp1,yp1,xp2,yp2,zpoint1(i),zpoint2(i),kindz(i),npart(i)
        write(unitouttraj,'(a)') compoint(i)(1:40)
13      continue

      return

998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
      write(*,*) ' #### trajectories.txt                         #### '
      write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
      write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
      write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
      stop

      end
