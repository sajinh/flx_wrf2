      subroutine partoutput_short(itime)
C                                   i
********************************************************************************
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subroutine assignland.        *
*                                                                              *
*     Dump all particle positions                                              *
*                                                                              *
*     Author: A. Stohl                                                         *
*     12 March 1999                                                            *
*                                                                              *
*     Dec 2005, R. Easter - Made this routine inactive.                        *
*                           Changed names of "*lon0*" & "*lat0*" variables     *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      double precision jul
      integer itime,i,j,jjjjmmdd,ihmmss,numshortout,numshortall
      integer ix,jy,ixp,jyp
      real xlon,ylat,zlim,dt1,dt2,dtt,ddx,ddy,rddx,rddy,p1,p2,p3,p4,topo
      character adate*8,atime*6

      integer*2 idump(3,maxpart)
      integer i4dump(maxpart)


c FLEXPART_WRF - currently this routine is inactive
      write(*,'(/a/)') 
     &   '*** partoutput_short is not implemented for FLEXPART_WRF'
      return


C Determine current calendar date, needed for the file name
***********************************************************

      jul=bdate+dble(float(itime))/86400.
      call caldate(jul,jjjjmmdd,ihmmss)
      write(adate,'(i8.8)') jjjjmmdd
      write(atime,'(i6.6)') ihmmss


C Some variables needed for temporal interpolation
**************************************************

      dt1=float(itime-memtime(1))
      dt2=float(memtime(2)-itime)
      dtt=1./(dt1+dt2)


C Loop about all particles
**************************

      numshortout=0
      numshortall=0
      do 10 i=1,numpart

C Take only valid particles
***************************

        if (itra1(i).eq.itime) then
          xlon=xmet0+xtra1(i)*dx
          ylat=ymet0+ytra1(i)*dy

**********************************************************************************
C Interpolate several variables (PV, specific humidity, etc.) to particle position
**********************************************************************************

          ix=xtra1(i)
          jy=ytra1(i)
          ixp=ix+1
          jyp=jy+1
          ddx=xtra1(i)-float(ix)
          ddy=ytra1(i)-float(jy)
          rddx=1.-ddx
          rddy=1.-ddy
          p1=rddx*rddy
          p2=ddx*rddy
          p3=rddx*ddy
          p4=ddx*ddy

C Topography
************

          topo=p1*oro(ix ,jy)
     +       + p2*oro(ixp,jy)
     +       + p3*oro(ix ,jyp)
     +       + p4*oro(ixp,jyp)


C Convert positions to integer*2 variables (from -32768 to 32767)
C Do this only for region of main interest, i.e. extended North Atlantic region,
C and for the tracer of interest, i.e. the North American one
********************************************************************************

          if (xlon.gt.180.) xlon=xlon-360.
          if (xlon.lt.-180.) xlon=xlon+360.

          numshortall=numshortall+1
          if ((xlon.gt.-140).and.(xlon.lt.60).and.(ylat.gt.10).and.
     +    (xmass1(i,1).gt.0.)) then
            numshortout=numshortout+1
            idump(1,numshortout)=nint(xlon*180.)
            idump(2,numshortout)=nint(ylat*360.)
            zlim=min(ztra1(i)+topo,32766.)
            idump(3,numshortout)=nint(zlim)
            i4dump(numshortout)=npoint(i)
          endif

        endif
10      continue


C Open output file and write the output
***************************************

      open(unitshortpart,file=path(2)(1:len(2))//'shortposit_'//adate//
     +atime,form='unformatted')

C Write current time to file
****************************

      write(unitshortpart) itime
      write(unitshortpart) numshortout
      write(unitshortpart)
     +(i4dump(i),(idump(j,i),j=1,3),i=1,numshortout)


      write(*,*) numshortout,numshortall

      close(unitshortpart)

      end
