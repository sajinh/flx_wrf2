      subroutine writeheader_nest()
********************************************************************************
*                                                                              *
*  Note:  This is the FLEXPART_WRF version of subroutine writeheder_nest.      *
*                                                                              *
*  This routine produces a file header containing basic information on the     *
*  settings of the FLEXPART run.                                               *
*  The header file is essential and must be read in by any postprocessing      *
*  program before reading in the output data.                                  *
*                                                                              *
*     Author: A. Stohl                                                         *
*     7 August 2002                                                            *
*                                                                              *
*     Dec 2005, R. Easter - changed names of "*lon0*" & "*lat0*" variables     *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
*                                                                              *
* xmet                   model x coordinate in grid-meters                     *
* ymet                   model x coordinate in grid-meters                     *
*                                                                              *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'
      
      integer jjjjmmdd,ihmmss,i,ix,jy,j
      real xp1,yp1,xp2,yp2


*************************
C Open header output file
*************************

      open(unitheader,file=path(2)(1:len(2))//'header_nest',
     +form='unformatted',err=998)


C Write the header information
******************************

      if (ldirect.eq.1) then
        write(unitheader) ibdate,ibtime,'FLEXPART V5.0'
      else
        write(unitheader) iedate,ietime,'FLEXPART V5.0'
      endif

C Write info on output interval, averaging time, sampling time
**************************************************************

      write(unitheader) loutstep,loutaver,loutsample

C Write information on output grid setup
****************************************

      write(unitheader) out_xm0n,out_ym0n,numxgridn,numygridn,
     +dxoutn,dyoutn
      write(unitheader) numzgrid,(outheight(i),i=1,numzgrid)

      call caldate(bdate,jjjjmmdd,ihmmss)
      write(unitheader) jjjjmmdd,ihmmss

C Write number of species, and name for each species (+extra name for depositions)
C Indicate the dimension of the fields (i.e., 1 for deposition fields, numzgrid for
C concentration fields
***********************************************************************************

      write(unitheader) 3*nspec
      do 12 i=1,nspec
        write(unitheader) 1,'WD_'//species(i)(1:7)
        write(unitheader) 1,'DD_'//species(i)(1:7)
12      write(unitheader) numzgrid,species(i)

C Write information on release points: total number, then for each point:
C start, end, coordinates, # of particles, name, mass
*************************************************************************

      write(unitheader) numpoint
      do 13 i=1,numpoint
        write(unitheader) ireleasestart(i),ireleaseend(i),kindz(i)
        xp1=xpoint1(i)*dx+xmet0
        yp1=ypoint1(i)*dy+ymet0
        xp2=xpoint2(i)*dx+xmet0
        yp2=ypoint2(i)*dy+ymet0
        write(unitheader) xp1,yp1,xp2,yp2,zpoint1(i),zpoint2(i)
        write(unitheader) npart(i),1
        write(unitheader) compoint(i)
        do 13 j=1,nspec
          write(unitheader) xmass(i,j)
          write(unitheader) xmass(i,j)
13        write(unitheader) xmass(i,j)

C Write information on some model switches
******************************************

      write(unitheader) method,lsubgrid,lconvection

C Write age class information
*****************************

      write(unitheader) nageclass,(lage(i),i=1,nageclass)


C Write topography to output file
*********************************

      do 30 ix=0,numxgridn-1
30      write(unitheader) (orooutn(ix,jy),jy=0,numygridn-1)
      close(unitheader)

      return


998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
      write(*,*) ' #### '//path(2)(1:len(2))//'header'//' #### '
      write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
      write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
      write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
      stop

      end
