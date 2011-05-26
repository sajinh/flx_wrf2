      subroutine writeheader()
********************************************************************************
*                                                                              *
*  Note:  This is the FLEXPART_WRF version of subroutine writeheader.          *
*                                                                              *
*  This routine produces a file header containing basic information on the     *
*  settings of the FLEXPART run.                                               *
*  The header file is essential and must be read in by any postprocessing      *
*  program before reading in the output data.                                  *
*                                                                              *
*     Author: A. Stohl                                                         *
*     7 August 2002                                                            *
*                                                                              *
*     Dec 2005, J. Fast & R. Easter -                                          *
*            Write formatted output when iouttype=1.                           *
*            Write iomode_xycoord.  Write xy coords. as lat-lon or             *
*            grid-meters depending on iomode_xycoord.                          *
*            changed names of "*lon0*" & "*lat0*" variables                    *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
*                                                                              *
* xmet                   model x coordinate in grid-meters                     *
* ymet                   model y coordinate in grid-meters                     *
*                                                                              *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'
      
      integer jjjjmmdd,ihmmss,i,ix,jy,j
      real dxtmp,dytmp
      real xp1,yp1,xp2,yp2


*************************
C Open header output file
*************************

      if (iouttype .eq. 0) then
         open(unitheader,file=path(2)(1:len(2))//'header',
     +   form='unformatted',err=998)
      else
         open(unitheader,file=path(2)(1:len(2))//'header',
     +   form='formatted',err=998)
      endif


C Write the header information
******************************

      if (ldirect.eq.1) then
         if (iouttype .eq. 0) then
           write(unitheader) ibdate,ibtime,'FLEXPWRF V5.0'
         else
           write(unitheader,*) ibdate,ibtime
           write(unitheader,*) 'FLEXPWRF V5.0'
         endif
      else
         if (iouttype .eq. 0) then
           write(unitheader) iedate,ietime,'FLEXPWRF V5.0'
         else
           write(unitheader,*) iedate,ietime
           write(unitheader,*) 'FLEXPWRF V5.0'
         endif
      endif

C Write info on output interval, averaging time, sampling time
C FLEXPART_WRF - also the iomode_xycoord
**************************************************************

      if (iouttype .eq. 0) then
         write(unitheader) loutstep,loutaver,loutsample,iomode_xycoord
      else
         write(unitheader,*) loutstep,loutaver,loutsample,iomode_xycoord
      endif

C Write information on output grid setup
****************************************

      call caldate(bdate,jjjjmmdd,ihmmss)

      if (iomode_xycoord .eq. iomode_xycoord_latlon) then
         xp1 = outgrid_swlon
         yp1 = outgrid_swlat
         dxtmp = (outgrid_nelon-outgrid_swlon)/numxgrid
         dytmp = (outgrid_nelat-outgrid_swlat)/numygrid
      else
         xp1 = out_xm0
         yp1 = out_ym0
         dxtmp = dxout
         dytmp = dyout
      endif

      if (iouttype .eq. 0) then
         write(unitheader) xp1,yp1,numxgrid,numygrid,dxtmp,dytmp
         write(unitheader) numzgrid,(outheight(i),i=1,numzgrid)
         write(unitheader) jjjjmmdd,ihmmss
      else
         write(unitheader,*) xp1,yp1,numxgrid,numygrid,dxtmp,dytmp
         write(unitheader,*) numzgrid,(outheight(i),i=1,numzgrid)
         write(unitheader,*) jjjjmmdd,ihmmss
      endif

C Write number of species, and name for each species (+extra name for depositions)
C Indicate the dimension of the fields (i.e., 1 for deposition fields, numzgrid for
C concentration fields
***********************************************************************************

      if (iouttype .eq. 0) then
         write(unitheader) 3*nspec, numreceptor, nageclass
      else
         write(unitheader,*) 3*nspec, numreceptor, nageclass
      endif

      do 12 i=1,nspec
         if (iouttype .eq. 0) then
           write(unitheader) 1,'WD_'//species(i)(1:7)
           write(unitheader) 1,'DD_'//species(i)(1:7)
           write(unitheader) numzgrid,species(i)
         else
           write(unitheader,*) 1
           write(unitheader,*) 'WD_'//species(i)(1:7)
           write(unitheader,*) 1
           write(unitheader,*) 'DD_'//species(i)(1:7)
           write(unitheader,*) numzgrid
           write(unitheader,*) species(i)
         endif
12    continue

C Write information on release points: total number, then for each point:
C start, end, coordinates, # of particles, name, mass
*************************************************************************

      if (iouttype .eq. 0) then
         write(unitheader) numpoint
      else
         write(unitheader,*) numpoint
      endif

      do 13 i=1,numpoint

        if (iomode_xycoord .eq. iomode_xycoord_latlon) then
           xp1=releases_swlon(i)
           yp1=releases_swlat(i)
           xp2=releases_nelon(i)
           yp2=releases_nelat(i)
        else
           xp1=xpoint1(i)*dx+xmet0
           yp1=ypoint1(i)*dy+ymet0
           xp2=xpoint2(i)*dx+xmet0
           yp2=ypoint2(i)*dy+ymet0
        endif

        if (iouttype .eq. 0) then
           write(unitheader) ireleasestart(i),ireleaseend(i),kindz(i)
           write(unitheader) xp1,yp1,xp2,yp2,zpoint1(i),zpoint2(i)
           write(unitheader) npart(i),1
           write(unitheader) compoint(i)
           do j=1,nspec
             write(unitheader) xmass(i,j)
             write(unitheader) xmass(i,j)
             write(unitheader) xmass(i,j)
           enddo
        else
           write(unitheader,*) ireleasestart(i),ireleaseend(i),kindz(i)
           write(unitheader,*) xp1,yp1,xp2,yp2,zpoint1(i),zpoint2(i)
           write(unitheader,*) npart(i),1
           write(unitheader,*) compoint(i)
           do j=1,nspec
             write(unitheader,*) xmass(i,j)
             write(unitheader,*) xmass(i,j)
             write(unitheader,*) xmass(i,j)
           enddo
        endif

13      continue

C Write information on some model switches
******************************************

      if (iouttype .eq. 0) then
         write(unitheader) method,lsubgrid,lconvection
      else
         write(unitheader,*) method,lsubgrid,lconvection
      endif

C Write age class information
*****************************

      if (iouttype .eq. 0) then
         write(unitheader) nageclass,(lage(i),i=1,nageclass)
      else
         write(unitheader,*) nageclass,(lage(i),i=1,nageclass)
      endif


C Write topography to output file
*********************************

      do 30 ix=0,numxgrid-1
         if (iouttype .eq. 0) then
           write(unitheader) (oroout(ix,jy),jy=0,numygrid-1)
         else
           write(unitheader,*) (oroout(ix,jy),jy=0,numygrid-1)
         endif
30    continue


      close(unitheader)

      return


998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
      write(*,*) ' #### '//path(2)(1:len(2))//'header'//' #### '
      write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
      write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
      write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
      stop

      end
