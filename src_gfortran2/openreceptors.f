      subroutine openreceptors()
********************************************************************************
*                                                                              *
*  Note:  This is the FLEXPART_WRF version of subroutine openreceptors.        *
*                                                                              *
*  This routine opens the receptor output files and writes out the receptor    *
*  names and the receptor locations. The receptor output files are not closed, *
*  but kept open throughout the simulation. Concentrations are continuously    *
*  dumped to these files.                                                      *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     7 August 2002                                                            *
*                                                                              *
*     Dec 2005, J. Fast - Output files can be either binary or ascii.          *
*                         Write iomode_xycoord to output files.                *
*                         Receptor positions can be lat-lon or grid-meters.    *
*     Dec 2005, R. Easter - changed names of "*lon0*" & "*lat0*" variables     *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* numreceptor            actual number of receptor points specified            *
* receptornames          names of the receptor points                          *
* xreceptor,yreceptor    coordinates of the receptor points                    *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'
      
      integer j
      real xtmp(maxreceptor),ytmp(maxreceptor)
      real xtmpb,ytmpb


C Open output file for receptor points and write out a short header
C containing receptor names and locations
*******************************************************************

      if (numreceptor.ge.1) then           ! do it only if receptors are specified

        do j = 1, numreceptor
            xtmp(j) = xreceptor(j)*dx + xmet0
            ytmp(j) = yreceptor(j)*dy + ymet0
            if (iomode_xycoord .eq. iomode_xycoord_latlon) then
               xtmpb = xtmp(j)
               ytmpb = ytmp(j)
               call xymeter_to_ll_wrf( xtmpb, ytmpb, xtmp(j), ytmp(j) )
            endif
        enddo

C Concentration output
**********************

        if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
          if (iouttype.eq.0) then 
          open(unitoutrecept,file=path(2)(1:len(2))//'receptor_conc',
     +    form='unformatted',err=997)
          write(unitoutrecept) (receptorname(j),j=1,numreceptor)
          write(unitoutrecept) (xtmp(j),ytmp(j),j=1,numreceptor),
     +       iomode_xycoord
          endif
          if (iouttype.eq.1) then 
          open(unitoutrecept,file=path(2)(1:len(2))//'receptor_conc',
     +    form='formatted',err=997)
          do j = 1, numreceptor
            write(unitoutrecept,*) receptorname(j)
          enddo
          write(unitoutrecept,*) (xtmp(j),ytmp(j),j=1,numreceptor),
     +       iomode_xycoord
          endif
        endif

C Mixing ratio output
*********************

        if ((iout.eq.2).or.(iout.eq.3)) then
          if (iouttype.eq.0) then 
          open(unitoutreceptppt,file=path(2)(1:len(2))//'receptor_pptv',
     +    form='unformatted',err=998)
          write(unitoutreceptppt) (receptorname(j),j=1,numreceptor)
          write(unitoutreceptppt) (xtmp(j),ytmp(j),j=1,numreceptor),
     +       iomode_xycoord
          endif
          if (iouttype.eq.1) then 
          open(unitoutreceptppt,file=path(2)(1:len(2))//'receptor_pptv',
     +    form='formatted',err=998)
          do j = 1, numreceptor
            write(unitoutreceptppt,*) receptorname(j)
          enddo
          write(unitoutreceptppt,*) (xtmp(j),ytmp(j),j=1,numreceptor),
     +       iomode_xycoord
          endif
        endif
      endif

      return


997   write(*,*) ' #### FLEXPART MODEL ERROR! THE FILE           #### '
      write(*,*) ' ####              receptor_conc               #### '
      write(*,*) ' #### CANNOT BE OPENED.                        #### '
      stop

998   write(*,*) ' #### FLEXPART MODEL ERROR! THE FILE           #### '
      write(*,*) ' ####              receptor_pptv               #### '
      write(*,*) ' #### CANNOT BE OPENED.                        #### '
      stop

      end
