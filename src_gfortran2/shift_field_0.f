      subroutine shift_field_0(field,nxf,nyf)
C                               i/o   i   i
********************************************************************************
*                                                                              *
*  This subroutine shifts global fields by nxshift grid cells, in order to     *
*  facilitate all sorts of nested wind fields, or output grids, which, without *
*  shifting, would overlap with the domain "boundary".                         *
*                                                                              *
*    Author: A. Stohl                                                          *
*                                                                              *
*    3 July 2002                                                               *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
*                                                                              *
* Constants:                                                                   *
*                                                                              *
********************************************************************************

      include 'includepar'

      integer nxf,nyf,ix,jy,ixs
      real field(0:nxmax-1,0:nymax-1),xshiftaux(0:nxmax-1)

C Loop over y and z
*******************

      do 10 jy=0,nyf-1

C Shift the data
****************

        if (nxshift.ne.0) then
          do 20 ix=0,nxf-1
            if (ix.ge.nxshift) then
              ixs=ix-nxshift
            else
              ixs=nxf-nxshift+ix
            endif
20          xshiftaux(ixs)=field(ix,jy)
          do 30 ix=0,nxf-1
30          field(ix,jy)=xshiftaux(ix)
        endif

C Repeat the westernmost grid cells at the easternmost domain "boundary"
************************************************************************

        field(nxf,jy)=field(0,jy)
10      continue

      return
      end
