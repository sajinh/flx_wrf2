      subroutine shift_field(field,nxf,nyf,nzfmax,nzf,nmax,n)
C                             i/o   i   i    i     i   i   i
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

      integer nxf,nyf,nzf,n,ix,jy,kz,ixs,nzfmax,nmax
      real field(0:nxmax-1,0:nymax-1,nzfmax,nmax),xshiftaux(0:nxmax-1)

C Loop over y and z
*******************

      do 10 kz=1,nzf
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
20            xshiftaux(ixs)=field(ix,jy,kz,n)
            do 30 ix=0,nxf-1
30            field(ix,jy,kz,n)=xshiftaux(ix)
          endif

C Repeat the westernmost grid cells at the easternmost domain "boundary"
************************************************************************

          field(nxf,jy,kz,n)=field(0,jy,kz,n)
10        continue

      end
