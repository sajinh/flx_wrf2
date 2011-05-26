      subroutine skplin(nlines,iunit)
C                         i      i
********************************************************************************
*                                                                              *
*     This routine reads nlines from unit iunit and discards them
*                                                                              *
*     Authors: Petra Seibert
*                                                                              *
*     31 Dec 1998
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:     
*
* iunit   unit number from which lines are to be skipped
* nlines  number of lines to be skipped
*
********************************************************************************

      integer i,iunit, nlines
      
      do 10 i=1,nlines
10      read(iunit,*)

      end
