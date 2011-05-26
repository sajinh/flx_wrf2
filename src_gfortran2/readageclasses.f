      subroutine readageclasses()
********************************************************************************
*                                                                              *
*     This routine reads the age classes to be used for the current model run. *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     20 March 2000                                                            *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
*                                                                              *
* Constants:                                                                   *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      integer i


C If age spectra calculation is switched off, set number of age classes
C to 1 and maximum age to a large number
***********************************************************************

      if (lagespectra.ne.1) then
        nageclass=1
        lage(nageclass)=999999999
        return
      endif


C If age spectra claculation is switched on,
C open the AGECLASSSES file and read user options
*************************************************

      open(unitageclasses,file=path(1)(1:len(1))//'AGECLASSES',
     +status='old',err=999)

      do 10 i=1,13
10      read(unitageclasses,*)
      read(unitageclasses,*) nageclass
      

      if (nageclass.gt.maxageclass) then
        write(*,*) ' #### FLEXPART MODEL ERROR! NUMBER OF AGE     #### ' 
        write(*,*) ' #### CLASSES GREATER THAN MAXIMUM ALLOWED.   #### '
        write(*,*) ' #### CHANGE SETTINGS IN FILE AGECLASSES OR   #### '
        write(*,*) ' #### RECOMPILE WITH LARGER MAXAGECLASS IN    #### '
        write(*,*) ' #### FILE INCLUDEPAR.                        #### '
        stop
      endif

      read(unitageclasses,*) lage(1)
      if (lage(1).le.0) then
        write(*,*) ' #### FLEXPART MODEL ERROR! AGE OF FIRST      #### ' 
        write(*,*) ' #### CLASS MUST BE GREATER THAN ZERO. CHANGE #### ' 
        write(*,*) ' #### SETTINGS IN FILE AGECLASSES.            #### ' 
        stop
      endif
 
      do 20 i=2,nageclass
        read(unitageclasses,*) lage(i)
        if (lage(i).le.lage(i-1)) then
          write(*,*) ' #### FLEXPART MODEL ERROR! AGE CLASSES     #### ' 
          write(*,*) ' #### MUST BE GIVEN IN TEMPORAL ORDER.      #### ' 
          write(*,*) ' #### CHANGE SETTINGS IN FILE AGECLASSES.   #### ' 
          stop
        endif
20      continue

      return    

999   write(*,*) ' #### FLEXPART MODEL ERROR! FILE "AGECLASSES" #### ' 
      write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
      write(*,'(a)') path(1)(1:len(1))
      stop

      end
