      subroutine readpaths()
********************************************************************************
*                                                                              *
*     Reads the pathnames, where input/output files are expected to be.        *
*     The file pathnames must be available in the current working directory.   *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     1 February 1994                                                          *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* len(numpath)       lengths of the path names                                 *
* path(numpath)      pathnames of input/output files                           *
*                                                                              *
* Constants:                                                                   *
* numpath            number of pathnames to be read in                         *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      integer i

C Read the pathname information stored in unitpath
**************************************************

      open(unitpath,file='pathnames',status='old',err=999)

      do 10 i=1,numpath
        read(unitpath,'(a)',err=998) path(i) 
10      len(i)=index(path(i),' ')-1

C Check whether any nested subdomains are to be used
****************************************************

      do 20 i=1,maxnests
        read(unitpath,'(a)') path(numpath+2*(i-1)+1) 
        read(unitpath,'(a)') path(numpath+2*(i-1)+2) 
        if (path(numpath+2*(i-1)+1)(1:5).eq.'=====') goto 30
        len(numpath+2*(i-1)+1)=index(path(numpath+2*(i-1)+1),' ')-1
20      len(numpath+2*(i-1)+2)=index(path(numpath+2*(i-1)+2),' ')-1


C Determine number of available nested domains
**********************************************

30    numbnests=i-1


      close(unitpath)
      return    

998   write(*,*) ' #### TRAJECTORY MODEL ERROR! ERROR WHILE     #### ' 
      write(*,*) ' #### READING FILE PATHNAMES.                 #### ' 
      stop

999   write(*,*) ' #### TRAJECTORY MODEL ERROR! FILE "pathnames"#### ' 
      write(*,*) ' #### CANNOT BE OPENED IN THE CURRENT WORKING #### '
      write(*,*) ' #### DIRECTORY.                              #### '
      stop

      end
