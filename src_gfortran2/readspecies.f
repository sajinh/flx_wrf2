      subroutine readspecies()
********************************************************************************
*                                                                              *
*     This routine reads names and physical constants of chemical species/     *
*     radionuclides available with FLEXPART.                                   *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     11 July 1996                                                             *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* decaytime(maxtable)  half time for radiological decay                        *
* specname(maxtable)   names of chemical species, radionuclides                *
* wetscava, wetscavb   Parameters for determining scavenging coefficient       *
*                                                                              *
* Constants:                                                                   *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      integer i


C Open the SPECIES file and read species names and properties
*************************************************************

      open(unitspecies,file=path(1)(1:len(1))//'SPECIES',status='old',
     +err=999)


      do 10 i=1,8
10      read(unitspecies,*)

      do 20 i=1,maxtable
        read(unitspecies,21,end=22) specname(i),decaytime(i),
     +  wetscava(i),wetscavb(i),drydiff(i),dryhenry(i),dryactiv(i),
     +  partrho(i),partmean(i),partsig(i),dryvelo(i),weightmol(i)

        if (partsig(i).eq.1.) partsig(i)=1.0001   ! avoid floating exception
        if (partsig(i).eq.0.) partsig(i)=1.0001   ! avoid floating exception

        if ((drydiff(i).gt.0.).and.(partrho(i).gt.0.)) then
          write(*,*) '#### FLEXPART MODEL ERROR! FILE "SPECIES"    ####'
          write(*,*) '#### IS CORRUPT. SPECIES CANNOT BE BOTH      ####'
          write(*,*) '#### PARTICLE AND GAS.                       ####'
          stop
        endif
20      continue

21    format(4x,a10,f10.1,e11.1,f6.2,f7.1,e9.1,f5.1,e10.1,2e8.1,2f8.2)


22    close(unitspecies)
      return


999   write(*,*) ' #### FLEXPART MODEL ERROR! FILE "SPECIES"    #### '
      write(*,*) ' #### CANNOT BE OPENED IN THE DIRECTORY       #### '
      write(*,'(a)') path(1)(1:len(1))
      stop
 
      end
