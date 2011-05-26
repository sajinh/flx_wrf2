      subroutine readlanduse()                                
*******************************************************************************
*                                                                             *
*     Note:  This is the FLEXPART_WRF version of subroutine readland.         *
*            CURRENTLY IT DOES NOTHING as there is no point reading           *
*            the ECMWF landuse inventory.                                     *
*            Eventually this will read the WRF landuse info.                  *
*                                                                             *
*      Reads the landuse inventory into memory and relates it to Leaf Area    *
*      Index and roughness length.                                            *
*                                                                             *
*      AUTHOR: Andreas Stohl, 10 January 1994                                 *
*                                                                             *
*******************************************************************************
*                                                                             *
* Variables:                                                                  *
* i                       loop indices                                        *
* xlandinvent(555,224,9)  area fractions of 9 landuse categories              *
* LEN(numpath)            length of the path names                            *
* PATH(numpath)           contains the path names                             *
* unitland                unit connected with landuse inventory               *
*                                                                             *
* LANDUSE CATEGORIES:                                                         *
*                                                                             *
* 1 Arable land                                                               *
* 2 Grassland                                                                 *
* 3 Permanent crops                                                           *
* 4 Inland water                                                              *
* 5 urban areas                                                               *
* 6 Other                                                                     *
* 7 Forest                                                                    *
* 8 Forest                                                                    *
* 9 Ocean                                                                     *
*                                                                             *
*******************************************************************************

      include 'includepar'
      include 'includecom'

      integer i,ix,jy,k


      return


C Read landuse inventory
************************

      open(unitland,file=path(1)(1:len(1))//'landuse.asc',status='old',
     +err=998)
      do 10 ix=1,555
        do 10 jy=1,224
10        read(unitland,*) (xlandinvent(ix,jy,k),k=1,numclass) 
      close(unitland)  


C Read relation landuse-LAI,z0
******************************

      open(unitsurfdata,file=path(1)(1:len(1))//'surfdata.t',
     +status='old',err=999)

      do 30 i=1,4
30      read(unitsurfdata,*)
      do 40 i=1,numclass
40      read(unitsurfdata,'(21x,f5.2)') z0(i)
      close(unitsurfdata)

      return

C Issue error messages
**********************

998   write(*,*) ' #### FLEXPART ERROR! FILE CONTAINING          ####'
      write(*,*) ' #### LANDUSE INVENTORY DOES NOT EXIST         ####'
      stop

999   write(*,*) ' #### FLEXPART ERROR! FILE CONTAINING          ####'
      write(*,*) ' #### RELATION LANDUSE-LAI,z0 DOES NOT EXIST   ####'
      stop

      end
