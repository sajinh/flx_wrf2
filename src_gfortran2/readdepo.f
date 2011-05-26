      subroutine readdepo()   
*******************************************************************************
*                                                                             *
*  Reads dry deposition parameters needed by the procedure of Wesely (1989).  *
*  Wesely (1989): Parameterization of surface resistances to gaseous          *
*  dry deposition in regional-scale numerical models.                         *
*  Atmos. Environ. 23, 1293-1304.                                             *
*                                                                             *
*                                                                             *
*      AUTHOR: Andreas Stohl, 19 May 1995                                     *
*                                                                             *
*******************************************************************************
*                                                                             *
C Variables:                                                                  *
*                                                                             *
* rcl(maxspec,5,9) [s/m] Lower canopy resistance                              *
* rgs(maxspec,5,9) [s/m] Ground resistance                                    *
* rlu(maxspec,5,9) [s/m] Leaf cuticular resistance                            *
* rm(maxspec) [s/m]      Mesophyll resistance                                 *
* ri(maxspec) [s/m]      Stomatal resistance                                  *
*                                                                             *
c Constants:                                                                  *
*                                                                             *
*******************************************************************************


      include 'includepar'
      include 'includecom'

C FOR THIS SUBROUTINE, numclass=9 IS ASSUMED
********************************************

      real rluh(5,numclass),rgssh(5,numclass),rgsoh(5,numclass)
      real rclsh(5,numclass),rcloh(5,numclass)
      integer i,j,ic


C Read deposition constants related with landuse and seasonal category
**********************************************************************

      open(unitwesely,file=path(1)(1:len(1))//'surfdepo.t',
     +status='old',err=999)

      do 30 i=1,16
30      read(unitwesely,*)
      do 40 i=1,5
        read(unitwesely,*)
        read(unitwesely,'(8x,9f8.0)') (ri(i,j),j=1,numclass)
        read(unitwesely,'(8x,9f8.0)') (rluh(i,j),j=1,numclass)
        read(unitwesely,'(8x,9f8.0)') (rac(i,j),j=1,numclass)
        read(unitwesely,'(8x,9f8.0)') (rgssh(i,j),j=1,numclass)
        read(unitwesely,'(8x,9f8.0)') (rgsoh(i,j),j=1,numclass)
        read(unitwesely,'(8x,9f8.0)') (rclsh(i,j),j=1,numclass)
40      read(unitwesely,'(8x,9f8.0)') (rcloh(i,j),j=1,numclass)

      do 42 i=1,5
        do 42 j=1,numclass
          ri(i,j)=max(ri(i,j),0.001)
          rluh(i,j)=max(rluh(i,j),0.001)
          rac(i,j)=max(rac(i,j),0.001)
          rgssh(i,j)=max(rgssh(i,j),0.001)
          rgsoh(i,j)=max(rgsoh(i,j),0.001)
          rclsh(i,j)=max(rclsh(i,j),0.001)
42        rcloh(i,j)=max(rcloh(i,j),0.001)
      close(unitwesely)


C Compute additional parameters
*******************************

      do 50 ic=1,nspec
        if (reldiff(ic).gt.0.) then     ! gas is dry deposited
          do 52 i=1,5
            do 52 j=1,numclass
              rlu(ic,i,j)=rluh(i,j)/(1.e-5*henry(ic)+f0(ic))
              rgs(ic,i,j)=1./(henry(ic)/(10.e5*rgssh(i,j))+f0(ic)/
     +        rgsoh(i,j))
52            rcl(ic,i,j)=1./(henry(ic)/(10.e5*rclsh(i,j))+f0(ic)/
     +        rcloh(i,j))
        endif
50      continue


      return


999   write(*,*) '### FLEXPART ERROR! FILE              ###'
      write(*,*) '### surfdepo.t DOES NOT EXIST.        ###'
      stop

      end
