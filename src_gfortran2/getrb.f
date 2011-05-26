      subroutine getrb(nc,ustar,nyl,diffh2o,reldiff,rb)   
c                      i    i    i     i       i    o
*******************************************************************************
*                                                                             *
*  Calculation of the quasilaminar sublayer resistance to dry deposition.     *
*                                                                             *
*      AUTHOR: Andreas Stohl, 20 May 1995                                     *
*                                                                             *
*******************************************************************************
*                                                                             *
C Variables:                                                                  *
* rb(ncmax)       sublayer resistance                                         *
* schmidt         Schmidt number                                              *
* ustar [m/s]     friction velocity                                           *
* diffh20 [m2/s]  diffusivity of water vapor in air                           *
* reldiff         diffusivity relative to H2O                                 *
*                                                                             *
C Constants:                                                                  *
* karman          von Karman constant                                         *
* pr              Prandtl number                                              *
*                                                                             *
*******************************************************************************


      include 'includepar'

      real ustar,diffh2o,rb(maxspec),schmidt,pr,nyl
      real reldiff(maxspec)
      integer ic,nc
      parameter(pr=0.72)

      do 10 ic=1,nc
        if (reldiff(ic).gt.0.) then
          schmidt=nyl/diffh2o*reldiff(ic)
          rb(ic)=2.0*(schmidt/pr)**0.67/(karman*ustar)
        endif
10      continue

      end
