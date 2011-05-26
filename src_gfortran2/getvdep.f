      subroutine getvdep(n,ix,jy,ust,temp,pa,L,gr,rh,rr,vdepo)
C                        i i  i   i   i   i  i i  i  i    o
********************************************************************************
*                                                                              *
*  This routine calculates the dry deposition velocities.                      *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     20 December 1996                                                         *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* gr [W/m2]         global radiation                                           *
* L [m]             Obukhov length                                             *
* nyl               kinematic viscosity                                        *
* pa [Pa]           surface air pressure                                       *
* ra [s/m]          aerodynamic resistance                                     *
* raquer [s/m]      average aerodynamic resistance                             *
* rh [0-1]          relative humidity                                          *
* rhoa              density of the air                                         *
* rr [mm/h]         precipitation rate                                         *
* temp [K]          2m temperature                                             *
* tc [C]            2m temperature                                             *
* ust [m/s]         friction velocity                                          *
*                                                                              *
* xlanduse          fractions of 9 landuses for each model grid point          *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      integer yyyymmdd,hhmmss,yyyy,mmdd,n,lseason,i,j,ix,jy
      real vdepo(maxspec),vd,rb(maxspec),rc(maxspec),raquer,eps
      real raerod,ra,ust,temp,tc,pa,L,gr,rh,rr,myl,nyl,rhoa,diffh2o
      parameter(eps=1.e-5)
      double precision jul



C Calculate month and determine the seasonal category
*****************************************************

      jul=bdate+dble(float(wftime(n))/86400.)

      call caldate(jul,yyyymmdd,hhmmss)
      yyyy=yyyymmdd/10000
      mmdd=yyyymmdd-10000*yyyy


      if ((mmdd.ge.1201).or.(mmdd.le.301)) then
        lseason=4
      else if ((mmdd.ge.1101).or.(mmdd.le.331)) then
        lseason=3
      else if ((mmdd.ge.401).and.(mmdd.le.515)) then
        lseason=5
      else if ((mmdd.ge.516).and.(mmdd.le.915)) then
        lseason=1
      else
        lseason=2
      endif
 

C Calculate diffusivity of water vapor
*************************************
 
      diffh2o=2.11e-5*(temp/273.15)**1.94*(101325/pa)
 
C Conversion of temperature from K to C
***************************************

      tc=temp-273.15

C Calculate dynamic viscosity
*****************************
 
      if (tc.lt.0) then
        myl=(1.718+0.0049*tc-1.2e-05*tc**2)*1.e-05
      else
        myl=(1.718+0.0049*tc)*1.e-05
      endif
 
C Calculate kinematic viscosity
*******************************
 
      rhoa=pa/(287.*temp)
      nyl=myl/rhoa
 

C 0. Set all deposition velocities zero
***************************************
 
      do 5 i=1,nspec
5       vdepo(i)=0.


C 1. Compute surface layer resistances rb
*****************************************

      call getrb(nspec,ust,nyl,diffh2o,reldiff,rb)

      raquer=0.
      do 10 j=1,numclass            ! loop over all landuse classes
 
        if (xlanduse(ix,jy,j).gt.eps)  then

C 2. Calculate aerodynamic resistance ra
****************************************

          ra=raerod(L,ust,z0(j))
          raquer=raquer+ra*xlanduse(ix,jy,j)

C 3. Calculate surface resistance for gases
*******************************************

          call getrc(nspec,lseason,j,tc,gr,rh,rr,rc)

C 4. Calculate deposition velocities for gases and ...
C 5. ... sum deposition velocities for all landuse classes
**********************************************************

          do 20 i=1,nspec
            if (reldiff(i).gt.0.) then
              if ((ra+rb(i)+rc(i)).gt.0.) then
                vd=1./(ra+rb(i)+rc(i))
              else
                vd=9.999
              endif
              vdepo(i)=vdepo(i)+vd*xlanduse(ix,jy,j)
            endif
20          continue
        endif
10      continue


C 6. Calculate deposition velocities for particles
**************************************************

      call partdep(nspec,density,fract,schmi,vset,raquer,ust,nyl,vdepo)


C 7. If no detailed parameterization available, take constant deposition
C    velocity if that is available
************************************************************************

      do 30 i=1,nspec
        if ((reldiff(i).lt.0.).and.(density(i).lt.0.).and.
     +  (dryvel(i).gt.0.)) then
          vdepo(i)=dryvel(i)
        endif
30      continue


      end
