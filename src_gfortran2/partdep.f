      subroutine partdep(nc,density,fract,schmi,vset,ra,ustar,nyl,vdep)
c                        i     i      i     i    i   i    i    i  i/o
*******************************************************************************
*                                                                             *
*      Calculation of the dry deposition velocities of particles.             *
*                                                                             *
*     AUTHOR: Andreas Stohl, 12 November 1993                                 *
*                            Update: 20 December 1996                         *
*                                                                             *
*                                                                             *
*     Literature:                                                             *
*     [1]  Hicks/Baldocchi/Meyers/Hosker/Matt (1987), A Preliminary           *
*             Multiple Resistance Routine for Deriving Dry Deposition         *
*             Velocities from Measured Quantities.                            *
*             Water, Air and Soil Pollution 36 (1987), pp.311-330.            *
*     [2]  Slinn (1982), Predictions for Particle Deposition to               *
*             Vegetative Canopies. Atm.Env.16-7 (1982), pp.1785-1794.         *
*     [3]  Slinn/Slinn (1980),  Predictions for Particle Deposition on        *
*             Natural Waters. Atm.Env.14 (1980), pp.1013-1016.                *
*     [4]  Scire/Yamartino/Carmichael/Chang (1989),                           *
*             CALGRID: A Mesoscale Photochemical Grid Model.                  *
*             Vol II: User's Guide. (Report No.A049-1, June, 1989)            *
*     [5]  Langer M. (1992): Ein einfaches Modell zur Abschaetzung der        *
*             Depositionsgeschwindigkeit von Teilchen und Gasen.              *
*             Internal report.                                                *
*******************************************************************************
*                                                                             *
c Variables:                                                                  *
* alpha                help variable                                          *
* fract(nc,ni)         mass fraction of each diameter interval                *
* lpdep(nc)            1 for particle deposition, 0 else                      *
* nc                   actual number of chemical components                   *
* ni                   number of diameter intervals, for which vdepj is calc. *
* rdp [s/m]            deposition layer resistance                            *
* ra [s/m]             aerodynamical resistance                               *
* schmi(nc,ni)         Schmidt number**2/3 of each diameter interval          *
* stokes               Stokes number                                          *
* ustar [m/s]          friction velocity                                      *
* vdep(nc) [m/s]       deposition velocities of all components                *
* vdepj [m/s]          help, deposition velocity of 1 interval                *
* vset(nc,ni)          gravitational settling velocity of each interval       *
*                                                                             *
c Constants:                                                                  *
* nc                   number of chemical species                             *
* ni                   number of diameter intervals, for which deposition     *
*                      is calculated                                          *
*                                                                             *
*******************************************************************************

      include 'includepar'

      real density(maxspec),schmi(maxspec,ni),fract(maxspec,ni)
      real vset(maxspec,ni)
      real vdep(maxspec),stokes,vdepj,rdp,ustar,alpha,ra,nyl,eps
      parameter(eps=1.e-5)
      integer ic,j,nc

 
      do 90 ic=1,nc                  ! loop over all species
        if (density(ic).gt.0.) then
          do 100 j=1,ni              ! loop over all diameter intervals
            if (ustar.gt.eps) then

C Stokes number for each diameter interval
******************************************

              stokes=vset(ic,j)/ga*ustar*ustar/nyl
              alpha=-3./stokes

C Deposition layer resistance
*****************************
 
              if (alpha.le.log10(eps)) then
                rdp=1./(schmi(ic,j)*ustar)
              else
               	rdp=1./((schmi(ic,j)+10.**alpha)*ustar)
              endif
              vdepj=vset(ic,j)+1./(ra+rdp+ra*rdp*vset(ic,j))
            else
              vdepj=vset(ic,j)
            endif

C deposition velocities of each interval are weighted with mass fraction
************************************************************************

100         vdep(ic)=vdep(ic)+vdepj*fract(ic,j)
        endif
90      continue

      end
