      subroutine richardson(psurf,ust,ttlev,qvlev,ulev,vlev,nuvz,
     +  pplev,hf,tt2,td2,h,wst,hmixplus,ierr,sfc_option)
c     subroutine richardson(psurf,ust,ttlev,qvlev,ulev,vlev,nuvz,
c    +akz,bkz,hf,tt2,td2,h,wst,hmixplus,ierr)

C                             i    i    i     i    i    i    i
C      i   i  i   i   i  o  o     o
*****************************************************************************
*                                                                           *
*     Note:  This is the FLEXPART_WRF version of subroutine richardson.     *
*                                                                           *
*     Calculation of mixing height based on the critical Richardson number. *
*     Calculation of convective time scale.                                 *
*     For unstable conditions, one iteration is performed. An excess        *
*     temperature (dependent on hf and wst) is calculated, added to the     *
*     temperature at the lowest model level. Then the procedure is repeated.*
*                                                                           *
*     Author: A. Stohl                                                      *
*                                                                           *
*     22 August 1996                                                        *
*                                                                           *
*     Literature:                                                           *
*     Vogelezang DHP and Holtslag AAM (1996): Evaluation and model impacts  *
*     of alternative boundary-layer height formulations. Boundary-Layer     *
*     Meteor. 81, 245-269.                                                  *
*                                                                           *
*     Update: 1999-02-01 by G. Wotawa                                       *
*                                                                           *
*     Two meter level (temperature, humidity) is taken as reference level   *
*     instead of first model level.                                         *
*     New input variables tt2, td2 introduced.                              *
*                                                                           *
*     17 Oct 2005 - R. Easter - added ierr status flag (0=ok, -1=fail)      *
*     15 Nov 2005 - R. Easter - use pplev instead of akz,bkz                *
*                                                                           *
*****************************************************************************
*                                                                           *
* Variables:                                                                *
* h                          mixing height [m]                              *
* hf                         sensible heat flux                             *
* psurf                      surface pressure at point (xt,yt) [Pa]         *
* tv                         virtual temperature                            *
* wst                        convective velocity scale                      *
*                                                                           *
* Constants:                                                                *
* ric                        critical Richardson number                     *
*                                                                           *
*****************************************************************************

      include 'includepar'

      integer ierr
      integer i,k,nuvz,itmax,iter
      real tv,tvold,zref,z,zold,pint,pold,theta,thetaref,ri
c     real akz(nuvz),bkz(nuvz),ulev(nuvz),vlev(nuvz),hf,wst,tt2,td2,ew
      real         pplev(nuvz),ulev(nuvz),vlev(nuvz),hf,wst,tt2,td2,ew
      real psurf,ust,ttlev(nuvz),qvlev(nuvz),h,const,ric,b,excess,bs
      real thetaold,zl,ul,vl,thetal,ril,hmixplus,wspeed,bvfsq,bvf
      real f_qvsat,rh,rhold,rhl,theta1,theta2,zl1,zl2,thetam
      real duma
      integer sfc_option
      parameter(const=r_air/ga,ric=0.25,b=100.,bs=8.5,itmax=3)

      excess=0.0
      iter=0
 
C Compute virtual temperature and virtual potential temperature at
C reference level (2 m)
******************************************************************
 
30    iter=iter+1
 
      pold=psurf
      tvold=tt2*(1.+0.378*ew(td2)/psurf)
      zold=2.0
      zref=zold
      rhold=ew(td2)/ew(tt2)
 
      thetaref=tvold*(100000./pold)**(r_air/cpa)+excess
      thetaold=thetaref
 
 
C Integrate z up to one level above zt
**************************************
 
      do 10 k=2,nuvz
c       pint=akz(k)+bkz(k)*psurf  ! pressure on model layers
        pint=pplev(k)             ! pressure on model layers
        tv=ttlev(k)*(1.+0.608*qvlev(k))
 
        if (abs(tv-tvold).gt.0.2) then
          z=zold+const*log(pold/pint)*(tv-tvold)/log(tv/tvold)
        else
          z=zold+const*log(pold/pint)*tv
        endif
 
        theta=tv*(100000./pint)**(r_air/cpa)
! Petra
        rh = qvlev(k) / f_qvsat( pint, ttlev(k) )
 
 
Calculate Richardson number at each level
*****************************************
 
        ri=ga/thetaref*(theta-thetaref)*(z-zref)/
     +  max(((ulev(k)-ulev(2))**2+(vlev(k)-vlev(2))**2+b*ust**2),0.1)
 
C  addition of second condition: MH should not be placed in an
C  unstable layer (PS / Feb 2000)
        if (ri.gt.ric .and. thetaold.lt.theta) goto 20
 
        tvold=tv
        pold=pint
        rhold=rh
        thetaold=theta
10      zold=z

        if (k .ge. nuvz) then
            write(*,*) 'richardson not working -- k = nuvz'
            ierr = -10
            goto 7000
        end if
      
20    continue
 
C Determine Richardson number between the critical levels
*********************************************************

      zl1=zold
      theta1=thetaold
      do 15 i=1,20
        zl=zold+float(i)/20.*(z-zold)
        ul=ulev(k-1)+float(i)/20.*(ulev(k)-ulev(k-1))
        vl=vlev(k-1)+float(i)/20.*(vlev(k)-vlev(k-1))
        thetal=thetaold+float(i)/20.*(theta-thetaold)
        rhl=rhold+float(i)/20.*(rh-rhold)
        ril=ga/thetaref*(thetal-thetaref)*(zl-zref)/
     +  max(((ul-ulev(2))**2+(vl-vlev(2))**2+b*ust**2),0.1)
        zl2=zl
        theta2=thetal
        if (ril.gt.ric) goto 25
        zl1=zl
        theta1=thetal
15      continue
 
25    continue
c if sfc_option = sfc_option_wrf, 
c pbl heights are read from WRF met. files and put into hmix (=h)
c     h=zl
      if(sfc_option .eq. sfc_option_diagnosed) h=zl
      if (h .le. 0.0) then
          write(*,*) 'richardson not working -- bad h =', h
          ierr = -20
          goto 7000
c     else if (h .lt. 10.0) then
c         write(*,*) 'richardson not working -- too small h =', h
c         ierr = +20
c         return
      end if
      
      thetam=0.5*(theta1+theta2)
      wspeed=sqrt(ul**2+vl**2)                    ! Wind speed at z=hmix
      bvfsq=(ga/thetam)*(theta2-theta1)/(zl2-zl1) ! Brunt-Vaisala frequency
                                                  ! at z=hmix

C Under stable conditions, limit the maximum effect of the subgrid-scale topography
C by the maximum lifting possible from the available kinetic energy
***********************************************************************************

      if(bvfsq.le.0.) then
        hmixplus=9999.
      else
        bvf=sqrt(bvfsq)
        hmixplus=wspeed/bvf*convke                ! keconv = kinetic energy
      endif                                       ! used for lifting

 
C Calculate convective velocity scale
*************************************
 
      if (hf.lt.0.) then
        wst=(-h*ga/thetaref*hf/cpa)**0.333
        excess=-bs*hf/cpa/wst
        if (iter.lt.itmax) goto 30
      else
        wst=0.
      endif

      ierr = 0
      return

C Fatal error -- print the inputs
7000  continue
      write(*,'(a         )') 'nuvz'
      write(*,'(i5        )')  nuvz
      write(*,'(a         )') 'psurf,ust,hf,tt2,td2,h,wst,hmixplus'
      write(*,'(1p,4e18.10)')  psurf,ust,hf,tt2,td2,h,wst,hmixplus
      write(*,'(a         )') 'ttlev'
      write(*,'(1p,4e18.10)')  ttlev
      write(*,'(a         )') 'qvlev'
      write(*,'(1p,4e18.10)')  qvlev
      write(*,'(a         )') 'ulev'
      write(*,'(1p,4e18.10)')  ulev
      write(*,'(a         )') 'vlev'
      write(*,'(1p,4e18.10)')  vlev
      write(*,'(a         )') 'pplev'
      write(*,'(1p,4e18.10)')  pplev
      return

      end
