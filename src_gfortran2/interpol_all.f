      subroutine interpol_all(itime,xt,yt,zt)
C                               i   i  i  i
********************************************************************************
*                                                                              *
*  This subroutine interpolates everything that is needed for calculating the  *
*  dispersion.                                                                 *
*                                                                              *
*    Author: A. Stohl                                                          *
*                                                                              *
*    16 December 1997                                                          *
*                                                                              *
*  Revision March 2005 by AST : all output variables in common block           *
*                               calculation of standard deviation done in this *
*                               routine rather than subroutine call in order   *
*                               to save computation time                       *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* itime [s]          current temporal position                                 *
* memtime(3) [s]     times of the wind fields in memory                        *
* xt,yt,zt           coordinates position for which wind data shall be calculat*
*                                                                              *
* Constants:                                                                   *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'
      include 'includeinterpol'
      include 'includehanna'

      integer itime
      real xt,yt,zt

C Auxiliary variables needed for interpolation
      real ust1(2),wst1(2),oli1(2),oliaux
      real y1(2),y2(2),y3(2),rho1(2),rhograd1(2)
      real usl,vsl,wsl,usq,vsq,wsq,xaux,eps
      integer i,m,n,indexh
      parameter(eps=1.0e-30)


*********************************************
C Multilinear interpolation in time and space
*********************************************

C Determine the lower left corner and its distance to the current position
**************************************************************************

      ddx=xt-float(ix)                     
      ddy=yt-float(jy)
      rddx=1.-ddx
      rddy=1.-ddy
      p1=rddx*rddy
      p2=ddx*rddy
      p3=rddx*ddy
      p4=ddx*ddy

C Calculate variables for time interpolation
********************************************

      dt1=float(itime-memtime(1))
      dt2=float(memtime(2)-itime)
      dtt=1./(dt1+dt2)


******************************************
C 1. Interpolate u*, w* and Obukhov length
******************************************

C a) Bilinear horizontal interpolation

      do 80 m=1,2
        indexh=memind(m)

        ust1(m)=p1*ustar(ix ,jy ,1,indexh)
     +        + p2*ustar(ixp,jy ,1,indexh)
     +        + p3*ustar(ix ,jyp,1,indexh)
     +        + p4*ustar(ixp,jyp,1,indexh)
        wst1(m)=p1*wstar(ix ,jy ,1,indexh)
     +        + p2*wstar(ixp,jy ,1,indexh)
     +        + p3*wstar(ix ,jyp,1,indexh)
     +        + p4*wstar(ixp,jyp,1,indexh)
80      oli1(m)=p1*oli(ix ,jy ,1,indexh)
     +        + p2*oli(ixp,jy ,1,indexh)
     +        + p3*oli(ix ,jyp,1,indexh)
     +        + p4*oli(ixp,jyp,1,indexh)

C b) Temporal interpolation

      ust=(ust1(1)*dt2+ust1(2)*dt1)*dtt
      wst=(wst1(1)*dt2+wst1(2)*dt1)*dtt
      oliaux=(oli1(1)*dt2+oli1(2)*dt1)*dtt

      if (oliaux.ne.0.) then
        ol=1./oliaux
      else
        ol=99999.
      endif


******************************************************
C 2. Interpolate vertical profiles of u,v,w,rho,drhodz
******************************************************


C Determine the level below the current position
************************************************

      do 5 i=2,nz
        if (height(i).gt.zt) then
          indz=i-1
          indzp=i
          goto 6
        endif
5       continue
6     continue

***************************************
C 1.) Bilinear horizontal interpolation
C 2.) Temporal interpolation (linear)
***************************************

C Loop over 2 time steps and indz levels
****************************************

      do 10 n=indz,indzp
        usl=0.
        vsl=0.
        wsl=0.
        usq=0.
        vsq=0.
        wsq=0.
        do 20 m=1,2
          indexh=memind(m)
          if (ngrid.lt.0) then
            y1(m)=p1*uupol(ix ,jy ,n,indexh)
     +           +p2*uupol(ixp,jy ,n,indexh)
     +           +p3*uupol(ix ,jyp,n,indexh)
     +           +p4*uupol(ixp,jyp,n,indexh)
            y2(m)=p1*vvpol(ix ,jy ,n,indexh)
     +           +p2*vvpol(ixp,jy ,n,indexh)
     +           +p3*vvpol(ix ,jyp,n,indexh)
     +           +p4*vvpol(ixp,jyp,n,indexh)
            usl=usl+uupol(ix ,jy ,n,indexh)+uupol(ixp,jy ,n,indexh)
     +      +uupol(ix ,jyp,n,indexh)+uupol(ixp,jyp,n,indexh)
            vsl=vsl+vvpol(ix ,jy ,n,indexh)+vvpol(ixp,jy ,n,indexh)
     +      +vvpol(ix ,jyp,n,indexh)+vvpol(ixp,jyp,n,indexh)

            usq=usq+uupol(ix ,jy ,n,indexh)*uupol(ix ,jy ,n,indexh)+
     +      uupol(ixp,jy ,n,indexh)*uupol(ixp,jy ,n,indexh)+
     +      uupol(ix ,jyp,n,indexh)*uupol(ix ,jyp,n,indexh)+
     +      uupol(ixp,jyp,n,indexh)*uupol(ixp,jyp,n,indexh)
            vsq=vsq+vvpol(ix ,jy ,n,indexh)*vvpol(ix ,jy ,n,indexh)+
     +      vvpol(ixp,jy ,n,indexh)*vvpol(ixp,jy ,n,indexh)+
     +      vvpol(ix ,jyp,n,indexh)*vvpol(ix ,jyp,n,indexh)+
     +      vvpol(ixp,jyp,n,indexh)*vvpol(ixp,jyp,n,indexh)
          else
            y1(m)=p1*uu(ix ,jy ,n,indexh)
     +           +p2*uu(ixp,jy ,n,indexh)
     +           +p3*uu(ix ,jyp,n,indexh)
     +           +p4*uu(ixp,jyp,n,indexh)
            y2(m)=p1*vv(ix ,jy ,n,indexh)
     +           +p2*vv(ixp,jy ,n,indexh)
     +           +p3*vv(ix ,jyp,n,indexh)
     +           +p4*vv(ixp,jyp,n,indexh)
            usl=usl+uu(ix ,jy ,n,indexh)+uu(ixp,jy ,n,indexh)
     +      +uu(ix ,jyp,n,indexh)+uu(ixp,jyp,n,indexh)
            vsl=vsl+vv(ix ,jy ,n,indexh)+vv(ixp,jy ,n,indexh)
     +      +vv(ix ,jyp,n,indexh)+vv(ixp,jyp,n,indexh)

            usq=usq+uu(ix ,jy ,n,indexh)*uu(ix ,jy ,n,indexh)+
     +      uu(ixp,jy ,n,indexh)*uu(ixp,jy ,n,indexh)+
     +      uu(ix ,jyp,n,indexh)*uu(ix ,jyp,n,indexh)+
     +      uu(ixp,jyp,n,indexh)*uu(ixp,jyp,n,indexh)
            vsq=vsq+vv(ix ,jy ,n,indexh)*vv(ix ,jy ,n,indexh)+
     +      vv(ixp,jy ,n,indexh)*vv(ixp,jy ,n,indexh)+
     +      vv(ix ,jyp,n,indexh)*vv(ix ,jyp,n,indexh)+
     +      vv(ixp,jyp,n,indexh)*vv(ixp,jyp,n,indexh)
          endif
          y3(m)=p1*ww(ix ,jy ,n,indexh)
     +         +p2*ww(ixp,jy ,n,indexh)
     +         +p3*ww(ix ,jyp,n,indexh)
     +         +p4*ww(ixp,jyp,n,indexh)
          rhograd1(m)=p1*drhodz(ix ,jy ,n,indexh)
     +               +p2*drhodz(ixp,jy ,n,indexh)
     +               +p3*drhodz(ix ,jyp,n,indexh)
     +               +p4*drhodz(ixp,jyp,n,indexh)
          rho1(m)=p1*rho(ix ,jy ,n,indexh)
     +           +p2*rho(ixp,jy ,n,indexh)
     +           +p3*rho(ix ,jyp,n,indexh)
     +           +p4*rho(ixp,jyp,n,indexh)
          wsl=wsl+ww(ix ,jy ,n,indexh)+ww(ixp,jy ,n,indexh)
     +    +ww(ix ,jyp,n,indexh)+ww(ixp,jyp,n,indexh)
          wsq=wsq+ww(ix ,jy ,n,indexh)*ww(ix ,jy ,n,indexh)+
     +    ww(ixp,jy ,n,indexh)*ww(ixp,jy ,n,indexh)+
     +    ww(ix ,jyp,n,indexh)*ww(ix ,jyp,n,indexh)+
     +    ww(ixp,jyp,n,indexh)*ww(ixp,jyp,n,indexh)
20        continue
        uprof(n)=(y1(1)*dt2+y1(2)*dt1)*dtt
        vprof(n)=(y2(1)*dt2+y2(2)*dt1)*dtt
        wprof(n)=(y3(1)*dt2+y3(2)*dt1)*dtt
        rhoprof(n)=(rho1(1)*dt2+rho1(2)*dt1)*dtt
        rhogradprof(n)=(rhograd1(1)*dt2+rhograd1(2)*dt1)*dtt
        indzindicator(n)=.false.

C Compute standard deviations
*****************************

        xaux=usq-usl*usl/8.
        if (xaux.lt.eps) then
          usigprof(n)=0.
        else
          usigprof(n)=sqrt(xaux/7.)
        endif

        xaux=vsq-vsl*vsl/8.
        if (xaux.lt.eps) then
          vsigprof(n)=0.
        else
          vsigprof(n)=sqrt(xaux/7.)
        endif


        xaux=wsq-wsl*wsl/8.
        if (xaux.lt.eps) then
          wsigprof(n)=0.
        else
          wsigprof(n)=sqrt(xaux/7.)
        endif

10      continue


      end
