      subroutine interpol_all_nests(itime,xt,yt,zt)
C                                     i   i  i  i
********************************************************************************
*                                                                              *
*  This subroutine interpolates everything that is needed for calculating the  *
*  dispersion.                                                                 *
*  Version for interpolating nested grids.                                     *
*                                                                              *
*    Author: A. Stohl                                                          *
*                                                                              *
*    9 February 1999                                                           *
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

        ust1(m)=p1*ustarn(ix ,jy ,1,indexh,ngrid)
     +        + p2*ustarn(ixp,jy ,1,indexh,ngrid)
     +        + p3*ustarn(ix ,jyp,1,indexh,ngrid)
     +        + p4*ustarn(ixp,jyp,1,indexh,ngrid)
        wst1(m)=p1*wstarn(ix ,jy ,1,indexh,ngrid)
     +        + p2*wstarn(ixp,jy ,1,indexh,ngrid)
     +        + p3*wstarn(ix ,jyp,1,indexh,ngrid)
     +        + p4*wstarn(ixp,jyp,1,indexh,ngrid)
80      oli1(m)=p1*olin(ix ,jy ,1,indexh,ngrid)
     +        + p2*olin(ixp,jy ,1,indexh,ngrid)
     +        + p3*olin(ix ,jyp,1,indexh,ngrid)
     +        + p4*olin(ixp,jyp,1,indexh,ngrid)

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

      do 10 n=indz,indz+1
        usl=0.
        vsl=0.
        wsl=0.
        usq=0.
        vsq=0.
        wsq=0.
        do 20 m=1,2
          indexh=memind(m)
          y1(m)=p1*uun(ix ,jy ,n,indexh,ngrid)
     +         +p2*uun(ixp,jy ,n,indexh,ngrid)
     +         +p3*uun(ix ,jyp,n,indexh,ngrid)
     +         +p4*uun(ixp,jyp,n,indexh,ngrid)
          y2(m)=p1*vvn(ix ,jy ,n,indexh,ngrid)
     +         +p2*vvn(ixp,jy ,n,indexh,ngrid)
     +         +p3*vvn(ix ,jyp,n,indexh,ngrid)
     +         +p4*vvn(ixp,jyp,n,indexh,ngrid)
          y3(m)=p1*wwn(ix ,jy ,n,indexh,ngrid)
     +         +p2*wwn(ixp,jy ,n,indexh,ngrid)
     +         +p3*wwn(ix ,jyp,n,indexh,ngrid)
     +         +p4*wwn(ixp,jyp,n,indexh,ngrid)
          rhograd1(m)=p1*drhodzn(ix ,jy ,n,indexh,ngrid)
     +               +p2*drhodzn(ixp,jy ,n,indexh,ngrid)
     +               +p3*drhodzn(ix ,jyp,n,indexh,ngrid)
     +               +p4*drhodzn(ixp,jyp,n,indexh,ngrid)
          rho1(m)=p1*rhon(ix ,jy ,n,indexh,ngrid)
     +           +p2*rhon(ixp,jy ,n,indexh,ngrid)
     +           +p3*rhon(ix ,jyp,n,indexh,ngrid)
     +           +p4*rhon(ixp,jyp,n,indexh,ngrid)

         usl=usl+uun(ix ,jy ,n,indexh,ngrid)+uun(ixp,jy ,n,indexh,ngrid)
     +    +uun(ix ,jyp,n,indexh,ngrid)+uun(ixp,jyp,n,indexh,ngrid)
         vsl=vsl+vvn(ix ,jy ,n,indexh,ngrid)+vvn(ixp,jy ,n,indexh,ngrid)
     +    +vvn(ix ,jyp,n,indexh,ngrid)+vvn(ixp,jyp,n,indexh,ngrid)
         wsl=wsl+wwn(ix ,jy ,n,indexh,ngrid)+wwn(ixp,jy ,n,indexh,ngrid)
     +    +wwn(ix ,jyp,n,indexh,ngrid)+wwn(ixp,jyp,n,indexh,ngrid)

        usq=usq+uun(ix ,jy ,n,indexh,ngrid)*uun(ix ,jy ,n,indexh,ngrid)+
     +    uun(ixp,jy ,n,indexh,ngrid)*uun(ixp,jy ,n,indexh,ngrid)+
     +    uun(ix ,jyp,n,indexh,ngrid)*uun(ix ,jyp,n,indexh,ngrid)+
     +    uun(ixp,jyp,n,indexh,ngrid)*uun(ixp,jyp,n,indexh,ngrid)
        vsq=vsq+vvn(ix ,jy ,n,indexh,ngrid)*vvn(ix ,jy ,n,indexh,ngrid)+
     +    vvn(ixp,jy ,n,indexh,ngrid)*vvn(ixp,jy ,n,indexh,ngrid)+
     +    vvn(ix ,jyp,n,indexh,ngrid)*vvn(ix ,jyp,n,indexh,ngrid)+
     +    vvn(ixp,jyp,n,indexh,ngrid)*vvn(ixp,jyp,n,indexh,ngrid)
        wsq=wsq+wwn(ix ,jy ,n,indexh,ngrid)*wwn(ix ,jy ,n,indexh,ngrid)+
     +    wwn(ixp,jy ,n,indexh,ngrid)*wwn(ixp,jy ,n,indexh,ngrid)+
     +    wwn(ix ,jyp,n,indexh,ngrid)*wwn(ix ,jyp,n,indexh,ngrid)+
     +    wwn(ixp,jyp,n,indexh,ngrid)*wwn(ixp,jyp,n,indexh,ngrid)
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
