      subroutine partoutput(itime)
C                             i
********************************************************************************
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subroutine partoutput.        *
*                                                                              *
*     Dump all particle positions                                              *
*                                                                              *
*     Author: A. Stohl                                                         *
*     12 March 1999                                                            *
*                                                                              *
*     Dec 2005, J. Fast - Output files can be either binary or ascii.          *
*                   Topo,pv,qv,... at particle positions are calculated        *
*                   using nested fields when (partoutput_use_nested .gt. 0)    *
*                   Particle xy coords can be either lat-lon or grid-meters.   *
*                   Changed names of "*lon0*" & "*lat0*" variables             *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      double precision jul
      integer itime,i,j,jjjjmmdd,ihmmss
      integer ix,jy,ixp,jyp,indexh,m,il,ind,indz,indzp
      integer numpart_out
      real xlon,ylat,xtmp,ytmp
      real dt1,dt2,dtt,ddx,ddy,rddx,rddy,p1,p2,p3,p4,dz1,dz2,dz
      real topo,hm(2),hmixi,pv1(2),pvprof(2),pvi,qv1(2),qvprof(2),qvi
      real tt1(2),ttprof(2),tti,rho1(2),rhoprof(2),rhoi
      real tr(2),tri
      character adate*8,atime*6

      integer k, ngrid
      real xtn, ytn


C Determine current calendar date, needed for the file name
***********************************************************

      jul=bdate+dble(float(itime))/86400.
      call caldate(jul,jjjjmmdd,ihmmss)
      write(adate,'(i8.8)') jjjjmmdd
      write(atime,'(i6.6)') ihmmss


C Some variables needed for temporal interpolation
**************************************************

      dt1=float(itime-memtime(1))
      dt2=float(memtime(2)-itime)
      dtt=1./(dt1+dt2)

C Open output file and write the output
***************************************

      if (ipout.eq.1) then
        if (iouttype.eq.0)
     +  open(unitpartout,file=path(2)(1:len(2))//'partposit_'//adate//
     +  atime,form='unformatted')
        if (iouttype.eq.1)
     +  open(unitpartout,file=path(2)(1:len(2))//'partposit_'//adate//
     +  atime,form='formatted')
      else
        if (iouttype.eq.0)
     +  open(unitpartout,file=path(2)(1:len(2))//'partposit_end',
     +  form='unformatted')
        if (iouttype.eq.1)
     +  open(unitpartout,file=path(2)(1:len(2))//'partposit_end',
     +  form='formatted')
      endif

C Write current time to file
C FLEXPART_WRF - write iomode_xycoord
****************************

      numpart_out = 0
      do i=1,numpart
        if (itra1(i).eq.itime) numpart_out = numpart_out + 1
      enddo

      if (iouttype.eq.0) write(unitpartout)   itime,
     +    numpart_out, iomode_xycoord
      if (iouttype.eq.1) write(unitpartout,*) itime,
     +    numpart_out, iomode_xycoord

      do 100 i=1,numpart

C Take only valid particles
***************************

        if (itra1(i).eq.itime) then
          xlon=xmet0+xtra1(i)*dx
          ylat=ymet0+ytra1(i)*dy

**********************************************************************************
C Interpolate several variables (PV, specific humidity, etc.) to particle position
**********************************************************************************

C If partoutput_use_nested=0, set ngrid=0, and use the outermost grid
C    for calculating topo, pv, qv, ... at the particle position
C Otherwise, determine the nest we are in
          ngrid=0
          if (partoutput_use_nested .gt. 0) then
          do 25 k=numbnests,1,-1
             if ((xtra1(i).gt.xln(k)).and.
     +           (xtra1(i).lt.xrn(k)).and.
     +           (ytra1(i).gt.yln(k)).and.
     +           (ytra1(i).lt.yrn(k))) then
                    ngrid=k
                    goto 26
             endif
25        continue
26        continue
          endif

          if (ngrid .le. 0) then
             ix=int(xtra1(i))
             jy=int(ytra1(i))
             ddy=ytra1(i)-float(jy)
             ddx=xtra1(i)-float(ix)
          else
             xtn=(xtra1(i)-xln(ngrid))*xresoln(ngrid)
             ytn=(ytra1(i)-yln(ngrid))*yresoln(ngrid)
             ix=int(xtn)
             jy=int(ytn)
             ddy=ytn-float(jy)
             ddx=xtn-float(ix)
          endif
          ixp=ix+1
          jyp=jy+1
          rddx=1.-ddx
          rddy=1.-ddy
          p1=rddx*rddy
          p2=ddx*rddy
          p3=rddx*ddy
          p4=ddx*ddy

C Topography
************
          if (ngrid .le. 0) then
             topo=p1*oro(ix ,jy)
     +          + p2*oro(ixp,jy)
     +          + p3*oro(ix ,jyp)
     +          + p4*oro(ixp,jyp)
          else
             topo=p1*oron(ix ,jy ,ngrid)
     +          + p2*oron(ixp,jy ,ngrid)
     +          + p3*oron(ix ,jyp,ngrid)
     +          + p4*oron(ixp,jyp,ngrid)
          endif

C Potential vorticity, specific humidity, temperature, and density
******************************************************************

          do 55 il=2,nz
            if (height(il).gt.ztra1(i)) then
              indz=il-1
              indzp=il
              goto 56
            endif
55          continue
56        continue

          dz1=ztra1(i)-height(indz)
          dz2=height(indzp)-ztra1(i)
          dz=1./(dz1+dz2)


          do 65 ind=indz,indzp
            do 60 m=1,2
              indexh=memind(m)

              if (ngrid .le. 0) then
C Potential vorticity
              pv1(m)=p1*pv(ix ,jy ,ind,indexh)
     +              +p2*pv(ixp,jy ,ind,indexh)
     +              +p3*pv(ix ,jyp,ind,indexh)
     +              +p4*pv(ixp,jyp,ind,indexh)
C Specific humidity
              qv1(m)=p1*qv(ix ,jy ,ind,indexh)
     +              +p2*qv(ixp,jy ,ind,indexh)
     +              +p3*qv(ix ,jyp,ind,indexh)
     +              +p4*qv(ixp,jyp,ind,indexh)
C Temperature
              tt1(m)=p1*tt(ix ,jy ,ind,indexh)
     +              +p2*tt(ixp,jy ,ind,indexh)
     +              +p3*tt(ix ,jyp,ind,indexh)
     +              +p4*tt(ixp,jyp,ind,indexh)
C Density
              rho1(m)=p1*rho(ix ,jy ,ind,indexh)
     +               +p2*rho(ixp,jy ,ind,indexh)
     +               +p3*rho(ix ,jyp,ind,indexh)
     +               +p4*rho(ixp,jyp,ind,indexh)

              else
              pv1(m)=p1*pvn(ix ,jy ,ind,indexh,ngrid)
     +              +p2*pvn(ixp,jy ,ind,indexh,ngrid)
     +              +p3*pvn(ix ,jyp,ind,indexh,ngrid)
     +              +p4*pvn(ixp,jyp,ind,indexh,ngrid)
              qv1(m)=p1*qvn(ix ,jy ,ind,indexh,ngrid)
     +              +p2*qvn(ixp,jy ,ind,indexh,ngrid)
     +              +p3*qvn(ix ,jyp,ind,indexh,ngrid)
     +              +p4*qvn(ixp,jyp,ind,indexh,ngrid)
              tt1(m)=p1*ttn(ix ,jy ,ind,indexh,ngrid)
     +              +p2*ttn(ixp,jy ,ind,indexh,ngrid)
     +              +p3*ttn(ix ,jyp,ind,indexh,ngrid)
     +              +p4*ttn(ixp,jyp,ind,indexh,ngrid)
              rho1(m)=p1*rhon(ix ,jy ,ind,indexh,ngrid)
     +               +p2*rhon(ixp,jy ,ind,indexh,ngrid)
     +               +p3*rhon(ix ,jyp,ind,indexh,ngrid)
     +               +p4*rhon(ixp,jyp,ind,indexh,ngrid)

              endif

60            continue
            pvprof(ind-indz+1)=(pv1(1)*dt2+pv1(2)*dt1)*dtt
            qvprof(ind-indz+1)=(qv1(1)*dt2+qv1(2)*dt1)*dtt
            ttprof(ind-indz+1)=(tt1(1)*dt2+tt1(2)*dt1)*dtt
            rhoprof(ind-indz+1)=(rho1(1)*dt2+rho1(2)*dt1)*dtt
65          continue
          pvi=(dz1*pvprof(2)+dz2*pvprof(1))*dz
          qvi=(dz1*qvprof(2)+dz2*qvprof(1))*dz
          tti=(dz1*ttprof(2)+dz2*ttprof(1))*dz
          rhoi=(dz1*rhoprof(2)+dz2*rhoprof(1))*dz

C Tropopause and PBL height
***************************

          do 75 m=1,2
            indexh=memind(m)

            if (ngrid .le. 0) then
C Tropopause
            tr(m)=p1*tropopause(ix ,jy ,1,indexh)
     +          + p2*tropopause(ixp,jy ,1,indexh)
     +          + p3*tropopause(ix ,jyp,1,indexh)
     +          + p4*tropopause(ixp,jyp,1,indexh)
C PBL height
            hm(m)=p1*hmix(ix ,jy ,1,indexh)
     +          + p2*hmix(ixp,jy ,1,indexh)
     +          + p3*hmix(ix ,jyp,1,indexh)
     +          + p4*hmix(ixp,jyp,1,indexh)

            else
            tr(m)=p1*tropopausen(ix ,jy ,1,indexh,ngrid)
     +          + p2*tropopausen(ixp,jy ,1,indexh,ngrid)
     +          + p3*tropopausen(ix ,jyp,1,indexh,ngrid)
     +          + p4*tropopausen(ixp,jyp,1,indexh,ngrid)
            hm(m)=p1*hmixn(ix ,jy ,1,indexh,ngrid)
     +          + p2*hmixn(ixp,jy ,1,indexh,ngrid)
     +          + p3*hmixn(ix ,jyp,1,indexh,ngrid)
     +          + p4*hmixn(ixp,jyp,1,indexh,ngrid)

            endif

75        continue

          hmixi=(hm(1)*dt2+hm(2)*dt1)*dtt
          tri=(tr(1)*dt2+tr(2)*dt1)*dtt


C Write the output
******************

          if (iomode_xycoord .eq. iomode_xycoord_latlon) then
             xtmp = xlon
             ytmp = ylat
             call xymeter_to_ll_wrf( xtmp, ytmp, xlon, ylat )
          endif
          if(iouttype.eq.0)    
     +    write(unitpartout) npoint(i),xlon,ylat,ztra1(i),
     +    itramem(i),topo,pvi,qvi,rhoi,hmixi,tri,tti,
     +    (xmass1(i,j),j=1,nspec)
          if(iouttype.eq.1)    
     +    write(unitpartout,101) npoint(i),itramem(i),xlon,ylat,
     +    ztra1(i),topo,pvi,qvi,rhoi,hmixi,tri,tti,
     +    (xmass1(i,j),j=1,nspec)
        endif

100      continue

      if(iouttype.eq.0) 
     +write(unitpartout) -99999,-9999.9,-9999.9,-9999.9,-99999,
     +-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,
     +(-9999.9,j=1,nspec)
      if(iouttype.eq.1) 
     +write(unitpartout,101) -99999,-99999,-9999.9,-9999.9,-9999.9,
     +-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,-9999.9,
     +(-9999.9,j=1,nspec)

101   format( 2i10, 1p, 21e14.6 )


      close(unitpartout)

      end
