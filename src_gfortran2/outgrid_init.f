      subroutine outgrid_init()
********************************************************************************
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subroutine outgrid_init.      *
*            The computational grid is the WRF x-y grid rather than lat-lon.   *
*                                                                              *
*  This routine calculates, for each grid cell of the output grid, the         *
*  volume, the surface area, and the areas of the northward and eastward       *
*  facing surfaces.                                                            *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     7 August 2002                                                            *
*                                                                              *
*    26 Oct 2005, R. Easter - changes in gridarea, areaeast, areanorth         *
*                             associated with WRF horizontal grid.             *
*     Dec 2005, R. Easter - changed names of "*lon0*" & "*lat0*" variables     *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
*                                                                              *
* area               surface area of all output grid cells                     *
* areaeast           eastward facing wall area of all output grid cells        *
* areanorth          northward facing wall area of all output grid cells       *
* volume             volumes of all output grid cells                          *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'
      
      integer ix,jy,kz,k,i,nage,l,iix,jjy,ixp,jyp,i1,j1,j,ngrid
c     real ylat,gridarea,ylatp,ylatm,hzone,cosfact,cosfactm,cosfactp
      real ymet,gridarea
      real xmet,xl,yl,ddx,ddy,rddx,rddy,p1,p2,p3,p4,xtn,ytn,oroh



C Compute surface area and volume of each grid cell: area, volume;
C and the areas of the northward and eastward facing walls: areaeast, areanorth
************************************************************************

      do 10 jy=0,numygrid-1

c        ylat=outlat0+(float(jy)+0.5)*dyout
c        ylatp=ylat+0.5*dyout
c        ylatm=ylat-0.5*dyout
c        if ((ylatm.lt.0).and.(ylatp.gt.0.)) then
c          hzone=dyout*r_earth*pi180
c        else
c
cC Calculate area of grid cell with formula M=2*pi*R*h*dx/360,
cC see Netz, Formeln der Mathematik, 5. Auflage (1983), p.90
c*************************************************************
c
c          cosfact=cos(ylat*pi180)*r_earth
c          cosfactp=cos(ylatp*pi180)*r_earth
c          cosfactm=cos(ylatm*pi180)*r_earth
c          if (cosfactp.lt.cosfactm) then
c            hzone=sqrt(r_earth**2-cosfactp**2)-
c     +      sqrt(r_earth**2-cosfactm**2)
c          else
c            hzone=sqrt(r_earth**2-cosfactm**2)-
c     +      sqrt(r_earth**2-cosfactp**2)
c          endif
c        endif
c
cC Surface are of a grid cell at a latitude ylat
c***********************************************
c
c        gridarea=2.*pi*r_earth*hzone*dxout/360.

c for FLEXPART_WRF, dx & dy are in meters, and no cos(lat) is needed
c ??? maybe should incorporate map factor here,
c     and also for areaeast & areanorth ???
        gridarea=dxout*dyout

        do 10 ix=0,numxgrid-1
          area(ix,jy)=gridarea

C Volume = area x box height
****************************

          volume(ix,jy,1)=area(ix,jy)*outheight(1)

c for FLEXPART_WRF, dx & dy are in meters, and no cos(lat) is needed
c          areaeast(ix,jy,1)=dyout*r_earth*pi180*outheight(1)
c          areanorth(ix,jy,1)=cos(ylat*pi180)*dxout*r_earth*pi180*
c     +    outheight(1)
          areaeast(ix,jy,1)=dyout*outheight(1)
          areanorth(ix,jy,1)=dxout*outheight(1)

          do 10 kz=2,numzgrid

c            areaeast(ix,jy,kz)=dyout*r_earth*pi180*
c     +      (outheight(kz)-outheight(kz-1))
c            areanorth(ix,jy,kz)=cos(ylat*pi180)*dxout*r_earth*pi180*
c     +      (outheight(kz)-outheight(kz-1))
            areaeast(ix,jy,kz)=dyout*(outheight(kz)-outheight(kz-1))
            areanorth(ix,jy,kz)=dxout*(outheight(kz)-outheight(kz-1))

10          volume(ix,jy,kz)=area(ix,jy)*(outheight(kz)-outheight(kz-1))




*******************************************************************
C Determine average height of model topography in output grid cells
*******************************************************************

C Loop over all output grid cells
*********************************

      do 20 jjy=0,numygrid-1
        do 20 iix=0,numxgrid-1
          oroh=0.

C Take 100 samples of the topography in every grid cell
*******************************************************

          do 22 j1=1,10
c for FLEXPART_WRF, x & y coords are in meters,
c and the lon & lat variables below are in meters.
            ymet=out_ym0+(float(jjy)+float(j1)/10.-0.05)*dyout
            yl=(ymet-ymet0)/dy
            do 22 i1=1,10
              xmet=out_xm0+(float(iix)+float(i1)/10.-0.05)*dxout
              xl=(xmet-xmet0)/dx

C Determine the nest we are in
******************************

              ngrid=0
              do 42 j=numbnests,1,-1
                if ((xl.gt.xln(j)).and.(xl.lt.xrn(j)).and.
     +          (yl.gt.yln(j)).and.(yl.lt.yrn(j))) then
                  ngrid=j
                  goto 43
                endif
42              continue
43            continue

C Determine (nested) grid coordinates and auxiliary parameters used for interpolation
*************************************************************************************

              if (ngrid.gt.0) then
                xtn=(xl-xln(ngrid))*xresoln(ngrid)
                ytn=(yl-yln(ngrid))*yresoln(ngrid)
                ix=int(xtn)
                jy=int(ytn)
                ddy=ytn-float(jy)
                ddx=xtn-float(ix)
              else
                ix=int(xl)
                jy=int(yl)
                ddy=yl-float(jy)
                ddx=xl-float(ix)
              endif
              ixp=ix+1
              jyp=jy+1
              rddx=1.-ddx
              rddy=1.-ddy
              p1=rddx*rddy
              p2=ddx*rddy
              p3=rddx*ddy
              p4=ddx*ddy

              if (ngrid.gt.0) then
                oroh=oroh+p1*oron(ix ,jy ,ngrid)
     +                  + p2*oron(ixp,jy ,ngrid)
     +                  + p3*oron(ix ,jyp,ngrid)
     +                  + p4*oron(ixp,jyp,ngrid)
              else
                oroh=oroh+p1*oro(ix ,jy)
     +                  + p2*oro(ixp,jy)
     +                  + p3*oro(ix ,jyp)
     +                  + p4*oro(ixp,jyp)
              endif
22            continue

C Divide by the number of samples taken
***************************************

          oroout(iix,jjy)=oroh/100.
20        continue



*************************
C Initialize output grids
*************************

      do 31 k=1,nspec
        do 33 i=1,numreceptor
C Receptor points
33        creceptor(i,k)=0.
        do 31 nage=1,nageclass
          do 31 jy=0,numygrid-1
            do 31 ix=0,numxgrid-1
              do 31 l=1,nclassunc
C Deposition fields
                wetgridunc(ix,jy,k,l,nage)=0.
                drygridunc(ix,jy,k,l,nage)=0.
                do 31 kz=1,numzgrid
C Flux fields
                  fluxe(ix,jy,kz,k,nage)=0.
                  fluxw(ix,jy,kz,k,nage)=0.
                  fluxn(ix,jy,kz,k,nage)=0.
                  fluxs(ix,jy,kz,k,nage)=0.
                  fluxu(ix,jy,kz,k,nage)=0.
                  fluxd(ix,jy,kz,k,nage)=0.
C Concentration fields
31                gridunc(ix,jy,kz,k,l,nage)=0.

 

      end
