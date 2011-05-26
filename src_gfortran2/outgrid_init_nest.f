      subroutine outgrid_init_nest()
********************************************************************************
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subroutine outgrid_init_nest. *
*            The computational grid is the WRF x-y grid rather than lat-lon.   *
*                                                                              *
*  This routine calculates, for each grid cell of the output nest, the         *
*  volume and the surface area.                                                *
*                                                                              *
*     Author: A. Stohl                                                         *
*    30 August 2004                                                            *
*                                                                              *
*    26 Oct 2005, R. Easter - changes to gridarea                              *
*                             associated with WRF horizontal grid.             *
*    Dec 2005, R. Easter - changed names of "*lon0*" & "*lat0*" variables      *
*                                                                              *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
*                                                                              *
* arean              surface area of all output nest cells                     *
* volumen            volumes of all output nest cells                          *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'
      
      integer ix,jy,kz,k,nage,l,iix,jjy,ixp,jyp,i1,j1,j,ngrid
c     real ylat,gridarea,ylatp,ylatm,hzone,cosfact,cosfactm,cosfactp
      real ymet,gridarea
      real xmet,xl,yl,ddx,ddy,rddx,rddy,p1,p2,p3,p4,xtn,ytn,oroh


      write(*,*)
      write(*,*) '*** stopping in outgrid_init_nest ***'
      write(*,*) 
     &    '    the wrf version of this routine is not yet implemented'
      write(*,*)
      stop


C Compute surface area and volume of each grid cell: area, volume;
C and the areas of the northward and eastward facing walls: areaeast, areanorth
************************************************************************

      do 10 jy=0,numygridn-1

c        ylat=outlat0n+(float(jy)+0.5)*dyoutn
c        ylatp=ylat+0.5*dyoutn
c        ylatm=ylat-0.5*dyoutn
c        if ((ylatm.lt.0).and.(ylatp.gt.0.)) then
c          hzone=dyoutn*r_earth*pi180
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

c for FLEXPART_WRF, dx & dy are in meters, and no cos(lat) is needed
c        gridarea=2.*pi*r_earth*hzone*dxoutn/360.
        gridarea=dxoutn*dyoutn

        do 10 ix=0,numxgridn-1
          arean(ix,jy)=gridarea

C Volume = area x box height
****************************

          volumen(ix,jy,1)=arean(ix,jy)*outheight(1)
          do 10 kz=2,numzgrid
10        volumen(ix,jy,kz)=arean(ix,jy)*(outheight(kz)-outheight(kz-1))


***************************************************************************
C Determine average height of model topography in nesteed output grid cells
***************************************************************************

C Loop over all output grid cells
*********************************

      do 20 jjy=0,numygridn-1
        do 20 iix=0,numxgridn-1
          oroh=0.

C Take 100 samples of the topography in every grid cell
*******************************************************

          do 22 j1=1,10
c for FLEXPART_WRF, x & y coords are in meters,
c and the lon & lat variables below are in meters.
            ymet=out_ym0n+(float(jjy)+float(j1)/10.-0.05)*dyoutn
            yl=(ymet-ymet0)/dy
            do 22 i1=1,10
              xmet=out_xm0n+(float(iix)+float(i1)/10.-0.05)*dxoutn
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

          orooutn(iix,jjy)=oroh/100.
20        continue



********************************
C Initialization of output grids
********************************

      do 31 k=1,nspec
        do 31 nage=1,nageclass
          do 31 jy=0,numygridn-1
            do 31 ix=0,numxgridn-1
              do 31 l=1,nclassunc
C Deposition fields
                wetgriduncn(ix,jy,k,l,nage)=0.
                drygriduncn(ix,jy,k,l,nage)=0.
C Concentration fields
                do 31 kz=1,numzgrid
31                griduncn(ix,jy,kz,k,l,nage)=0.

 

      end
