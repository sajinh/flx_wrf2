      subroutine assignland()
********************************************************************************
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subroutine assignland.        *
*            The computational grid is the WRF x-y grid rather than lat-lon.   *
*                                                                              *
*     This routine assigns fractions of the 8 landuse classes to each ECMWF    *
*     grid point.                                                              *
*     The landuse inventory of                                                 *
*                                                                              *
*     van de Velde R.J., Faber W.S., van Katwijk V.F., Kuylenstierna J.C.I.,   *
*     Scholten H.J., Thewessen T.J.M., Verspuij M., Zevenbergen M. (1994):     *
*     The Preparation of a European Land Use Database. National Institute of   *
*     Public Health and Environmental Protection, Report nr 712401001,         *
*     Bilthoven, The Netherlands.                                              *
*                                                                              *
*     is used to create a detailed landuse inventory for Europe.               *
*                                                                              *
*     Outside of Europe, the ECMWF land/sea mask is used to distinguish        *
*     between sea (-> ocean) and land (-> grasslands).                         *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     5 December 1996                                                          *
*     8 February 1999 Additional use of nests, A. Stohl                        *
*                                                                              *
*    14 October  2005 R. Easter -- modified for WRF.                           *
*                     The landuse inventory is not used at all.                *
*                     The land/sea mask is used everywhere.                    *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* xlanduse          fractions of numclass landuses for each model grid point   *
* xlandinvent       landuse inventory (1/6 deg resolution)                     *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      integer ix,jy,i,j,k,n,l
      real x,y,xlon,ylat

      do 10 ix=0,nxmin1
c       xlon=float(ix)*dx+xlon0         ! long.
        do 10 jy=0,nymin1
c         ylat=float(jy)*dy+ylat0       ! and lat. of each gridpoint
c FLEXPART_WRF - use this routine to get lat,lon
          call xyindex_to_ll_wrf( 0, float(ix), float(jy), xlon, ylat )
          do 20 k=1,numclass
20          xlanduse(ix,jy,k)=0.

          n=0

C Sample landuse inventory around grid point with 1/6 degree resolution
***********************************************************************
c          do 30 x=xlon-dx/2.+0.083,xlon+dx/2.-0.083,1./6.
c            do 30 y=ylat-dy/2.+0.083,ylat+dy/2.-0.083,1./6.
c              if ((x.gt.-24.66666).and.(x.lt.67.66666).and.
c     +        (y.gt.34.83333).and.(y.lt.72.0)) then           ! in domain
c                n=n+1
c                i=int((x+24.66666)*6.)+1
c                j=int((y-34.83333)*6.)+1
c                do 40 k=1,numclass
c40                xlanduse(ix,jy,k)=xlanduse(ix,jy,k)+xlandinvent(i,j,k)
c              endif
c30            continue

          if (n.ge.1) then                     ! detailed landuse available
            do 50 k=1,numclass
50            xlanduse(ix,jy,k)=xlanduse(ix,jy,k)/float(n)
          else                                    ! check land/sea mask
            if (lsm(ix,jy).lt.0.1) then           ! over sea  -> ocean
              xlanduse(ix,jy,9)=1.
            else                                  ! over land -> grasslands
              xlanduse(ix,jy,1)=1.
            endif
          endif
10        continue


*****************************************
C Same as above, but for the nested grids
*****************************************


      do 60 l=1,numbnests

        do 60 ix=0,nxn(l)-1
c         xlon=float(ix)*dxn(l)+xlon0n(l)         ! long. and
          do 60 jy=0,nyn(l)-1
c           ylat=float(jy)*dyn(l)+ylat0n(l)       ! lat. of each gridpoint
c FLEXPART_WRF - use this routine to get lat,lon
            call xyindex_to_ll_wrf( l, float(ix), float(jy), 
     &                              xlon, ylat )
            do 70 k=1,numclass
70            xlandusen(ix,jy,k,l)=0.
            n=0

C Sample landuse inventory around grid point with 1/6 degree resolution
***********************************************************************
c            do 80 x=xlon-dx/2.+0.083,xlon+dx/2.-0.083,1./6.
c              do 80 y=ylat-dy/2.+0.083,ylat+dy/2.-0.083,1./6.
c                if ((x.gt.-24.66666).and.(x.lt.67.66666).and.
c     +          (y.gt.34.83333).and.(y.lt.72.0)) then           ! in domain
c                  n=n+1
c                  i=int((x+24.66666)*6.)+1
c                  j=int((y-34.83333)*6.)+1
c                  do 90 k=1,numclass
c90                  xlandusen(ix,jy,k,l)=xlandusen(ix,jy,k,l)+
c     +              xlandinvent(i,j,k)
c                endif
c80              continue

            if (n.ge.1) then                     ! detailed landuse available
              do 100 k=1,numclass
100             xlandusen(ix,jy,k,l)=xlandusen(ix,jy,k,l)/float(n)
            else                                    ! check land/sea mask
              if (lsmn(ix,jy,l).lt.0.1) then   ! over sea  -> ocean
                xlandusen(ix,jy,9,l)=1.
              else                        ! over land -> grasslands
                xlandusen(ix,jy,1,l)=1.
              endif
            endif
60          continue

      end
