      subroutine clustering(xl,yl,zl,n,xclust,yclust,zclust,fclust,rms,
     +rmsclust,zrms)
C                           i  i  i  i   o      o      o      o     o
C        o      o
********************************************************************************
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subroutine clustering.        *
*            The computational grid is the WRF x-y grid rather than lat-lon.   *
*                                                                              *
*   This routine clusters the particle position into ncluster custers.         *
*   Input are the longitudes (xl) and latitudes (yl) of the individual         *
*   points, output are the cluster mean positions (xclust,yclust).             *
*   Vertical positions are not directly used for the clustering.               *
*                                                                              *
*   For clustering, the procedure described in Dorling et al. (1992) is used.  *
*                                                                              *
*   Dorling, S.R., Davies, T.D. and Pierce, C.E. (1992):                       *
*   Cluster analysis: a technique for estimating the synoptic meteorological   *
*   controls on air and precipitation chemistry - method and applications.     *
*   Atmospheric Environment 26A, 2575-2581.                                    *
*                                                                              *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     1 February 2002                                                          *
*                                                                              *
*    26 Oct 2005, R. Easter - changes associated with WRF horizontal grid.     *
*                 x and y coordinates are in m, so the clustering              *
*                 calculations are simpler, with no coordinate conversions.    *
c    10 Mar 2006, R. Easter - bug fix at (new) lines 131-2                     *
c                 change "yclust(j)" to "yclust(nclust(i))", same for xclust   *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* fclust          fraction of particles belonging to each cluster              *
* ncluster        number of clusters to be used                                *
* rms             total horizontal rms distance after clustering               *
* rmsclust        horizontal rms distance for each individual cluster          *
* zrms            total vertical rms distance after clustering                 *
* xclust,yclust,  Cluster centroid positions                                   *
* zclust                                                                       *
* xl,yl,zl        particle positions                                           *
*                                                                              *
*                                                                              *
********************************************************************************

      include 'includepar'

      integer n,i,j,l,nclust(maxpart),numb(ncluster),ncl
      real xl(n),yl(n),zl(n),xclust(ncluster),yclust(ncluster),x,y,z
      real zclust(ncluster),distance2,distances,distancemin,rms,rmsold
      real xav(ncluster),yav(ncluster),zav(ncluster),fclust(ncluster)
      real rmsclust(ncluster)
      real zdist,zrms



      if (n.lt.ncluster) return
      rmsold=-5.

C Convert longitude and latitude from degrees to radians
********************************************************

      do 5 i=1,n
        nclust(i)=i

c for FLEXPART_WRF, x & y coords are in meters
c        xl(i)=xl(i)*pi180
c5       yl(i)=yl(i)*pi180
5       continue


C Generate a seed for each cluster
**********************************

      do 10 j=1,ncluster
        zclust(j)=0.
        xclust(j)=xl(j*n/ncluster)
10      yclust(j)=yl(j*n/ncluster)


C Iterative loop to compute the cluster means
*********************************************

      do 100 l=1,100

C Assign each particle to a cluster: criterion minimum distance to the
C cluster mean position
**********************************************************************

      
        do 20 i=1,n
          distancemin=10.**10.
          do 30 j=1,ncluster

c for FLEXPART_WRF, x & y coords are in meters, so calc distance directly
c           distances=distance2(yl(i),xl(i),yclust(j),xclust(j))
            distances=sqrt( (yl(i)-yclust(j))**2 + 
     +                      (xl(i)-xclust(j))**2 )

            if (distances.lt.distancemin) then
              distancemin=distances
              ncl=j
            endif
30          continue
          nclust(i)=ncl
20        continue


C Recalculate the cluster centroid position: convert to 3D Cartesian coordinates,
C calculate mean position, and re-project this point onto the Earth's surface
*********************************************************************************

        do 40 j=1,ncluster
          xav(j)=0.
          yav(j)=0.
          zav(j)=0.
          rmsclust(j)=0.
40        numb(j)=0
        rms=0.

        do 50 i=1,n
          numb(nclust(i))=numb(nclust(i))+1

c for FLEXPART_WRF, x & y coords are in meters, so calc distance directly
c          distances=distance2(yl(i),xl(i),
c     +    yclust(nclust(i)),xclust(nclust(i)))
c 10-mar-2006 rce - bug fix - change "yclust(j)" to 
c    "yclust(nclust(i))", same for xclust
          distances=sqrt( (yl(i)-yclust(nclust(i)))**2 + 
     +                    (xl(i)-xclust(nclust(i)))**2 )

C rms is the total rms of all particles
C rmsclust is the rms for a particular cluster
**********************************************

          rms=rms+distances*distances
          rmsclust(nclust(i))=rmsclust(nclust(i))+distances*distances

C Calculate Cartesian 3D coordinates from longitude and latitude
****************************************************************

c for FLEXPART_WRF, x & y coords are in meters, 
c so no conversion is needed
c          x = cos(yl(i))*sin(xl(i))
c          y = -1.*cos(yl(i))*cos(xl(i))
c          z = sin(yl(i))
c          xav(nclust(i))=xav(nclust(i))+x
c          yav(nclust(i))=yav(nclust(i))+y
c50        zav(nclust(i))=zav(nclust(i))+z
          xav(nclust(i))=xav(nclust(i))+xl(i)
          yav(nclust(i))=yav(nclust(i))+yl(i)
50        zav(nclust(i))=0.0

        rms=sqrt(rms/float(n))


C Find the mean location in Cartesian coordinates
*************************************************

        do 60 j=1,ncluster
          if (numb(j).gt.0) then
            rmsclust(j)=sqrt(rmsclust(j)/float(numb(j)))
            xav(j)=xav(j)/float(numb(j))
            yav(j)=yav(j)/float(numb(j))
            zav(j)=zav(j)/float(numb(j))

C Project the point back onto Earth's surface
*********************************************

c for FLEXPART_WRF, x & y coords are in meters, 
c so no conversion is needed
c            xclust(j)=atan2(xav(j),-1.*yav(j))
c            yclust(j)=atan2(zav(j),sqrt(xav(j)*xav(j)+yav(j)*yav(j)))
            xclust(j)=xav(j)
            yclust(j)=yav(j)
          endif
60        continue


C Leave the loop if the RMS distance decreases only slightly between 2 iterations
*********************************************************************************

        if ((l.gt.1).and.(abs(rms-rmsold)/rmsold.lt.0.005)) goto 99
        rmsold=rms

100     continue

99    continue

C Convert longitude and latitude from radians to degrees
********************************************************

      do 6 i=1,n
c for FLEXPART_WRF, x & y coords are in meters
c        xl(i)=xl(i)/pi180
c        yl(i)=yl(i)/pi180
6       zclust(nclust(i))=zclust(nclust(i))+zl(i)

      do 7 j=1,ncluster
c for FLEXPART_WRF, x & y coords are in meters
c        xclust(j)=xclust(j)/pi180
c        yclust(j)=yclust(j)/pi180
        if (numb(j).gt.0) zclust(j)=zclust(j)/float(numb(j))
        fclust(j)=100.*float(numb(j))/float(n)
7       continue

C Determine total vertical RMS deviation
****************************************

      zrms=0.
      do 8 i=1,n
        zdist=zl(i)-zclust(nclust(i))
8       zrms=zrms+zdist*zdist
      if (zrms.gt.0.) zrms=sqrt(zrms/float(n))

      end
