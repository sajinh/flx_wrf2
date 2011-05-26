      subroutine plumetraj(itime)
C                            i
********************************************************************************
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subroutine plumetraj.         *
*            The computational grid is the WRF x-y grid rather than lat-lon.   *
*                                                                              *
* Determines a plume centroid trajectory for each release site, and manages    *
* clustering of particle locations. Certain parameters (average PV,            *
* tropopause height, etc., are provided along the plume trajectories.          *
* At the end, output is written to file 'trajectories.txt'.                    *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     24 January 2002                                                          *
*                                                                              *
*    26 Oct 2005, R. Easter - changes associated with WRF horizontal grid.     *
*                 Calculate the distance between 2 points directly             *
*                 instead of using the distance function.                      *
*     Dec 2005, R. Easter - changed names of "*lon0*" & "*lat0*" variables     *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* fclust          fraction of particles belonging to each cluster              *
* hmixcenter      mean mixing height for all particles                         *
* ncluster        number of clusters to be used                                *
* pvcenter        mean PV for all particles                                    *
* pvfract         fraction of particles with PV<2pvu                           *
* rms             total horizontal rms distance after clustering               *
* rmsdist         total horizontal rms distance before clustering              *
* rmsclust        horizontal rms distance for each individual cluster          *
* topocenter      mean topography underlying all particles                     *
* tropocenter     mean tropopause height at the positions of particles         *
* tropofract      fraction of particles within the troposphere                 *
* zrms            total vertical rms distance after clustering                 *
* zrmsdist        total vertical rms distance before clustering                *
* xclust,yclust,  Cluster centroid positions                                   *
* zclust                                                                       *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      integer itime,ix,jy,ixp,jyp,indexh,i,j,k,m,n,il,ind,indz,indzp
      real xl(maxpart),yl(maxpart),zl(maxpart)
      real xcenter,ycenter,zcenter,dist,distance,rmsdist,zrmsdist

      real xclust(ncluster),yclust(ncluster),zclust(ncluster)
      real fclust(ncluster),rms,rmsclust(ncluster),zrms

      real dt1,dt2,dtt,ddx,ddy,rddx,rddy,p1,p2,p3,p4,dz1,dz2,dz
      real topo,topocenter,hm(2),hmixi,hmixfract,hmixcenter
      real pv1(2),pvprof(2),pvi,pvcenter,pvfract,tr(2),tri,tropofract
      real tropocenter


      dt1=float(itime-memtime(1))
      dt2=float(memtime(2)-itime)
      dtt=1./(dt1+dt2)


C Loop about all release points
*******************************

      do 10 j=1,numpoint
        if (abs(ireleasestart(j)-itime).gt.lage(nageclass)) goto 10
        topocenter=0.
        hmixcenter=0.
        hmixfract=0.
        tropocenter=0.
        tropofract=0.
        pvfract=0.
        pvcenter=0.
        rmsdist=0.
        zrmsdist=0.

        n=0
        do 20 i=1,numpart
          if (itra1(i).ne.itime) goto 20
          if (npoint(i).ne.j) goto 20
          n=n+1
          xl(n)=xmet0+xtra1(i)*dx
          yl(n)=ymet0+ytra1(i)*dy
          zl(n)=ztra1(i)


C Interpolate PBL height, PV, and tropopause height to each
C particle position in order to determine fraction of particles
C within the PBL, above tropopause height, and average PV.     
C Interpolate topography, too, and convert to altitude asl
***************************************************************

          ix=int(xtra1(i))
          jy=int(ytra1(i))
          ixp=ix+1
          jyp=jy+1
          ddx=xtra1(i)-float(ix)
          ddy=ytra1(i)-float(jy)
          rddx=1.-ddx
          rddy=1.-ddy
          p1=rddx*rddy
          p2=ddx*rddy
          p3=rddx*ddy
          p4=ddx*ddy

C Topography
************

          topo=p1*oro(ix ,jy)
     +       + p2*oro(ixp,jy)
     +       + p3*oro(ix ,jyp)
     +       + p4*oro(ixp,jyp)
          topocenter=topocenter+topo

C Potential vorticity
*********************

          do 5 il=2,nz
            if (height(il).gt.zl(n)) then
              indz=il-1
              indzp=il
              goto 6
            endif
5           continue
6         continue

          dz1=zl(n)-height(indz)
          dz2=height(indzp)-zl(n)
          dz=1./(dz1+dz2)


          do 70 ind=indz,indzp
            do 80 m=1,2
              indexh=memind(m)
              pv1(m)=p1*pv(ix ,jy ,ind,indexh)
     +              +p2*pv(ixp,jy ,ind,indexh)
     +              +p3*pv(ix ,jyp,ind,indexh)
     +              +p4*pv(ixp,jyp,ind,indexh)
80            continue
            pvprof(ind-indz+1)=(pv1(1)*dt2+pv1(2)*dt1)*dtt
70          continue
          pvi=(dz1*pvprof(2)+dz2*pvprof(1))*dz
          pvcenter=pvcenter+pvi
          if (yl(n).gt.0.) then
            if (pvi.lt.2.) pvfract=pvfract+1.
          else
            if (pvi.gt.-2.) pvfract=pvfract+1.
          endif


C Tropopause and PBL height
***************************

          do 40 m=1,2
            indexh=memind(m)

            tr(m)=p1*tropopause(ix ,jy ,1,indexh)
     +          + p2*tropopause(ixp,jy ,1,indexh)
     +          + p3*tropopause(ix ,jyp,1,indexh)
     +          + p4*tropopause(ixp,jyp,1,indexh)

40          hm(m)=p1*hmix(ix ,jy ,1,indexh)
     +          + p2*hmix(ixp,jy ,1,indexh)
     +          + p3*hmix(ix ,jyp,1,indexh)
     +          + p4*hmix(ixp,jyp,1,indexh)

          hmixi=(hm(1)*dt2+hm(2)*dt1)*dtt
          tri=(tr(1)*dt2+tr(2)*dt1)*dtt
          if (zl(n).lt.tri) tropofract=tropofract+1.
          tropocenter=tropocenter+tri+topo
          if (zl(n).lt.hmixi) hmixfract=hmixfract+1.
          zl(n)=zl(n)+topo        ! convert to height asl
          hmixcenter=hmixcenter+hmixi


20        continue


C Make statistics for all plumes with n>0 particles
***************************************************

        if (n.gt.0) then
          topocenter=topocenter/float(n)
          hmixcenter=hmixcenter/float(n)
          pvcenter=pvcenter/float(n)
          tropocenter=tropocenter/float(n)
          hmixfract=100.*hmixfract/float(n)
          pvfract=100.*pvfract/float(n)
          tropofract=100.*tropofract/float(n)

C Cluster the particle positions
********************************

          call clustering(xl,yl,zl,n,xclust,yclust,zclust,fclust,rms,
     +    rmsclust,zrms)


C Determine center of mass position on earth and average height
***************************************************************

          call centerofmass(xl,yl,n,xcenter,ycenter)
          call mean(zl,zcenter,zrmsdist,n)

C Root mean square distance from center of mass
***********************************************

          do 30 k=1,n
c for FLEXPART_WRF, x,y coords are in meters, so xl,yl are in meters
c           dist=distance(yl(k),xl(k),ycenter,xcenter)
            dist=sqrt( (yl(k)-ycenter)**2 + (xl(k)-xcenter)**2 )

            rmsdist=rmsdist+dist*dist
30          continue
          if (rmsdist.gt.0.) rmsdist=sqrt(rmsdist/float(n))

C Write out results in trajectory data file
*******************************************

          write(unitouttraj,'(i5,i8,2f9.4,4f8.1,f8.2,4f8.1,3f6.1,
     +    5(2f8.3,f7.0,f6.1,f8.1))')
     +    j,itime-(ireleasestart(j)+ireleaseend(j))/2,
     +    xcenter,ycenter,zcenter,topocenter,hmixcenter,tropocenter,
     +    pvcenter,rmsdist,rms,zrmsdist,zrms,hmixfract,pvfract,
     +    tropofract,
     +    (xclust(k),yclust(k),zclust(k),fclust(k),rmsclust(k),
     +    k=1,ncluster)
        endif
          

10      continue


      end
