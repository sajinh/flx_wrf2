      subroutine calcfluxes(nage,jpart,xold,yold,zold)
C                            i     i    i    i    i
********************************************************************************
*                                                                              *
*     Calculation of the gross fluxes across horizontal, eastward and          *
*     northward facing surfaces. The routine calculates the mass flux          *
*     due to the motion of only one particle. The fluxes of subsequent calls   *
*     to this subroutine are accumulated until the next output is due.         *
*     Upon output, flux fields are re-set to zero in subroutine fluxoutput.f.  *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     04 April 2000                                                            *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
*                                                                              *
* nage                  Age class of the particle considered                   *
* jpart                 Index of the particle considered                       *
* xold,yold,zold        "Memorized" old positions of the particle              *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      integer jpart,nage,ixave,jyave,kz,kzave
      integer k,k1,k2,ix,ix1,ix2,ixs,jy,jy1,jy2
      real xold,yold,zold,xmean,ymean


C Determine average positions
*****************************

      xmean=(xold+xtra1(jpart))/2.
      ymean=(yold+ytra1(jpart))/2.

      ixave=int((xmean*dx+xoutshift)/dxout)
      jyave=int((ymean*dy+youtshift)/dyout)
      do 15 kz=1,numzgrid                ! determine height of cell
        if (outheight(kz).gt.ztra1(jpart)) goto 16
15      continue
16      kzave=kz


C Determine vertical fluxes
***************************

      if ((ixave.ge.0).and.(jyave.ge.0).and.(ixave.le.numxgrid-1).and.
     +(jyave.le.numygrid-1)) then
        do 10 kz=1,numzgrid                ! determine height of cell
          if (outheighthalf(kz).gt.zold) goto 11
10        continue
11      k1=min(numzgrid,kz)
        do 20 kz=1,numzgrid                ! determine height of cell
          if (outheighthalf(kz).gt.ztra1(jpart)) goto 21
20        continue
21      k2=min(numzgrid,kz)

        do 40 k=1,nspec
          do 41 kz=k1,k2-1
41          fluxu(ixave,jyave,kz,k,nage)=fluxu(ixave,jyave,kz,k,nage)+
     +      xmass1(jpart,k)
          do 40 kz=k2,k1-1
40          fluxd(ixave,jyave,kz,k,nage)=fluxd(ixave,jyave,kz,k,nage)+
     +      xmass1(jpart,k)
      endif


C Determine west-east fluxes (fluxw) and east-west fluxes (fluxe)
*****************************************************************

      if ((kzave.le.numzgrid).and.(jyave.ge.0).and.
     +(jyave.le.numygrid-1)) then

C 1) Particle does not cross domain boundary

        if (abs(xold-xtra1(jpart)).lt.float(nx)/2.) then
          ix1=int((xold*dx+xoutshift)/dxout+0.5)
          ix2=int((xtra1(jpart)*dx+xoutshift)/dxout+0.5)
          do 50 k=1,nspec
            do 51 ix=ix1,ix2-1
              if ((ix.ge.0).and.(ix.le.numxgrid-1)) then
                fluxw(ix,jyave,kzave,k,nage)=
     +          fluxw(ix,jyave,kzave,k,nage)+xmass1(jpart,k)
              endif
51            continue
            do 50 ix=ix2,ix1-1
              if ((ix.ge.0).and.(ix.le.numxgrid-1)) then
                fluxe(ix,jyave,kzave,k,nage)=
     +          fluxe(ix,jyave,kzave,k,nage)+xmass1(jpart,k)
              endif
50            continue

C 2) Particle crosses domain boundary: use cyclic boundary condition
C    and attribute flux to easternmost grid row only (approximation valid
C    for relatively slow motions compared to output grid cell size)

        else
          ixs=int(((float(nxmin1)-1.e5)*dx+xoutshift)/dxout)
          if ((ixs.ge.0).and.(ixs.le.numxgrid-1)) then
            if (xold.gt.xtra1(jpart)) then       ! west-east flux
              do 52 k=1,nspec
52              fluxw(ixs,jyave,kzave,k,nage)=
     +          fluxw(ixs,jyave,kzave,k,nage)+xmass1(jpart,k)
            else                                 ! east-west flux
              do 53 k=1,nspec
53              fluxe(ixs,jyave,kzave,k,nage)=
     +          fluxe(ixs,jyave,kzave,k,nage)+xmass1(jpart,k)
            endif
          endif
        endif
      endif


C Determine south-north fluxes (fluxs) and north-south fluxes (fluxn)
*********************************************************************

      if ((kzave.le.numzgrid).and.(ixave.ge.0).and.
     +(ixave.le.numxgrid-1)) then
        jy1=int((yold*dy+youtshift)/dyout+0.5)
        jy2=int((ytra1(jpart)*dy+youtshift)/dyout+0.5)

        do 60 k=1,nspec
          do 61 jy=jy1,jy2-1
            if ((jy.ge.0).and.(jy.le.numygrid-1)) then
              fluxs(ixave,jy,kzave,k,nage)=
     +        fluxs(ixave,jy,kzave,k,nage)+xmass1(jpart,k)
            endif
61          continue
          do 60 jy=jy2,jy1-1
            if ((jy.ge.0).and.(jy.le.numygrid-1)) then
              fluxn(ixave,jy,kzave,k,nage)=
     +        fluxn(ixave,jy,kzave,k,nage)+xmass1(jpart,k)
            endif
60          continue
      endif

      end
