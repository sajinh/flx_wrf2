      subroutine wetdepokernel(nunc,deposit,x,y,nage)
C                               i      i    i i  i
********************************************************************************
*                                                                              *
*     Attribution of the deposition from an individual particle to the         *
*     deposition fields using a uniform kernel with bandwidths dxout and dyout.*
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     26 December 1996                                                         *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
*                                                                              *
* nunc             uncertainty class of the respective particle                *
* nage             age class of the respective particle                        *
* deposit          amount (kg) to be deposited                                 *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'


      real x,y,deposit(maxspec),ddx,ddy,xl,yl,wx,wy,w
      integer ix,jy,ixp,jyp,k,nunc,nage
   

      xl=(x*dx+xoutshift)/dxout
      yl=(y*dy+youtshift)/dyout
      ix=int(xl)
      jy=int(yl)
      ddx=xl-float(ix)                   ! distance to left cell border
      ddy=yl-float(jy)                   ! distance to lower cell border

      if (ddx.gt.0.5) then
        ixp=ix+1
        wx=1.5-ddx
      else
        ixp=ix-1
        wx=0.5+ddx
      endif

      if (ddy.gt.0.5) then
        jyp=jy+1
        wy=1.5-ddy
      else
        jyp=jy-1
        wy=0.5+ddy
      endif


C Determine mass fractions for four grid points
***********************************************

      if ((ix.ge.0).and.(jy.ge.0).and.(ix.le.numxgrid-1).and.
     +(jy.le.numygrid-1)) then
        w=wx*wy
        do 22 k=1,nspec
22        wetgridunc(ix,jy,k,nunc,nage)=
     +    wetgridunc(ix,jy,k,nunc,nage)+deposit(k)*w
      endif

      if ((ixp.ge.0).and.(jyp.ge.0).and.(ixp.le.numxgrid-1).and.
     +(jyp.le.numygrid-1)) then
        w=(1.-wx)*(1.-wy)
        do 23 k=1,nspec
23        wetgridunc(ixp,jyp,k,nunc,nage)=
     +    wetgridunc(ixp,jyp,k,nunc,nage)+deposit(k)*w
      endif

      if ((ixp.ge.0).and.(jy.ge.0).and.(ixp.le.numxgrid-1).and.
     +(jy.le.numygrid-1)) then
        w=(1.-wx)*wy
        do 24 k=1,nspec
24        wetgridunc(ixp,jy,k,nunc,nage)=
     +    wetgridunc(ixp,jy,k,nunc,nage)+deposit(k)*w
      endif

      if ((ix.ge.0).and.(jyp.ge.0).and.(ix.le.numxgrid-1).and.
     +(jyp.le.numygrid-1)) then
        w=wx*(1.-wy)
        do 25 k=1,nspec
25        wetgridunc(ix,jyp,k,nunc,nage)=
     +    wetgridunc(ix,jyp,k,nunc,nage)+deposit(k)*w
      endif

      end
