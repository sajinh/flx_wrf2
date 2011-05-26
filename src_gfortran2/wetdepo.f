      subroutine wetdepo(itime,ltsample,loutnext)
C                          i      i        i
********************************************************************************
*                                                                              *
* Calculation of wet deposition using the concept of scavenging coefficients.  *
* For lack of detailed information, washout and rainout are jointly treated.   *
* It is assumed that precipitation does not occur uniformly within the whole   *
* grid cell, but that only a fraction of the grid cell experiences rainfall.   *
* This fraction is parameterized from total cloud cover and rates of large     *
* scale and convective precipitation.                                          *
*                                                                              *
*    Author: A. Stohl                                                          *
*                                                                              *
*    1 December 1996                                                           *
*                                                                              *
* Correction by Petra Seibert, Sept 2002:                                      *
* use centred precipitation data for integration                               *
* Code may not be correct for decay of deposition!                             *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* cc [0-1]           total cloud cover                                         *
* convp [mm/h]       convective precipitation rate                             *
* fraction [0-1]     fraction of grid, for which precipitation occurs          *
* ix,jy              indices of output grid cell for each particle             *
* itime [s]          actual simulation time [s]                                *
* jpart              particle index                                            *
* ldeltat [s]        interval since radioactive decay was computed             *
* lfr, cfr           area fraction covered by precipitation for large scale    *
*                    and convective precipitation (dependent on prec. rate)    *
* loutnext [s]       time for which gridded deposition is next output          *
* loutstep [s]       interval at which gridded deposition is output            *
* lsp [mm/h]         large scale precipitation rate                            *
* ltsample [s]       interval over which mass is deposited                     *
* prec [mm/h]        precipitation rate in subgrid, where precipitation occurs *
* wetdeposit         mass that is wet deposited                                *
* wetgrid            accumulated deposited mass on output grid                 *
* wetscav            scavenging coefficient                                    *
*                                                                              *
* Constants:                                                                   *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      integer jpart,itime,ltsample,loutnext,ldeltat,i,j,k,ix,jy
      integer ngrid,itage,nage
      real xtn,ytn,lsp,convp,cc,fraction,prec,wetscav
      real wetdeposit(maxspec),lfr(5),cfr(5),restmass,smallnum
      save lfr,cfr,smallnum

      data lfr/0.5,0.65,0.8,0.9,0.95/
      data cfr/0.4,0.55,0.7,0.8,0.9/
      data smallnum /1.e-38/ ! should be smallest number that can be handled

C Compute interval since radioactive decay of deposited mass was computed
*************************************************************************

      if (itime.le.loutnext) then
        ldeltat=itime-(loutnext-loutstep)
      else                                  ! first half of next interval
        ldeltat=itime-loutnext
      endif


C Loop over all particles
*************************

      do 20 jpart=1,numpart
cAF        if (itra1(jpart).eq.-999999999) goto 20
         if(ldirect.eq.1)then
            if (itra1(jpart).gt.itime) goto 20
         else
            if (itra1(jpart).lt.itime) goto 20
         endif   
C Determine age class of the particle
        itage=abs(itra1(jpart)-itramem(jpart))
        do 32 nage=1,nageclass
          if (itage.lt.lage(nage)) goto 33
32        continue
33      continue


C Determine which nesting level to be used
******************************************

        ngrid=0
        do 22 j=numbnests,1,-1
          if ((xtra1(jpart).gt.xln(j)).and.(xtra1(jpart).lt.xrn(j)).and.
     +    (ytra1(jpart).gt.yln(j)).and.(ytra1(jpart).lt.yrn(j))) then
            ngrid=j
            goto 23
          endif
22        continue
23      continue


C Determine nested grid coordinates
***********************************

        if (ngrid.gt.0) then
          xtn=(xtra1(jpart)-xln(ngrid))*xresoln(ngrid)
          ytn=(ytra1(jpart)-yln(ngrid))*yresoln(ngrid)
          ix=int(xtn)
          jy=int(ytn)
        else
          ix=int(xtra1(jpart))
          jy=int(ytra1(jpart))
        endif


C Interpolate large scale precipitation, convective precipitation and
C total cloud cover
C Note that interpolated time refers to itime-0.5*ltsample [PS]
*********************************************************************

        if (ngrid.eq.0) then
          call interpol_rain(lsprec,convprec,tcc,nxmax,nymax,
     +    1,nx,ny,memind,sngl(xtra1(jpart)),sngl(ytra1(jpart)),1,
     +    memtime(1),memtime(2),nint(itime-0.5*ltsample),lsp,convp,cc)
        else
          call interpol_rain_nests(lsprecn,convprecn,tccn,
     +    nxmaxn,nymaxn,1,maxnests,ngrid,nxn,nyn,memind,xtn,ytn,1,
     +    memtime(1),memtime(2),nint(itime-0.5*ltsample),lsp,convp,cc)
        endif

        if ((lsp.lt.0.01).and.(convp.lt.0.01)) goto 20

C 1) Parameterization of the the area fraction of the grid cell where the
C    precipitation occurs: the absolute limit is the total cloud cover, but
C    for low precipitation rates, an even smaller fraction of the grid cell
C    is used. Large scale precipitation occurs over larger areas than
C    convective precipitation.
***************************************************************************

        if (lsp.gt.20.) then
          i=5
        else if (lsp.gt.8.) then
          i=4
        else if (lsp.gt.3.) then
          i=3
        else if (lsp.gt.1.) then
          i=2
        else
          i=1
        endif

        if (convp.gt.20.) then
          j=5
        else if (convp.gt.8.) then
          j=4
        else if (convp.gt.3.) then
          j=3
        else if (convp.gt.1.) then
          j=2
        else
          j=1
        endif

        fraction=max(0.05,cc*(lsp*lfr(i)+convp*cfr(j))/(lsp+convp))

C 2) Computation of precipitation rate in sub-grid cell
*******************************************************

        prec=(lsp+convp)/fraction


C 3) Computation of scavenging coefficients for all species
C    Computation of wet deposition
***********************************************************

        do 10 k=1,nspec                                  ! loop over species
          if (weta(k).gt.0.) then
            wetscav=weta(k)*prec**wetb(k)                ! scavenging coeff.
            wetdeposit(k)=xmass1(jpart,k)*
     +        (1.-exp(-wetscav*abs(ltsample)))*fraction  ! wet deposition
C           new particle mass:
            restmass = xmass1(jpart,k)-wetdeposit(k)
            if (restmass .gt. smallnum) then
              xmass1(jpart,k)=restmass
            else
              xmass1(jpart,k)=0.
            endif


C Correct deposited mass to the last time step when radioactive decay of 
C gridded deposited mass was calculated
************************************************************************
            if (decay(k).gt.0.) then
              wetdeposit(k)=wetdeposit(k)*exp(abs(ldeltat)*decay(k))
            endif
          else
            wetdeposit(k)=0.
          endif
10        continue

C Add the wet deposition to accumulated amount on output grid and nested output grid
************************************************************************************

        call wetdepokernel(nclass(jpart),wetdeposit,sngl(xtra1(jpart)),
     +  sngl(ytra1(jpart)),nage)
        if (nested_output.eq.1) call wetdepokernel_nest(nclass(jpart),
     +  wetdeposit,sngl(xtra1(jpart)),sngl(ytra1(jpart)),nage)

20      continue

      end
