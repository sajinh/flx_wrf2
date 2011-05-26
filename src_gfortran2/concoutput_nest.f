      subroutine concoutput_nest(itime,outnum)
C                                  i     i
********************************************************************************
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subroutine concoutput_nest.   *
*                                                                              *
*     Output of the concentration nest, adapted from concoutput.f              *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     24 May 1995                                                              *
*                                                                              *
*     13 April 1999, Major update: if output size is smaller, dump output      *
*                    in sparse matrix format; additional output of uncertainty *
*                                                                              *
*     05 April 2000, Major update: output of age classes; output for backward  *
*                    runs is time spent in grid cell times total mass of       *
*                    species.                                                  *
*                                                                              *
*     17 February 2002, Appropriate dimensions for backward and forward runs   *
*                       are now specified in file includepar                   *
*                                                                              *
*      2 September 2004, Adaptation from concoutput.f to output nested data.   *
*                                                                              *
*     Dec 2005, J. Fast - Output files can be either binary or ascii.          *
*                         Sparse output option is turned off.                  *
*     Dec 2005, R. Easter - changed names of "*lon0*" & "*lat0*" variables     *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* outnum          number of samples                                            *
* ncells          number of cells with non-zero concentrations                 *
* sparse          .true. if in sparse matrix format, else .false.              *
* nspeciesdim     either nspec (forward runs), or numpoint (backward runs)     *
* tot_mu          1 for forward, initial mass mixing ration for backw. runs    *
* maxpointspec    maxspec for forward runs, maxpoint for backward runs         *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      double precision jul
      integer itime,ix,jy,kz,k,l,iix,jjy,kzz,nage,jjjjmmdd,ihmmss
      integer ncells(maxpointspec,maxageclass)
      integer ncellsd(maxpointspec,maxageclass)
      integer ncellsw(maxpointspec,maxageclass),nspeciesdim
      real outnum,weightair,xl,yl
      real densityoutgridn(0:maxxgridn-1,0:maxygridn-1,maxzgrid)
      real gridn(0:maxxgridn-1,0:maxygridn-1,maxzgrid,
     +maxpointspec,maxageclass)
      real wetgridn(0:maxxgridn-1,0:maxygridn-1,maxpointspec,
     +maxageclass)
      real drygridn(0:maxxgridn-1,0:maxygridn-1,maxpointspec,
     +maxageclass)
      real gridsigman(0:maxxgridn-1,0:maxygridn-1,maxzgrid,maxpointspec,
     +maxageclass)
      real drygridsigman(0:maxxgridn-1,0:maxygridn-1,maxpointspec,
     +maxageclass)
      real wetgridsigman(0:maxxgridn-1,0:maxygridn-1,maxpointspec,
     +maxageclass)
      real auxgrid(nclassunc)
      real factorn(0:maxxgridn-1,0:maxygridn-1,maxzgrid)
      real halfheight,dz,dz1,dz2,tot_mu(maxpointspec)
      parameter(weightair=28.97)
      logical sparse(maxpointspec,maxageclass)
      logical sparsed(maxpointspec,maxageclass)
      logical sparsew(maxpointspec,maxageclass)
      character adate*8,atime*6


C Determine current calendar date, needed for the file name
***********************************************************

      jul=bdate+dble(float(itime))/86400.
      call caldate(jul,jjjjmmdd,ihmmss)
      write(adate,'(i8.8)') jjjjmmdd
      write(atime,'(i6.6)') ihmmss
c     write(unitdates,'(a)') adate//atime


C For forward simulations, output fields have dimension MAXSPEC,
C for backward simulations, output fields have dimension MAXPOINT.
C Thus, make loops either about nspec, or about numpoint
******************************************************************

      if (ldirect.eq.1) then
        nspeciesdim=nspec
        do 7 k=1,nspeciesdim
7         tot_mu(k)=1.
      else
        nspeciesdim=numpoint
        do 8 k=1,nspeciesdim
8         tot_mu(k)=xmass(k,1)
      endif


**********************************************************************
C Determine the standard deviation of the mean concentration or mixing
C ratio (uncertainty of the output) and the dry and wet deposition
**********************************************************************

      do 120 k=1,nspeciesdim
        do 120 jy=0,numygridn-1
          do 120 ix=0,numxgridn-1
C WET DEPOSITION
            if (WETDEP) then
              do 124 nage=1,nageclass
                do 122 l=1,nclassunc
122               auxgrid(l)=wetgriduncn(ix,jy,k,l,nage)
                call mean(auxgrid,wetgridn(ix,jy,k,nage),
     +          wetgridsigman(ix,jy,k,nage),nclassunc)
C Multiply by number of classes to get total concentration
                wetgridn(ix,jy,k,nage)=wetgridn(ix,jy,k,nage)*nclassunc
C Calculate standard deviation of the mean
124             wetgridsigman(ix,jy,k,nage)=wetgridsigman(ix,jy,k,nage)*
     +          sqrt(float(nclassunc))
            else
              do 144 nage=1,nageclass
144             wetgridn(ix,jy,k,nage)=0.
            endif

C DRY DEPOSITION
            if (DRYDEP) then
              do 125 nage=1,nageclass
                do 123 l=1,nclassunc
123               auxgrid(l)=drygriduncn(ix,jy,k,l,nage)
                call mean(auxgrid,drygridn(ix,jy,k,nage),
     +          drygridsigman(ix,jy,k,nage),nclassunc)
C Multiply by number of classes to get total concentration
                drygridn(ix,jy,k,nage)=drygridn(ix,jy,k,nage)*nclassunc
C Calculate standard deviation of the mean
125             drygridsigman(ix,jy,k,nage)=drygridsigman(ix,jy,k,nage)*
     +          sqrt(float(nclassunc))
            else
              do 145 nage=1,nageclass
145             drygridn(ix,jy,k,nage)=0.
            endif

C CONCENTRATION OR MIXING RATIO
            do 120 kz=1,numzgrid
              do 120 nage=1,nageclass
                do 121 l=1,nclassunc
121               auxgrid(l)=griduncn(ix,jy,kz,k,l,nage)
                call mean(auxgrid,gridn(ix,jy,kz,k,nage),
     +          gridsigman(ix,jy,kz,k,nage),nclassunc)
C Multiply by number of classes to get total concentration
                gridn(ix,jy,kz,k,nage)=gridn(ix,jy,kz,k,nage)*nclassunc
C Calculate standard deviation of the mean
                gridsigman(ix,jy,kz,k,nage)=gridsigman(ix,jy,kz,k,nage)*
     +          sqrt(float(nclassunc))
120          continue


***************************************************************
C Check, whether output of full grid or sparse matrix format is
C more efficient in terms of storage space. This is checked for
C every species and for every age class
***************************************************************

      do 10 k=1,nspeciesdim
        do 10 nage=1,nageclass
          ncellsw(k,nage)=0
          ncellsd(k,nage)=0
10        ncells(k,nage)=0

      do 20 k=1,nspeciesdim
        do 20 nage=1,nageclass
          do 20 jy=0,numygridn-1
            do 20 ix=0,numxgridn-1
              if (wetgridn(ix,jy,k,nage).gt.0) ncellsw(k,nage)=
     +        ncellsw(k,nage)+1
              if (drygridn(ix,jy,k,nage).gt.0) ncellsd(k,nage)=
     +        ncellsd(k,nage)+1
              do 20 kz=1,numzgrid
                if (gridn(ix,jy,kz,k,nage).gt.0) ncells(k,nage)=
     +          ncells(k,nage)+1
20            continue

C Output in sparse matrix format is more efficient, if less than
C 2/5 of all cells contains concentrations>0, because one line
C of sparse matrix output contains (3+2)*4 byte, whereas the equivalent
C in full grid format occupies 2*4 byte. However, note that this is
C compiler-(option)dependent. With some compiler(s)(settings), sparse
C matrix output would require only 3*4 byte per entry (no extra bytes per
C line). Then, 2.5 below should be changed to 1.5 for optimum efficiency.
*************************************************************************

      do 15 k=1,nspeciesdim
        do 15 nage=1,nageclass
          if (4.0*ncellsw(k,nage).lt.numxgridn*numygridn) then
            sparsew(k,nage)=.true.
          else
            sparsew(k,nage)=.false.
          endif
          if (4.0*ncellsd(k,nage).lt.numxgridn*numygridn) then
            sparsed(k,nage)=.true.
          else
            sparsed(k,nage)=.false.
          endif
          if (4.0*ncells(k,nage).lt.numxgridn*numygridn*numzgrid) then
            sparse(k,nage)=.true.
          else
            sparse(k,nage)=.false.
          endif
15        continue


********************************************************************
C Compute air density: sufficiently accurate to take it
C from coarse grid at some time
C Determine center altitude of output layer, and interpolate density
C data to that altitude
********************************************************************

      do 44 kz=1,numzgrid
        if (kz.eq.1) then
          halfheight=outheight(1)/2.
        else
          halfheight=(outheight(kz)+outheight(kz-1))/2.
        endif
        do 45 kzz=2,nz
          if ((height(kzz-1).lt.halfheight).and.
     +    (height(kzz).gt.halfheight)) goto 46
45        continue
46      kzz=max(min(kzz,nz),2)
        dz1=halfheight-height(kzz-1)
        dz2=height(kzz)-halfheight
        dz=dz1+dz2
        do 44 jy=0,numygridn-1
          do 44 ix=0,numxgridn-1
            xl=out_xm0n+float(ix)*dxoutn
            yl=out_ym0n+float(jy)*dyoutn
            xl=(xl-xmet0)/dx
            yl=(yl-ymet0)/dx
            iix=max(min(nint(xl),nxmin1),0)
            jjy=max(min(nint(yl),nymin1),0)
44          densityoutgridn(ix,jy,kz)=(rho(iix,jjy,kzz,2)*dz1+
     +      rho(iix,jjy,kzz-1,2)*dz2)/dz


********************************************************************
C Generate output: may be in concentration (ng/m3) or in mixing
C ratio (ppt) or both
C Output either in full grid dump or sparse matrix format
C For backward simulations, the unit is seconds, stored in grid_conc
********************************************************************

C Concentration output
**********************

      if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then
        if (ldirect.eq.1) then
        if(iouttype.eq.0) then
          open(unitoutgrid,file=path(2)(1:len(2))//'grid_conc_nest_'//
     +    adate//atime,form='unformatted')
        endif
        if(iouttype.eq.1) then
          open(unitoutgrid,file=path(2)(1:len(2))//'grid_conc_nest_'//
     +    adate//atime,form='formatted')
        endif
        else
        if(iouttype.eq.0) then
          open(unitoutgrid,file=path(2)(1:len(2))//'grid_time_nest_'//
     +    adate//atime,form='unformatted')
        endif
        if(iouttype.eq.1) then
          open(unitoutgrid,file=path(2)(1:len(2))//'grid_time_nest_'//
     +    adate//atime,form='formatted')
        endif
        endif


C Output is different for forward and backward simulations
        do 27 kz=1,numzgrid
          do 27 jy=0,numygridn-1
            do 27 ix=0,numxgridn-1
              if (ldirect.eq.1) then
                factorn(ix,jy,kz)=1.e12/volumen(ix,jy,kz)/outnum
              else
                factorn(ix,jy,kz)=float(abs(loutaver))/outnum
              endif
27            continue

        if(iouttype.eq.0) write(unitoutgrid) itime
        if(iouttype.eq.1) write(unitoutgrid,*) itime
        do 30 k=1,nspeciesdim
          do 30 nage=1,nageclass

C Wet deposition
            if (iouttype .eq. 0) then 
            if (sparsew(k,nage)) then
              write(unitoutgrid) 1
              do 31 jy=0,numygridn-1
                do 31 ix=0,numxgridn-1
                  if (wetgridn(ix,jy,k,nage).gt.0.) write(unitoutgrid)
     +            ix+jy*numxgridn,1.e12*wetgridn(ix,jy,k,nage)/
     +            arean(ix,jy)
cc   +            ,1.e12*wetgridsigman(ix,jy,k,nage)/arean(ix,jy)
31                continue
              write(unitoutgrid) -999,999.,999.
            else
              write(unitoutgrid) 2
              do 32 ix=0,numxgridn-1
32              write(unitoutgrid) (1.e12*wetgridn(ix,jy,k,nage)/
     +          arean(ix,jy)
cc   +          ,1.e12*wetgridsigman(ix,jy,k,nage)/arean(ix,jy)
     +          ,jy=0,numygridn-1)
            endif
            endif
            if (iouttype .eq. 1) then 
            if (sparsew(k,nage)) then
              write(unitoutgrid,*) 1
              do 311 jy=0,numygridn-1
                do 311 ix=0,numxgridn-1
                  if (wetgridn(ix,jy,k,nage).gt.0.) write(unitoutgrid,*)
     +            ix+jy*numxgridn,1.e12*wetgridn(ix,jy,k,nage)/
     +            arean(ix,jy)
cc   +            ,1.e12*wetgridsigman(ix,jy,k,nage)/arean(ix,jy)
311               continue
              write(unitoutgrid,*) -999,999.,999.
            else
              write(unitoutgrid,*) 2
              do 321 jy=0,numygridn-1
              do 321 ix=0,numxgridn-1
321             write(unitoutgrid,*) (1.e12*wetgridn(ix,jy,k,nage)/
     +          arean(ix,jy))
cc   +          ,1.e12*wetgridsigman(ix,jy,k,nage)/arean(ix,jy)
            endif
            endif

C Dry deposition
            if (iouttype.eq.0) then
            if (sparsed(k,nage)) then
              write(unitoutgrid) 1
              do 33 jy=0,numygridn-1
                do 33 ix=0,numxgridn-1
                  if (drygridn(ix,jy,k,nage).gt.0.) write(unitoutgrid)
     +            ix+jy*numxgridn,1.e12*drygridn(ix,jy,k,nage)/
     +            arean(ix,jy)
cc   +            ,1.e12*drygridsigman(ix,jy,k,nage)/arean(ix,jy)
33                continue
              write(unitoutgrid) -999,999.,999.
            else
              write(unitoutgrid) 2
              do 34 ix=0,numxgridn-1
34              write(unitoutgrid) (1.e12*drygridn(ix,jy,k,nage)/
     +          arean(ix,jy)
cc   +          ,1.e12*drygridsigman(ix,jy,k,nage)/arean(ix,jy)
     +          ,jy=0,numygridn-1)
            endif
            endif
            if (iouttype.eq.1) then
            if (sparsed(k,nage)) then
              write(unitoutgrid,*) 1
              do 331 jy=0,numygridn-1
                do 331 ix=0,numxgridn-1
                  if (drygridn(ix,jy,k,nage).gt.0.) write(unitoutgrid,*)
     +            ix+jy*numxgridn,1.e12*drygridn(ix,jy,k,nage)/
     +            arean(ix,jy)
cc   +            ,1.e12*drygridsigman(ix,jy,k,nage)/arean(ix,jy)
331               continue
              write(unitoutgrid,*) -999,999.,999.
            else
              write(unitoutgrid,*) 2
              do 341 jy=0,numygridn-1
              do 341 ix=0,numxgridn-1
341             write(unitoutgrid,*) (1.e12*drygridn(ix,jy,k,nage)/
     +          arean(ix,jy))
cc   +          ,1.e12*drygridsigman(ix,jy,k,nage)/arean(ix,jy)
            endif
            endif

C Concentrations
            if (iouttype.eq.0) then 
            if (sparse(k,nage)) then
              write(unitoutgrid) 1
              do 35 kz=1,numzgrid
                do 35 jy=0,numygridn-1
                  do 35 ix=0,numxgridn-1
                    if (gridn(ix,jy,kz,k,nage).gt.0.) write(unitoutgrid)
     +              ix+jy*numxgridn+kz*numxgridn*numygridn,
     +              gridn(ix,jy,kz,k,nage)*factorn(ix,jy,kz)/tot_mu(k) 
cc   +              ,gridsigman(ix,jy,kz,k,nage)*factorn(ix,jy,kz)/
cc   +              tot_mu(k)
35                  continue
              write(unitoutgrid) -999,999.,999.
            else
              write(unitoutgrid) 2
              do 36 kz=1,numzgrid
                do 36 ix=0,numxgridn-1
36                write(unitoutgrid) (gridn(ix,jy,kz,k,nage)*
     +            factorn(ix,jy,kz)/tot_mu(k)
cc   +            ,gridsigman(ix,jy,kz,k,nage)*factorn(ix,jy,kz)/
cc   +            tot_mu(k)
     +            ,jy=0,numygridn-1)
            endif
            endif
            if (iouttype.eq.1) then 
            if (sparse(k,nage)) then
              write(unitoutgrid,*) 1
              do 351 kz=1,numzgrid
                do 351 jy=0,numygridn-1
                  do 351 ix=0,numxgridn-1
                  if (gridn(ix,jy,kz,k,nage).gt.0.) write(unitoutgrid,*)
     +              ix+jy*numxgridn+kz*numxgridn*numygridn,
     +              gridn(ix,jy,kz,k,nage)*factorn(ix,jy,kz)/tot_mu(k) 
cc   +              ,gridsigman(ix,jy,kz,k,nage)*factorn(ix,jy,kz)/
cc   +              tot_mu(k)
351                 continue
              write(unitoutgrid,*) -999,999.,999.
            else
              write(unitoutgrid,*) 2
              do 361 kz=1,numzgrid
                do 361 jy=0,numygridn-1
                do 361 ix=0,numxgridn-1
361               write(unitoutgrid,*) (gridn(ix,jy,kz,k,nage)*
     +            factorn(ix,jy,kz)/tot_mu(k))
cc   +            ,gridsigman(ix,jy,kz,k,nage)*factorn(ix,jy,kz)/
cc   +            tot_mu(k)
            endif
            endif
30          continue

        close(unitoutgrid)

      endif


C Mixing ratio output
*********************

      if ((iout.eq.2).or.(iout.eq.3)) then      ! mixing ratio
        if (iouttype.eq.0) 
     +  open(unitoutgridppt,file=path(2)(1:len(2))//'grid_pptv_nest_'//
     +  adate//atime,form='unformatted')
        if (iouttype.eq.1) 
     +  open(unitoutgridppt,file=path(2)(1:len(2))//'grid_pptv_nest_'//
     +  adate//atime,form='formatted')


        if (iouttype.eq.0) write(unitoutgridppt) itime
        if (iouttype.eq.1) write(unitoutgridppt) itime
        do 130 k=1,nspeciesdim
          do 130 nage=1,nageclass

C Wet deposition
            if (iouttype.eq.0) then
            if (sparsew(k,nage)) then
              write(unitoutgridppt) 1
              do 131 jy=0,numygridn-1
                do 131 ix=0,numxgridn-1
                 if (wetgridn(ix,jy,k,nage).gt.0.) write(unitoutgridppt)
     +            ix+jy*numxgridn,1.e12*wetgridn(ix,jy,k,nage)/
     +            arean(ix,jy)
cc   +            ,1.e12*wetgridsigman(ix,jy,k,nage)/arean(ix,jy)
131               continue
              write(unitoutgridppt) -999,999.,999.
            else
              write(unitoutgridppt) 2
              do 132 ix=0,numxgridn-1
132             write(unitoutgridppt) (1.e12*wetgridn(ix,jy,k,nage)/
     +          arean(ix,jy)
cc   +          ,1.e12*wetgridsigman(ix,jy,k,nage)/arean(ix,jy)
     +          ,jy=0,numygridn-1)
            endif
            endif
            if (iouttype.eq.1) then
            if (sparsew(k,nage)) then
              write(unitoutgridppt,*) 1
              do 1311 jy=0,numygridn-1
                do 1311 ix=0,numxgridn-1
               if (wetgridn(ix,jy,k,nage).gt.0.) write(unitoutgridppt,*)
     +            ix+jy*numxgridn,1.e12*wetgridn(ix,jy,k,nage)/
     +            arean(ix,jy)
cc   +            ,1.e12*wetgridsigman(ix,jy,k,nage)/arean(ix,jy)
1311              continue
              write(unitoutgridppt,*) -999,999.,999.
            else
              write(unitoutgridppt,*) 2
              do 1321 jy=0,numygridn-1
              do 1321 ix=0,numxgridn-1
1321            write(unitoutgridppt,*) (1.e12*wetgridn(ix,jy,k,nage)/
     +          arean(ix,jy))
cc   +          ,1.e12*wetgridsigman(ix,jy,k,nage)/arean(ix,jy)
            endif
            endif

C Dry deposition
            if (iouttype.eq.0) then
            if (sparsed(k,nage)) then
              write(unitoutgridppt) 1
              do 133 jy=0,numygridn-1
                do 133 ix=0,numxgridn-1
                 if (drygridn(ix,jy,k,nage).gt.0.) write(unitoutgridppt)
     +            ix+jy*numxgridn,1.e12*drygridn(ix,jy,k,nage)/
     +            arean(ix,jy)
cc   +            ,1.e12*drygridsigman(ix,jy,k,nage)/arean(ix,jy)
133               continue
              write(unitoutgridppt) -999,999.,999.
            else
              write(unitoutgridppt) 2
              do 134 ix=0,numxgridn-1
134             write(unitoutgridppt) (1.e12*drygridn(ix,jy,k,nage)/
     +          arean(ix,jy)
cc   +          ,1.e12*drygridsigman(ix,jy,k,nage)/arean(ix,jy)
     +          ,jy=0,numygridn-1)
            endif
            endif
            if (iouttype.eq.1) then
            if (sparsed(k,nage)) then
              write(unitoutgridppt,*) 1
              do 1331 jy=0,numygridn-1
                do 1331 ix=0,numxgridn-1
               if (drygridn(ix,jy,k,nage).gt.0.) write(unitoutgridppt,*)
     +            ix+jy*numxgridn,1.e12*drygridn(ix,jy,k,nage)/
     +            arean(ix,jy)
cc   +            ,1.e12*drygridsigman(ix,jy,k,nage)/arean(ix,jy)
1331              continue
              write(unitoutgridppt,*) -999,999.,999.
            else
              write(unitoutgridppt,*) 2
              do 1341 jy=0,numygridn-1
              do 1341 ix=0,numxgridn-1
1341            write(unitoutgridppt,*) (1.e12*drygridn(ix,jy,k,nage)/
     +          arean(ix,jy))
cc   +          ,1.e12*drygridsigman(ix,jy,k,nage)/arean(ix,jy)
            endif
            endif

C Mixing ratios
            if (iouttype.eq.0) then
            if (sparse(k,nage)) then
              write(unitoutgridppt) 1
              do 135 kz=1,numzgrid
                do 135 jy=0,numygridn-1
                  do 135 ix=0,numxgridn-1
                    if (gridn(ix,jy,kz,k,nage).gt.0.)
     +              write(unitoutgridppt)
     +              ix+jy*numxgridn+kz*numxgridn*numygridn,
     +              1.e12*gridn(ix,jy,kz,k,nage)/volumen(ix,jy,kz)/
     +              outnum*
     +              weightair/weightmolar(k)/densityoutgridn(ix,jy,kz)
cc   +             ,1.e12*gridsigman(ix,jy,kz,k,nage)/volumen(ix,jy,kz)/
cc   +              outnum*weightair/weightmolar(k)/
cc   +              densityoutgridn(ix,jy,kz)
135                 continue
              write(unitoutgridppt) -999,999.,999.
            else
              write(unitoutgridppt) 2
              do 136 kz=1,numzgrid
                do 136 ix=0,numxgridn-1
136               write(unitoutgridppt) (1.e12*gridn(ix,jy,kz,k,nage)/
     +            volumen(ix,jy,kz)/outnum*weightair/weightmolar(k)/
     +            densityoutgridn(ix,jy,kz) 
cc   +            ,1.e12*gridsigman(ix,jy,kz,k,nage)/
cc   +            volumen(ix,jy,kz)/outnum*weightair/weightmolar(k)/
cc   +            densityoutgridn(ix,jy,kz)
     +            ,jy=0,numygridn-1)
            endif
            endif
            if (iouttype.eq.1) then
            if (sparse(k,nage)) then
              write(unitoutgridppt,*) 1
              do 1351 kz=1,numzgrid
                do 1351 jy=0,numygridn-1
                  do 1351 ix=0,numxgridn-1
                    if (gridn(ix,jy,kz,k,nage).gt.0.)
     +              write(unitoutgridppt,*)
     +              ix+jy*numxgridn+kz*numxgridn*numygridn,
     +              1.e12*gridn(ix,jy,kz,k,nage)/volumen(ix,jy,kz)/
     +              outnum*
     +              weightair/weightmolar(k)/densityoutgridn(ix,jy,kz)
cc   +             ,1.e12*gridsigman(ix,jy,kz,k,nage)/volumen(ix,jy,kz)/
cc   +              outnum*weightair/weightmolar(k)/
cc   +              densityoutgridn(ix,jy,kz)
1351                continue
              write(unitoutgridppt,*) -999,999.,999.
            else
              write(unitoutgridppt,*) 2
              do 1361 kz=1,numzgrid
                do 1361 jy=0,numygridn-1
                do 1361 ix=0,numxgridn-1
1361              write(unitoutgridppt,*) (1.e12*gridn(ix,jy,kz,k,nage)/
     +            volumen(ix,jy,kz)/outnum*weightair/weightmolar(k)/
     +            densityoutgridn(ix,jy,kz))
cc   +            ,1.e12*gridsigman(ix,jy,kz,k,nage)/
cc   +            volumen(ix,jy,kz)/outnum*weightair/weightmolar(k)/
cc   +            densityoutgridn(ix,jy,kz)
            endif
            endif
130         continue

        close(unitoutgridppt)

      endif


C Reinitialization of grid
**************************

      do 40 k=1,nspeciesdim
        do 40 jy=0,numygridn-1
          do 40 ix=0,numxgridn-1
            do 40 l=1,nclassunc
              do 40 kz=1,numzgrid
                do 40 nage=1,nageclass
40                griduncn(ix,jy,kz,k,l,nage)=0.


      end
