      subroutine fluxoutput(itime)
C                             i
********************************************************************************
*                                                                              *
*     Output of the gridded fluxes.                                            *
*     Eastward, westward, northward, southward, upward and downward gross      *
*     fluxes are written to output file in either sparse matrix or grid dump   *
*     format, whichever is more efficient.                                     *
*                                                                              *
*     Author: A. Stohl                                                         *
*                                                                              *
*     04 April 2000                                                            *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
* ncellse         number of cells with non-zero values for eastward fluxes     *
* sparsee         .true. if in sparse matrix format, else .false.              *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'

      double precision jul
      integer itime,ix,jy,kz,k,nage,jjjjmmdd,ihmmss
      integer ncellse(maxspec,maxageclass),ncellsw(maxspec,maxageclass)
      integer ncellss(maxspec,maxageclass),ncellsn(maxspec,maxageclass)
      integer ncellsu(maxspec,maxageclass),ncellsd(maxspec,maxageclass)
      logical sparsee(maxspec,maxageclass),sparsew(maxspec,maxageclass)
      logical sparses(maxspec,maxageclass),sparsen(maxspec,maxageclass)
      logical sparseu(maxspec,maxageclass),sparsed(maxspec,maxageclass)
      character adate*8,atime*6


C Determine current calendar date, needed for the file name
***********************************************************

      jul=bdate+dble(float(itime))/86400.
      call caldate(jul,jjjjmmdd,ihmmss)
      write(adate,'(i8.8)') jjjjmmdd
      write(atime,'(i6.6)') ihmmss


      open(unitflux,file=path(2)(1:len(2))//'grid_flux_'//adate//
     +atime,form='unformatted')

***************************************************************
C Check, whether output of full grid or sparse matrix format is
C more efficient in terms of storage space. This is checked for
C every species and for every age class
***************************************************************

      do 10 k=1,nspec
        do 10 nage=1,nageclass
          ncellse(k,nage)=0
          ncellsw(k,nage)=0
          ncellsn(k,nage)=0
          ncellss(k,nage)=0
          ncellsu(k,nage)=0
10        ncellsd(k,nage)=0

      do 20 k=1,nspec
        do 20 nage=1,nageclass
          do 20 jy=0,numygrid-1
            do 20 ix=0,numxgrid-1
              do 20 kz=1,numzgrid
                if (fluxe(ix,jy,kz,k,nage).gt.0) ncellse(k,nage)=
     +          ncellse(k,nage)+1
                if (fluxw(ix,jy,kz,k,nage).gt.0) ncellsw(k,nage)=
     +          ncellsw(k,nage)+1
                if (fluxn(ix,jy,kz,k,nage).gt.0) ncellsn(k,nage)=
     +          ncellsn(k,nage)+1
                if (fluxs(ix,jy,kz,k,nage).gt.0) ncellss(k,nage)=
     +          ncellss(k,nage)+1
                if (fluxu(ix,jy,kz,k,nage).gt.0) ncellsu(k,nage)=
     +          ncellsu(k,nage)+1
                if (fluxd(ix,jy,kz,k,nage).gt.0) ncellsd(k,nage)=
     +          ncellsd(k,nage)+1
20            continue

C Output in sparse matrix format more efficient, if less than
C 2/5 of all cells contains concentrations>0
*************************************************************

      do 15 k=1,nspec
        do 15 nage=1,nageclass
          if (4*ncellse(k,nage).lt.numxgrid*numygrid*numzgrid) then
            sparsee(k,nage)=.true.
          else
            sparsee(k,nage)=.false.
          endif
          if (4*ncellsw(k,nage).lt.numxgrid*numygrid*numzgrid) then
            sparsew(k,nage)=.true.
          else
            sparsew(k,nage)=.false.
          endif
          if (4*ncellsn(k,nage).lt.numxgrid*numygrid*numzgrid) then
            sparsen(k,nage)=.true.
          else
            sparsen(k,nage)=.false.
          endif
          if (4*ncellss(k,nage).lt.numxgrid*numygrid*numzgrid) then
            sparses(k,nage)=.true.
          else
            sparses(k,nage)=.false.
          endif
          if (4*ncellsu(k,nage).lt.numxgrid*numygrid*numzgrid) then
            sparseu(k,nage)=.true.
          else
            sparseu(k,nage)=.false.
          endif
          if (4*ncellsd(k,nage).lt.numxgrid*numygrid*numzgrid) then
            sparsed(k,nage)=.true.
          else
            sparsed(k,nage)=.false.
          endif
15        continue



C Flux output: divide by area and time to get flux in ng/m2/s
*************************************************************

      write(unitflux) itime
      do 30 k=1,nspec
        do 30 nage=1,nageclass

          if (sparsee(k,nage)) then
            write(unitflux) 1
            do 35 kz=1,numzgrid
              do 35 jy=0,numygrid-1
                do 35 ix=0,numxgrid-1
                  if (fluxe(ix,jy,kz,k,nage).gt.0.) write(unitflux)
     +              ix+jy*numxgrid+kz*numxgrid*numygrid,1.e12*
     +              fluxe(ix,jy,kz,k,nage)/areaeast(ix,jy,kz)/outstep
35                continue
            write(unitflux) -999,999.
          else
            write(unitflux) 2
            do 36 kz=1,numzgrid
              do 36 ix=0,numxgrid-1
36              write(unitflux) (1.e12*fluxe(ix,jy,kz,k,nage)/
     +          areaeast(ix,jy,kz)/outstep,jy=0,numygrid-1)
          endif

          if (sparsew(k,nage)) then
            write(unitflux) 1
            do 45 kz=1,numzgrid
              do 45 jy=0,numygrid-1
                do 45 ix=0,numxgrid-1
                  if (fluxw(ix,jy,kz,k,nage).gt.0.) write(unitflux)
     +              ix+jy*numxgrid+kz*numxgrid*numygrid,1.e12*
     +              fluxw(ix,jy,kz,k,nage)/areaeast(ix,jy,kz)/outstep
45                continue
            write(unitflux) -999,999.
          else
            write(unitflux) 2
            do 46 kz=1,numzgrid
              do 46 ix=0,numxgrid-1
46              write(unitflux) (1.e12*fluxw(ix,jy,kz,k,nage)/
     +          areaeast(ix,jy,kz)/outstep,jy=0,numygrid-1)
          endif

          if (sparses(k,nage)) then
            write(unitflux) 1
            do 65 kz=1,numzgrid
              do 65 jy=0,numygrid-1
                do 65 ix=0,numxgrid-1
                  if (fluxs(ix,jy,kz,k,nage).gt.0.) write(unitflux)
     +              ix+jy*numxgrid+kz*numxgrid*numygrid,1.e12*
     +              fluxs(ix,jy,kz,k,nage)/areanorth(ix,jy,kz)/outstep
65                continue
            write(unitflux) -999,999.
          else
            write(unitflux) 2
            do 66 kz=1,numzgrid
              do 66 ix=0,numxgrid-1
66              write(unitflux) (1.e12*fluxs(ix,jy,kz,k,nage)/
     +          areanorth(ix,jy,kz)/outstep,jy=0,numygrid-1)
          endif

          if (sparsen(k,nage)) then
            write(unitflux) 1
            do 55 kz=1,numzgrid
              do 55 jy=0,numygrid-1
                do 55 ix=0,numxgrid-1
                  if (fluxn(ix,jy,kz,k,nage).gt.0.) write(unitflux)
     +              ix+jy*numxgrid+kz*numxgrid*numygrid,1.e12*
     +              fluxn(ix,jy,kz,k,nage)/areanorth(ix,jy,kz)/outstep
55                continue
            write(unitflux) -999,999.
          else
            write(unitflux) 2
            do 56 kz=1,numzgrid
              do 56 ix=0,numxgrid-1
56              write(unitflux) (1.e12*fluxn(ix,jy,kz,k,nage)/
     +          areanorth(ix,jy,kz)/outstep,jy=0,numygrid-1)
          endif

          if (sparseu(k,nage)) then
            write(unitflux) 1
            do 75 kz=1,numzgrid
              do 75 jy=0,numygrid-1
                do 75 ix=0,numxgrid-1
                  if (fluxu(ix,jy,kz,k,nage).gt.0.) write(unitflux)
     +              ix+jy*numxgrid+kz*numxgrid*numygrid,1.e12*
     +              fluxu(ix,jy,kz,k,nage)/area(ix,jy)/outstep
75                continue
            write(unitflux) -999,999.
          else
            write(unitflux) 2
            do 76 kz=1,numzgrid
              do 76 ix=0,numxgrid-1
76              write(unitflux) (1.e12*fluxu(ix,jy,kz,k,nage)/
     +          area(ix,jy)/outstep,jy=0,numygrid-1)
          endif

          if (sparsed(k,nage)) then
            write(unitflux) 1
            do 85 kz=1,numzgrid
              do 85 jy=0,numygrid-1
                do 85 ix=0,numxgrid-1
                  if (fluxd(ix,jy,kz,k,nage).gt.0.) write(unitflux)
     +              ix+jy*numxgrid+kz*numxgrid*numygrid,1.e12*
     +              fluxd(ix,jy,kz,k,nage)/area(ix,jy)/outstep
85                continue
            write(unitflux) -999,999.
          else
            write(unitflux) 2
            do 86 kz=1,numzgrid
              do 86 ix=0,numxgrid-1
86              write(unitflux) (1.e12*fluxd(ix,jy,kz,k,nage)/
     +          area(ix,jy)/outstep,jy=0,numygrid-1)
          endif

30        continue


      close(unitflux)


C Reinitialization of grid
**************************

      do 40 k=1,nspec
        do 40 jy=0,numygrid-1
          do 40 ix=0,numxgrid-1
              do 40 kz=1,numzgrid
                do 40 nage=1,nageclass
                  fluxe(ix,jy,kz,k,nage)=0.
                  fluxw(ix,jy,kz,k,nage)=0.
                  fluxn(ix,jy,kz,k,nage)=0.
                  fluxs(ix,jy,kz,k,nage)=0.
                  fluxu(ix,jy,kz,k,nage)=0.
40                fluxd(ix,jy,kz,k,nage)=0.


      end
