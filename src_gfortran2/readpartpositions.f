      subroutine readpartpositions()
********************************************************************************
*                                                                              *
*   Note:  This is the FLEXPART_WRF version of subroutine readpartpositions.   *
*                                                                              *
*   This routine opens the particle dump file and reads all the particle       *
*   positions from a previous run to initialize the current run.               *
*                                                                              *
*                                                                              *
*     Author: A. Stohl                                                         *
*     24 March 2000                                                            *
*                                                                              *
*     Dec 2005, R. Easter                                                      *
*             Changed names of "*lon0*" & "*lat0*" variables                   *
*             Reads either binary or ascii output files from previous run.     *
*             Particle positions may be in lat-lon or grid-meter units.        *
*                                                                              *
********************************************************************************
*                                                                              *
* Variables:                                                                   *
*                                                                              *
********************************************************************************

      include 'includepar'
      include 'includecom'
      
      integer ibdatein,ibtimein,nspecin,itimein,numpointin,i,j,idummy,ix
      integer iomode_xycoord_in, numpart_in
      integer itmp,ntmp, numxgridin,numygridin
      real xlonin,ylatin,ran1,topo,hmixi,pvi,qvi,rhoi,tri,tti
      real xtmp
      character specin*7
      character ctmp*1
      double precision julin,julpartin,juldate

      data idummy/-8/


C Open and read header file of dumped particle data
******************************************

      if (iouttype .eq. 0) then

      open(unitpartin,file=path(2)(1:len(2))//'header',
     +form='unformatted',err=998)

      read(unitpartin) ibdatein,ibtimein
      read(unitpartin)
      read(unitpartin)

      read(unitpartin)
      read(unitpartin)
      read(unitpartin) nspecin
      nspecin=nspecin/3
      if ((ldirect.eq.1).and.(nspec.ne.nspecin)) goto 997

      do i=1,nspecin
        read(unitpartin)
        read(unitpartin)
        read(unitpartin) j,specin
        if ((ldirect.eq.1).and.(species(i)(1:7).ne.specin)) goto 996
      enddo
 
      read(unitpartin) numpointin
      if (numpointin.ne.numpoint) goto 995
      do i=1,numpointin
        read(unitpartin)
        read(unitpartin)
        read(unitpartin)
        read(unitpartin)
        do j=1,nspec
          read(unitpartin)
          read(unitpartin)
          read(unitpartin)
        enddo
      enddo
      read(unitpartin)
      read(unitpartin)

      do ix=0,numxgrid-1
        read(unitpartin)
      enddo

      close(unitpartin)

      else    ! (iouttype .eq. 1)

      open(unitpartin,file=path(2)(1:len(2))//'header',
     +form='formatted',err=998)

c with formatted header file, need to read every variable
c because some of the "writes" cover multiple lines
c (and the number of lines is compiler dependent)
      read(unitpartin,*) ibdatein,ibtimein
      read(unitpartin,*) ctmp   ! version id
      read(unitpartin,*) itmp,itmp,itmp,itmp   ! loutstep,loutaver,...
      read(unitpartin,*) xtmp,xtmp,numxgridin,numygridin,xtmp,xtmp   ! outgrid x,y info

      read(unitpartin,*) ntmp,(xtmp, i=1,ntmp)   ! numzgrid,(outheight(i),i=1,numzgrid)
      read(unitpartin,*) itmp,itmp   ! jjjjmmdd,ihmmss
      read(unitpartin,*) nspecin,itmp,itmp
      nspecin=nspecin/3
      if ((ldirect.eq.1).and.(nspec.ne.nspecin)) goto 997

      do i=1,nspecin
        read(unitpartin,*) itmp   ! 1
        read(unitpartin,*) ctmp   ! "WD" name
        read(unitpartin,*) itmp   ! 1
        read(unitpartin,*) ctmp   ! "DD" name
        read(unitpartin,*) j
        read(unitpartin,*) specin
        if ((ldirect.eq.1).and.(species(i)(1:7).ne.specin)) goto 996
      enddo
 
      read(unitpartin,*) numpointin
      if (numpointin.ne.numpoint) goto 995
      do i=1,numpointin
        read(unitpartin,*) itmp,itmp,itmp   ! release start,end,kindz
        read(unitpartin,*) xtmp,xtmp,xtmp,xtmp,xtmp,xtmp   ! release x,y,z info
        read(unitpartin,*) itmp,itmp   ! npart(i), 1
        read(unitpartin,*) ctmp   ! compoint(i)
        do j=1,nspec
          read(unitpartin,*) xtmp   ! xmass(i,j)
          read(unitpartin,*) xtmp   ! xmass(i,j)
          read(unitpartin,*) xtmp   ! xmass(i,j)
        enddo
      enddo
      read(unitpartin,*) itmp,itmp,itmp   ! method,lsubgrid,lconvection
      read(unitpartin,*) ntmp,(itmp, i=1,ntmp)   ! nageclass,(lage(i),i=1,nageclass)

      do ix=0,numxgridin-1
        read(unitpartin,*) (xtmp, j=0,numygridin-1)    ! oroout
      enddo

      close(unitpartin)

      endif   ! (iouttype .eq. 0/1)


C Open and read data file of dumped particle data
****************************************

      if (iouttype .eq. 0) then
        open(unitpartin,file=path(2)(1:len(2))//'partposit_end',
     +    form='unformatted',err=998)
      else
        open(unitpartin,file=path(2)(1:len(2))//'partposit_end',
     +    form='formatted',err=998)
      endif

100   continue
      if (iouttype .eq. 0) then
        read(unitpartin,end=99) itimein,numpart_in,
     +    iomode_xycoord_in
      else
        read(unitpartin,*,end=99) itimein,numpart_in,
     +    iomode_xycoord_in
      endif

c iomode_xycoord of previous & current runs must match
      if (iomode_xycoord_in .ne. iomode_xycoord) then
         write(*,'(/a/a/)') '*** readpartpositions fatal error',
     +     'iomode_xycoord from previous & current runs differ'
         stop
      end if

      i=0
200   i=i+1
      if (iouttype .eq. 0) then
        read(unitpartin) npoint(i),xlonin,ylatin,ztra1(i),itramem(i),
     +    topo,pvi,qvi,rhoi,hmixi,tri,tti,(xmass1(i,j),j=1,nspec)
      else
        read(unitpartin,*) npoint(i),itramem(i),xlonin,ylatin,ztra1(i),
     +    topo,pvi,qvi,rhoi,hmixi,tri,tti,(xmass1(i,j),j=1,nspec)
      endif
         
      if (xlonin.eq.-9999.9) goto 100

      if (iomode_xycoord_in .eq. iomode_xycoord_latlon) then
c convert from lat-lon to grid-index coordinates
        call ll_to_xyindex_wrf( xlonin, ylatin, xtra1(i), ytra1(i) )
      else
c convert from grid-meter to grid-index coordinates
        xtra1(i)=(xlonin-xmet0)/dx
        ytra1(i)=(ylatin-ymet0)/dy
      endif
      goto 200

99    numpart=i-1

      close(unitpartin)


C Set nclass, idt, itra1, itramem, itrasplit to be consistent
C with current run
****************************************

      julin=juldate(ibdatein,ibtimein)+dble(float(itimein)/86400.)
      if (abs(julin-bdate).gt.1.e-5) goto 994
      do 50 i=1,numpart
        julpartin=juldate(ibdatein,ibtimein)+
     +  dble(float(itramem(i))/86400.)
        nclass(i)=min(int(ran1(idummy)*float(nclassunc))+1,
     +  nclassunc)
        idt(i)=mintime
        itra1(i)=0
        itramem(i)=nint((julpartin-bdate)*86400.)
50      itrasplit(i)=ldirect*itsplit

      return


994   write(*,*) ' #### FLEXPART MODEL ERROR IN READPARTPOSITIONS#### '
      write(*,*) ' #### ENDING TIME OF PREVIOUS MODEL RUN DOES   #### '
      write(*,*) ' #### NOT AGREE WITH STARTING TIME OF THIS RUN.#### '
      write(*,*) 'julin: ',julin
      write(*,*) 'bdate: ',bdate
      stop

995   write(*,*) ' #### FLEXPART MODEL ERROR IN READPARTPOSITIONS#### '
      write(*,*) ' #### NUMBER OF RELEASE LOCATIONS DOES NOT     #### '
      write(*,*) ' #### AGREE WITH CURRENT SETTINGS!             #### '
      stop

996   write(*,*) ' #### FLEXPART MODEL ERROR IN READPARTPOSITIONS#### '
      write(*,*) ' #### SPECIES NAMES TO BE READ IN DO NOT       #### '
      write(*,*) ' #### AGREE WITH CURRENT SETTINGS!             #### '
      stop

997   write(*,*) ' #### FLEXPART MODEL ERROR IN READPARTPOSITIONS#### '
      write(*,*) ' #### THE NUMBER OF SPECIES TO BE READ IN DOES #### '
      write(*,*) ' #### NOT AGREE WITH CURRENT SETTINGS!         #### '
      stop

998   write(*,*) ' #### FLEXPART MODEL ERROR!   THE FILE         #### '
      write(*,*) ' #### '//path(2)(1:len(2))//'grid'//' #### '
      write(*,*) ' #### CANNOT BE OPENED. IF A FILE WITH THIS    #### '
      write(*,*) ' #### NAME ALREADY EXISTS, DELETE IT AND START #### '
      write(*,*) ' #### THE PROGRAM AGAIN.                       #### '
      stop

      end
