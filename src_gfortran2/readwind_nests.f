      subroutine readwind_nests(indj,n,uuhn,vvhn,wwhn)
C                                i   i  o    o    o
********************************************************************************
*                                                                              *
*     Note:  This is the FLEXPART_WRF version of subroutine readwind_nests.    *
*            The met fields are read from WRF netcdf output files.             *
*            There are many differences from the FLEXPART version.             *
*                                                                              *
*     This routine reads the wind fields for the nested model domains.         *
*     It is similar to subroutine readwind, which reads the mother domain.     *
*                                                                              *
*     Authors: A. Stohl, G. Wotawa                                             *
*                                                                              *
*     8 February 1999                                                          *
*     Last update: 17 October 2000, A. Stohl                                   *
*     Changes, Bernd C. Krueger, Feb. 2001:                                    *
*        Variables tthn and qvhn (on eta coordinates) in common block          *
*                                                                              *
*     Oct-Dec 2005, R. Easter -- Major changes for WRF.                        *
*                                                                              *
********************************************************************************


      include 'includepar'
      include 'includecom'

c subr arguments
      integer indj,n
      real uuhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
      real vvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,maxnests)
      real wwhn(0:nxmaxn-1,0:nymaxn-1,nwzmax,maxnests)

c local variables
      integer ndims_max
      parameter (ndims_max = 4)

      integer i, idiagaa, ierr, itime
      integer iduma
      integer j, jhhmmss, jyyyymmdd
      integer k, kbgn
      integer l
      integer lendim(ndims_max), lendim_exp(ndims_max), 
     &    lendim_max(ndims_max)
      integer m
      integer ndims, ndims_exp
      integer n_west_east, n_south_north, n_bottom_top
      integer m_grid_id_dum, m_parent_grid_id_dum, 
     &  m_parent_grid_ratio_dum, 
     &  i_parent_start_dum, j_parent_start_dum,
     &  map_proj_id_dum

      real dx_met, dy_met
      real duma, dumb, dumc, dumd, dume
      real dumdz
      real dumarray_aa(nwzmax+1)
      real dumarray_pp(0:nxmaxn-1,0:nymaxn-1,nwzmax+1)
      real ewater_mb, esatwater_mb
      real ew      ! this is an external function
      real map_stdlon_dum, map_truelat1_dum, map_truelat2_dum
      real toler

      real ewss(0:nxmaxn-1,0:nymaxn-1),nsss(0:nxmaxn-1,0:nymaxn-1)
      real plev1,pmean,tv,fu,hlev1,ff10m,fflev1

      double precision jul
      double precision juldate    ! juldate is a function

      character*160 fnamenc, varname

      logical hflswitch,strswitch


c
c main loop -- process each nest
c
      do 100 l=1,numbnests

c
c   get grid info from the wrf netcdf file
c   and check it for consistency against values from gridcheck
c
      m = numpath+2*(l-1)+1
      fnamenc = path(m)(1:len(m)) // wfnamen(l,indj)

      idiagaa = 0

      call read_ncwrfout_gridinfo( ierr, idiagaa, fnamenc,
     &  n_west_east, n_south_north, n_bottom_top, 
     &  dx_met, dy_met, 
     &  m_grid_id_dum, m_parent_grid_id_dum, m_parent_grid_ratio_dum, 
     &  i_parent_start_dum, j_parent_start_dum,
     &  map_proj_id_dum, map_stdlon_dum, 
     &  map_truelat1_dum, map_truelat2_dum )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error getting gridinfo for met file', 
     &        fnamenc
          stop
      end if

c subtract 1 here because i & j indexing in flexpart always starts at 0
      i_parent_start_dum = i_parent_start_dum-1
      j_parent_start_dum = j_parent_start_dum-1

9100	format( / '*** readwind_nests, l=', i2, ' -- ', 
     &      a / 'file = ', a )
9110	format( / '*** readwind_nests, l=', i2, ' -- ', 
     &      a, 1x, i8 / 'file = ', a )
9120	format( / '*** readwind_nests, l=', i2, ' -- ', 
     &      a, 2(1x,i8) / 'file = ', a )
9130	format( / '*** readwind_nests, l=', i2, ' -- ', 
     &      a, 3(1x,i8) / 'file = ', a )
9115	format( / '*** readwind_nests, l=', i2, ' -- ', 
     &      a / a, 1x, i8 / 'file = ', a )
9125	format( / '*** readwind_nests, l=', i2, ' -- ', 
     &      a / a, 2(1x,i8) / 'file = ', a )
9135	format( / '*** readwind_nests, l=', i2, ' -- ', 
     &      a / a, 3(1x,i8) / 'file = ', a )

      toler = 2.0e-7

      if (nxn(l) .ne. n_west_east) then
          write(*,9100) l, 'nx not consistent', fnamenc
          stop
      end if
      if (nyn(l) .ne. n_south_north) then
          write(*,9100) l, 'ny not consistent', fnamenc
          stop
      end if
      if (nlev_ec .ne. n_bottom_top) then
          write(*,9100) l, 'nlev_ec not consistent', fnamenc
          stop
      end if
      if (nwz .ne. n_bottom_top+1) then
          write(*,9100) l, 'nwz not consistent', fnamenc
          stop
      end if
      if (nuvz .ne. n_bottom_top+add_sfc_level) then
          write(*,9100) l, 'nuvz not consistent', fnamenc
          stop
      end if

      if (m_grid_id(l) .ne. m_grid_id_dum) then
          write(*,9100) l, 'm_grid_id not consistent', fnamenc
          write(*,*) m_grid_id(l), m_grid_id_dum
          stop
      end if
      if (m_parent_grid_id(l) .ne. m_parent_grid_id_dum) then
          write(*,9100) l, 'm_parent_grid_id not consistent', fnamenc
          stop
      end if
      if (m_parent_grid_ratio(l) .ne. m_parent_grid_ratio_dum) then
          write(*,9100) l, 'm_parent_grid_ratio not consistent', fnamenc
          stop
      end if
      if (i_parent_start(l) .ne. i_parent_start_dum) then
          write(*,9100) l, 'i_parent_start not consistent', fnamenc
          stop
      end if
      if (j_parent_start(l) .ne. j_parent_start_dum) then
          write(*,9100) l, 'j_parent_start not consistent', fnamenc
          stop
      end if

      if (abs(dxn(l) - dx_met) .gt. toler*abs(dxn(l))) then
          write(*,9100) l, 'dx not consistent', fnamenc
          stop
      end if
      if (abs(dyn(l) - dy_met) .gt. toler*abs(dyn(l))) then
          write(*,9100) l, 'dy not consistent', fnamenc
          stop
      end if

c locate the date/time in the file
      itime = 0
1100  itime = itime + 1
      call read_ncwrfout_1datetime( ierr, fnamenc,
     &    itime, jyyyymmdd, jhhmmss )
      if (ierr .eq. -1) then
          write(*,9100) l, 'error reading time from met file', fnamenc
          stop
      else if (ierr .ne. 0) then
          write(*,9125) l, 'unable to locate date/time in met file', 
     &        'indj, itime =', indj, itime, fnamenc
          stop
      else 
          jul = juldate( jyyyymmdd, jhhmmss )
          duma = (jul-bdate)*86400.
          iduma = nint(duma)
          if (iduma .ne. wftime(indj)) goto 1100
      end if
      write(*,*) 
      write(*,*) 'readwind_nests processing wrfout file ='
      write(*,*) fnamenc
      write(*,*) 'itime, ymd, hms =', itime, jyyyymmdd, jhhmmss


c read eta_w_wrf, eta_u_wrf, p_top_wrf, ylat2d, xlon2d from the 
c netcdf wrfout file and compare to those from the 1st met. file
      varname = 'ZNW'
      do i = 1, ndims_max
          lendim_exp(i) = 0
          lendim_max(i) = 1
      end do
      lendim_exp(1) = nwz
      lendim_max(1) = nwzmax
      ndims_exp = 2
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, dumarray_aa,
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of ZNW', fnamenc
          stop
      end if
      do k = 1, nwz
          if (abs(eta_w_wrf(k) - dumarray_aa(k)) 
     &            .gt. toler*abs(eta_w_wrf(k))) then
              write(*,9100) l, 'eta_w_wrf not consistent', fnamenc
              stop
          end if
      end do

      varname = 'ZNU'
      lendim_exp(1) = nwz-1
      lendim_max(1) = nwzmax
      ndims_exp = 2
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, dumarray_aa,
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of ZNU', fnamenc
          stop
      end if
      do k = 1, nwz-1
          if (abs(eta_u_wrf(k) - dumarray_aa(k)) 
     &            .gt. toler*abs(eta_u_wrf(k))) then
              write(*,9100) l, 'eta_u_wrf not consistent', fnamenc
              stop
          end if
      end do

      varname = 'P_TOP'
      lendim_exp(1) = 1
      lendim_max(1) = 1
! changed from ndims_exp=2 ; saji
      ndims_exp = 1
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, duma,
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of P_TOP', fnamenc
          stop
      end if
      if (abs(p_top_wrf - duma) .gt. toler*abs(p_top_wrf)) then
          write(*,9100) l, 'p_top_wrf not consistent', fnamenc
          stop
      end if

      varname = 'XLAT'
      lendim_exp(1) = nxn(l)
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)
      lendim_max(2) = nymaxn
      ndims_exp = 3
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, dumarray_pp,
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,*)
          write(*,9100) 'error doing ncread of XLAT', fnamenc
          stop
      end if
      toler = 1.0e-6
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
          if (abs(ylat2dn(i,j,l) - dumarray_pp(i,j,1)) .gt. 
     &                        toler*abs(ylat2dn(i,j,l))) then
              write(*,9100) l, 'ylat2dn not consistent', fnamenc
              write(*,'(a,2i5,2f16.6)') 'i,j,ylats =', i, j,
     &                ylat2dn(i,j,l), dumarray_pp(i,j,1)
              stop
          end if
      end do
      end do

      varname = 'XLONG'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, dumarray_pp,
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,*)
          write(*,9100) 'error doing ncread of XLONG', fnamenc
          stop
      end if
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
          if (abs(xlon2dn(i,j,l) - dumarray_pp(i,j,1)) .gt. 
     &                        toler*abs(xlon2dn(i,j,l))) then
              write(*,9100) l, 'xlon2dn not consistent', fnamenc
              write(*,'(a,2i5,2f16.6)') 'i,j,xlons =', i, j,
     &                xlon2dn(i,j,l), dumarray_pp(i,j,1)
              stop
          end if
      end do
      end do


c
c
c now read the data fields for current time
c the following are read from ecmwf met files
c       U VELOCITY
c       V VELOCITY
c       W VELOCITY
c       TEMPERATURE
c       SPEC. HUMIDITY  
c       SURF. PRESS.
c       SEA LEVEL PRESS.
c       10 M U VELOCITY
c       10 M V VELOCITY
c       2 M TEMPERATURE
c       2 M DEW POINT  
c       SNOW DEPTH
c       CLOUD COVER
c       LARGE SCALE PREC.
c       CONVECTIVE PREC.
c       SENS. HEAT FLUX
c       SOLAR RADIATION
c       EW SURFACE STRESS
c       NS SURFACE STRESS
c       ECMWF OROGRAPHY
c       STANDARD DEVIATION OF OROGRAPHY
c       ECMWF LAND SEA MASK
c
c
      hflswitch=.false.
      strswitch=.false.

      kbgn = 1 + add_sfc_level
c at this point
c   if add_sfc_level=1, then nuvz=nwz   and kbgn=2
c   if add_sfc_level=0, then nuvz=nwz-1 and kbgn=1

c u wind velocity
c   the wrf output file contains (nuvz-add_sfc_level) levels
c   read the data into k=kbgn,nuvz
c   (interpolate it from "U-grid" to "T-grid" later)
      varname = 'U'
      do i = 1, ndims_max
          lendim_exp(i) = 0
          lendim_max(i) = 1
      end do
      lendim_exp(1) = nxn(l)+1
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)
      lendim_max(2) = nymaxn
      lendim_exp(3) = nuvz-add_sfc_level
      lendim_max(3) = nuvzmax
      ndims_exp = 4
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, uuhn(0,0,kbgn,l),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of U', fnamenc
          stop
      end if


c v wind velocity
c   the wrf output file contains (nuvz-add_sfc_level) levels
c   read the data into k=kbgn,nuvz
c   (interpolate it from "V-grid" to "T-grid" later)
      varname = 'V'
      lendim_exp(1) = nxn(l)
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)+1
      lendim_max(2) = nymaxn
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, vvhn(0,0,kbgn,l),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of V', fnamenc
          stop
      end if


c w wind velocity
c   this is on the "W-grid", and 
c   the wrf output file contains nwz levels, so no shifting needed
      varname = 'W'
      lendim_exp(1) = nxn(l)
      lendim_max(1) = nxmaxn
      lendim_exp(2) = nyn(l)
      lendim_max(2) = nymaxn
      lendim_exp(3) = nwz
      lendim_max(3) = nwzmax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, wwhn(0,0,1,l),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of W', fnamenc
          stop
      end if


c pressure - read base state and perturbation pressure,
c     then combine
c   the wrf output file contains (nuvz-add_sfc_level) levels
c   read the data into k=kbgn,nuvz
      varname = 'PB'
      lendim_exp(3) = nuvz-1
      lendim_max(3) = nuvzmax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, pphn(0,0,kbgn,n,l),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of PB', fnamenc
          stop
      end if

      varname = 'P'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, dumarray_pp(0,0,kbgn),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of P', fnamenc
          stop
      end if

      do k = kbgn, nuvz
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
          pphn(i,j,k,n,l) = pphn(i,j,k,n,l) + dumarray_pp(i,j,k)
      end do
      end do
      end do


c height - read base state and perturbation geopotential,
c     then combine and divide by gravity
c   these are on the "W-grid", and 
c     the wrf output file contains nwz levels
c   shift them also so they will be consistent with pph
      varname = 'PHB'
      lendim_exp(3) = nwz
      lendim_max(3) = nwzmax+1
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, zzhn(0,0,kbgn,n,l),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of PB', fnamenc
          stop
      end if

      varname = 'PH'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, dumarray_pp(0,0,kbgn),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of P', fnamenc
          stop
      end if

      do k = kbgn, nwz+add_sfc_level
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
          zzhn(i,j,k,n,l) = 
     &            (zzhn(i,j,k,n,l) + dumarray_pp(i,j,k))/9.81
      end do
      end do
      end do

c now use dumarray_pp to store 1/density for stress calculation below
      if(sfc_option .eq. sfc_option_wrf) then

      varname = 'ALT'
      lendim_exp(3) = nuvz-1
      lendim_max(3) = nwzmax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &    varname, dumarray_pp(0,0,kbgn),
     &    itime,
     &    ndims, ndims_exp, ndims_max,
     &    lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of ALT', fnamenc
          stop
      end if

      end if

c temperature - read perturbation potential temperature,
c     add 300 K (base value), then add and convert
c   the wrf output file contains (nuvz-add_sfc_level) levels
c   read the data into k=kbgn,nuvz
      varname = 'T'
      lendim_exp(3) = nuvz-1
      lendim_max(3) = nuvzmax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, tthn(0,0,kbgn,n,l),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of T', fnamenc
          stop
      end if

      do k = kbgn, nuvz
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
          tthn(i,j,k,n,l) = (tthn(i,j,k,n,l) + 300.0) *
     &            (pphn(i,j,k,n,l)/1.0e5)**0.286
      end do
      end do
      end do


c specific humidity - read mixing ratio (kg-water-vapor/kg-dry-air),
c     then convert to (kg-water-vapor/kg-moist-air)
c   the wrf output file contains (nuvz-add_sfc_level) levels
c   read the data into k=kbgn,nuvz
      varname = 'QVAPOR'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, qvhn(0,0,kbgn,n,l),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of QVAPOR', fnamenc
          stop
      end if

      do k = kbgn, nuvz
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
          qvhn(i,j,k,n,l) = max( qvhn(i,j,k,n,l), 0.0 )
          qvhn(i,j,k,n,l) = qvhn(i,j,k,n,l)/(1.0 + qvhn(i,j,k,n,l))
      end do
      end do
      end do


c surface pressure
      varname = 'PSFC'
      lendim_exp(3) = 0
      lendim_max(3) = 1
      ndims_exp = 3
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, psn(0,0,1,n,l),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of PSFC', fnamenc
          stop
      end if

c for the mexico city grid 3 simulation, the surface and
c   level 1 pressures are not as consistent as one would like,
c   with the problems occuring near the domain boundaries.
c so diagnose surface pressure from other variables
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1

c better fix 
c   -- calculate surface pressure from lowest level pressure, temp, height
c   -- use wrf pressures (pph array) wherever possible
c      (avoid using surface pressure and the akz/bkz, akm/bkm)
          duma = psn(i,j,1,n,l)
          dumdz = 0.5*(zzhn(i,j,kbgn+1,n,l) - zzhn(i,j,kbgn,n,l))
          tv = tthn(i,j,kbgn,n,l)*(1.+0.61*qvhn(i,j,kbgn,n,l))
          psn(i,j,1,n,l) = pphn(i,j,kbgn,n,l)*exp( dumdz*ga/(r_air*tv) )

      end do
      end do


c 10 meter u velocity
c   note:  u10 is on the "T-grid" already
      varname = 'U10'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, u10n(0,0,1,n,l),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of U10', fnamenc
          stop
      end if


c 10 meter v velocity
c   note:  v10 is on the "T-grid" already
      varname = 'V10'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, v10n(0,0,1,n,l),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of V10', fnamenc
          stop
      end if


c 2 meter temperature
      varname = 'T2'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, tt2n(0,0,1,n,l),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of T2', fnamenc
          stop
      end if


c 2 meter dew point - read 2 meter water vapor mixing ratio
c   then calculate the dew point
      varname = 'Q2'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, td2n(0,0,1,n,l),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of Q2', fnamenc
          do j = 0, nyn(l)-1
          do i = 0, nxn(l)-1
c 29-nov-2005 - changed qvhn(i,j,1,n,l) to qvhn(i,j,kbgn,n,l) here
              td2n(i,j,1,n,l) = qvhn(i,j,kbgn,n,l)
          end do
          end do
      end if

c calculate water vapor pressure in mb, from sfc pressure
c   and 2 m mixing ratio
      iduma = 0
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
c 29-nov-2005 - added this to catch occasional tt2n=0.0 values
          duma = max( 100.0, tthn(i,j,kbgn,n,l)-50.0 )
          if (tt2n(i,j,1,n,l) .le. duma) then
              iduma = iduma + 1
              if (iduma .eq. 1) then
                  write(*,*) 'readwind_nests - bad tt2n at'
                  write(*,*) 'l, i, j, tt2n =', l, i, j, tt2n(i,j,1,n,l)
              end if
              tt2n(i,j,1,n,l) = tthn(i,j,kbgn,n,l)
              td2n(i,j,1,n,l) = qvhn(i,j,kbgn,n,l)
          end if
          duma = td2n(i,j,1,n,l)/0.622
          ewater_mb = 0.01*( 0.99976*psn(i,j,1,n,l)*duma/(1.0+duma) )
          esatwater_mb = 0.01*ew(tt2n(i,j,1,n,l))
          ewater_mb = max( 1.0e-10, min( esatwater_mb, ewater_mb ) )
c then use the following, which is from an old 1970's report
c   (reference not available, but the formula works)
c   tdew(in C) = (4318.76/(19.5166 - ln(ewater(in mb)))) - 243.893
          td2n(i,j,1,n,l) = 273.16 +
     +           (4318.76/(19.5166 - log(ewater_mb))) - 243.893
      end do
      end do
      if (iduma .gt. 0) write(*,*)
     &    'readwind_nests - bad tt2n count =', iduma


c sea level pressure - calculate it from surface pressure and 
c    ground elevation using standard atmosphere relations
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
          msln(i,j,1,n,l) = psn(i,j,1,n,l)/
     &            ((1.0 - 6.5e-3*oron(i,j,l)/288.0)**5.2553)
      end do
      end do


c large scale precipitation
c convective  precipitation
c   the wrf output files contain these as "accumulated totals"
c   I need to find out if these are accumulated over the output
c       file frequency, or over the total run.
c   For now, set to zero
c total cloud cover
c   Doesn't appear to be any 2-d cloud cover field in the
c       wrf output.
c   For now, set to zero
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
          lsprecn(i,j,1,n,l) = 0.0
          convprecn(i,j,1,n,l) = 0.0
          tccn(i,j,1,n,l) = 0.0
      end do
      end do


c snow depth
      varname = 'SNOWH'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, sdn(0,0,1,n,l),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of SNOWH', fnamenc
          do j = 0, nyn(l)-1
          do i = 0, nxn(l)-1
              sdn(i,j,1,n,l) = 0.0
          end do
          end do
      end if

      if(sfc_option .eq. sfc_option_wrf) then

c surface sensible heat flux (positive <--> upwards)
      varname = 'HFX'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, sshfn(0,0,1,n,l),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of HFX', fnamenc
          do j = 0, nyn(l)-1
          do i = 0, nxn(l)-1
              sshfn(i,j,1,n,l) = 0.0
          end do
          end do
          hflswitch=.false.    ! Heat flux is not available
      else
          hflswitch=.true.     ! Heat flux is available
c limit to values to bounds originally used by flexpart?
c         do 1502 j=0,nyn(l)-1
c         do 1502 i=0,nxn(l)-1
c            if(sshfn(i,j,1,n,l).gt.200.) sshfn(i,j,1,n,l)=200.
c            if(sshfn(i,j,1,n,l).lt.-400.) sshfn(i,j,1,n,l)=-400.
c1502     continue
      end if

c ustar
      varname = 'UST'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &    varname, ustarn(0,0,1,n,l),
     &    itime,
     &    ndims, ndims_exp, ndims_max,
     &    lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of UST', fnamenc
          do j = 0, nyn(l)
          do i = 0, nxn(l)
              ustarn(i,j,1,n,l) = 0.0
          end do
          end do
          strswitch=.false.    ! ustar is not available
      else
          strswitch=.true.     ! ustar is available
          do 1501 j=0,nyn(l)
          do 1501 i=0,nxn(l)
            surfstrn(i,j,1,n,l)=ustarn(i,j,1,n,l)/dumarray_pp(i,j,kbgn)
1501      continue
      end if

c pblh
      varname = 'PBLH'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &    varname, hmixn(0,0,1,n,l),
     &    itime,
     &    ndims, ndims_exp, ndims_max,
     &    lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of PBLH', fnamenc
          stop
      endif

      endif

c surface solar radiation flux (positive <--> downwards)
      varname = 'SWDOWN'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, ssrn(0,0,1,n,l),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) l, 'error doing ncread of SWDOWN', fnamenc
          do j = 0, nyn(l)-1
          do i = 0, nxn(l)-1
              ssrn(i,j,1,n,l) = 0.0
          end do
          end do
      else
          do j = 0, nyn(l)-1
          do i = 0, nxn(l)-1
              ssrn(i,j,1,n,l) = max( ssrn(i,j,1,n,l), 0.0 )
          end do
          end do
      end if


c ew & ns surface stress
c   Doesn't appear to be any 2-d cloud cover field in the
c       wrf output.
c   For now, set to zero
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
          ewss(i,j) = 0.0
          nsss(i,j) = 0.0
      end do
      end do
c     strswitch=.false.    ! Surface stress is not available


c orography
c standard deviation of orography
c land sea mask
c    these should be fixed during a simulation
c    so there is no reason to do them again ??


c *** done with reading the wrf output file ***



c interpolate uuh from the "U-grid" to the "T-grid"
c interpolate vvh from the "V-grid" to the "T-grid"
      do k = kbgn, nuvz
      do j = 0, nyn(l)-1
      do i = 0, nxn(l)-1
          uuhn(i,j,k,l) = 0.5*(uuhn(i,j,k,l) + uuhn(i+1,j,k,l))
          vvhn(i,j,k,l) = 0.5*(vvhn(i,j,k,l) + vvhn(i,j+1,k,l))
      end do
      end do
      end do



C CALCULATE SURFSTR
      if(sfc_option .eq. sfc_option_diagnosed) then
        do 65 j=0,nyn(l)-1
        do 65 i=0,nxn(l)-1
65        surfstrn(i,j,1,n,l)=sqrt(ewss(i,j)**2+nsss(i,j)**2)
        strswitch=.false.    ! Surface stress is not available
      endif

      if ((.not.hflswitch).or.(.not.strswitch)) then
        write(*,*) 'WARNING: No (or incomplete) flux data ' // 
     +  'contained in WRF output file ',
     +  wfname(indj)
 

C CALCULATE USTAR AND SSHF USING THE PROFILE METHOD
C    As ECMWF has increased the model resolution, such that now the first model
C    level is at about 10 m (where 10-m wind is given), use the 2nd ECMWF level
C    (3rd model level in FLEXPART) for the profile method 
C
C FLEXPART_WRF - use k=(2+add_sfc_level) here instead of k=3
****************************************************************************
        k = 2 + add_sfc_level
        do 1500 j=0,nyn(l)-1
          do 1500 i=0,nxn(l)-1
c           plev1=akz(3)+bkz(3)*psn(i,j,1,n,l)
            plev1=pphn(i,j,k,n,l)
            pmean=0.5*(psn(i,j,1,n,l)+plev1)
            tv=tthn(i,j,k,n,l)*(1.+0.61*qvhn(i,j,k,n,l))
            fu=-r_air*tv/ga/pmean
            hlev1=fu*(plev1-psn(i,j,1,n,l))   ! HEIGTH OF FIRST MODEL LAYER
            ff10m= sqrt(u10n(i,j,1,n,l)**2+v10n(i,j,1,n,l)**2)
            fflev1=sqrt(uuhn(i,j,k,l)**2+vvhn(i,j,k,l)**2)
            call pbl_profile(psn(i,j,1,n,l),td2n(i,j,1,n,l),hlev1,
     &                       tt2n(i,j,1,n,l),tthn(i,j,k,n,l),
     &                       ff10m,fflev1,
     &                       surfstrn(i,j,1,n,l),sshfn(i,j,1,n,l))
            if(sshfn(i,j,1,n,l).gt.200.) sshfn(i,j,1,n,l)=200.
            if(sshfn(i,j,1,n,l).lt.-400.) sshfn(i,j,1,n,l)=-400.
1500        continue
      endif


C Assign 10 m wind to model level at eta=1.0 to have one additional model
C     level at the ground
C Specific humidity is taken the same as at one level above
C Temperature is taken as 2 m temperature         
C
C Note that the uuh, vvh, tth, & qvh data have already been shifted
c     upwards by one level, when they were read in.
***************************************************************************

      if (add_sfc_level .eq. 1) then
      do 70 j = 0, nyn(l)-1
      do 70 i = 0, nxn(l)-1
          uuhn(i,j,1,l)   = u10n(i,j,1,n,l)
          vvhn(i,j,1,l)   = v10n(i,j,1,n,l)
          tthn(i,j,1,n,l) = tt2n(i,j,1,n,l)
          qvhn(i,j,1,n,l) = qvhn(i,j,2,n,l)
c pressure at 2 m AGL
          pphn(i,j,1,n,l) = 0.99976*psn(i,j,1,n,l)
c height (MSL) at ground level (shift it down)
          zzhn(i,j,1,n,l) = zzhn(i,j,2,n,l)
c height (MSL) at top of the added level
          zzhn(i,j,2,n,l) = zzhn(i,j,1,n,l) + 4.0
70    continue
      end if


100   continue


      return    
      end
