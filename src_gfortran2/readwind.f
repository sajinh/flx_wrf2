      subroutine readwind(indj,n,uuh,vvh,wwh)
***********************************************************************
*                                                                     * 
*             TRAJECTORY MODEL SUBROUTINE READWIND                    *
*                                                                     *
***********************************************************************
*                                                                     * 
* AUTHOR:      G. WOTAWA                                              *
* DATE:        1997-08-05                                             *
* LAST UPDATE: 2000-10-17, Andreas Stohl                              *
*                                                                     * 
* Bernd C. Krueger, Feb. 2001:  Variables tth and qvh                 *
*                               (on eta coordinates) in common block  *
*                                                                     * 
* Oct-Dec, 2005: R. Easter.  Major changes for WRF.                   *
*    06-nov-2005 rce - change uuh,vvh dimension back to original      * 
*    16-nov-2005 rce - zzh is shifted like pph,tth                    * 
*                                                                     * 
***********************************************************************
*                                                                     *
* Note:  This is the FLEXPART_WRF version of subroutine readwind.     *
*    The met fields are read from WRF netcdf output files.            *
*    There are many differences from the FLEXPART version.            *
*                                                                     *
* DESCRIPTION:                                                        *
*                                                                     *
* READING OF ECMWF METEOROLOGICAL FIELDS FROM INPUT DATA FILES. THE   *
* INPUT DATA FILES ARE EXPECTED TO BE AVAILABLE IN GRIB CODE          *
*                                                                     *
* INPUT:                                                              *
* indj               indicates number of the wind field to be read in *
* n                  temporal index for meteorological fields (1 to 3)*
*                                                                     *
* IMPORTANT VARIABLES FROM COMMON BLOCK:                              *
*                                                                     *
* wfname             File name of data to be read in                  *
* nx,ny,nuvz,nwz     expected field dimensions                        *
* nlev_ec            number of "T-grid" vertical levels wwf model     *
*                    (the unstaggered "bottom_top" dimension)         *
* uu,vv,ww           wind fields                                      *
* tt,qv              temperature and specific humidity                *
* ps                 surface pressure                                 *
*                                                                     *
***********************************************************************
*

      include 'includepar'
      include 'includecom'

c subr arguments
      integer indj, n

      real uuh(0:nxmax-1,0:nymax-1,nuvzmax)
      real vvh(0:nxmax-1,0:nymax-1,nuvzmax)
      real wwh(0:nxmax-1,0:nymax-1,nwzmax)

c local variables
      integer ndims_max
      parameter (ndims_max = 4)

      integer i, idiagaa, ierr, ifn, itime
      integer iduma
      integer j, jhhmmss, jyyyymmdd
      integer k, kbgn
      integer lendim(ndims_max), lendim_exp(ndims_max), 
     &    lendim_max(ndims_max)
      integer levdiff2
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
      real dumarray_pp(0:nxmax-1,0:nymax-1,nwzmax+1)
      real ewater_mb, esatwater_mb
      real ew      ! this is an external function
      real map_stdlon_dum, map_truelat1_dum, map_truelat2_dum
      real pint
      real toler

      real ewss(0:nxmax-1,0:nymax-1),nsss(0:nxmax-1,0:nymax-1)
      real plev1,pmean,tv,fu,hlev1,ff10m,fflev1

      double precision jul
      double precision juldate    ! juldate is a function

      character*160 fnamenc, varname

      logical hflswitch,strswitch

c
c   get grid info from the wrf netcdf file
c   and check it for consistency against values from gridcheck
c
      fnamenc = path(3)(1:len(3))//wfname(indj)
      idiagaa = 0

      call read_ncwrfout_gridinfo( ierr, idiagaa, fnamenc,
     &  n_west_east, n_south_north, n_bottom_top, 
     &  dx_met, dy_met, 
     &  m_grid_id_dum, m_parent_grid_id_dum, m_parent_grid_ratio_dum, 
     &  i_parent_start_dum, j_parent_start_dum,
     &  map_proj_id_dum, map_stdlon_dum, 
     &  map_truelat1_dum, map_truelat2_dum )
      if (ierr .ne. 0) then
          write(*,9100) 'error getting gridinfor for met file', fnamenc
          stop
      end if

9100	format( / '*** readwind -- ', a /
     &	'file = ', a )
9110	format( / '*** readwind -- ', a, 1x, i8 /
     &	'file = ', a )
9120	format( / '*** readwind -- ', a, 2(1x,i8) /
     &	'file = ', a )
9130	format( / '*** readwind -- ', a, 3(1x,i8) /
     &	'file = ', a )
9115	format( / '*** readwind -- ', a / a, 1x, i8 /
     &	'file = ', a )
9125	format( / '*** readwind -- ', a / a, 2(1x,i8) /
     &	'file = ', a )
9135	format( / '*** readwind -- ', a / a, 3(1x,i8) /
     &	'file = ', a )

      toler = 2.0e-7

      if (nx .ne. n_west_east) then
          write(*,9100) 'nx not consistent', fnamenc
          stop
      end if
      if (ny .ne. n_south_north) then
          write(*,9100) 'ny not consistent', fnamenc
          stop
      end if
      if (nlev_ec .ne. n_bottom_top) then
          write(*,9100) 'nlev_ec not consistent', fnamenc
          stop
      end if
      if (nwz .ne. n_bottom_top+1) then
          write(*,9100) 'nwz not consistent', fnamenc
          stop
      end if
      if (nuvz .ne. n_bottom_top+1) then
          write(*,9100) 'nuvz not consistent', fnamenc
          stop
      end if

      if (m_grid_id(0) .ne. m_grid_id_dum) then
          write(*,9100) 'm_grid_id not consistent', fnamenc
          write(*,*) m_grid_id(0), m_grid_id_dum
          stop
      end if
      if (m_parent_grid_id(0) .ne. m_parent_grid_id_dum) then
          write(*,9100) 'm_parent_grid_id not consistent', fnamenc
          stop
      end if
      if (m_parent_grid_ratio(0) .ne. m_parent_grid_ratio_dum) then
          write(*,9100) 'm_parent_grid_ratio not consistent', fnamenc
          stop
      end if
      if (i_parent_start(0) .ne. i_parent_start_dum) then
          write(*,9100) 'i_parent_start not consistent', fnamenc
          stop
      end if
      if (j_parent_start(0) .ne. j_parent_start_dum) then
          write(*,9100) 'j_parent_start not consistent', fnamenc
          stop
      end if

      if (abs(dx - dx_met) .gt. toler*abs(dx)) then
          write(*,9100) 'dx not consistent', fnamenc
          stop
      end if
      if (abs(dy - dy_met) .gt. toler*abs(dy)) then
          write(*,9100) 'dy not consistent', fnamenc
          stop
      end if

c locate the date/time in the file
      itime = 0
1100  itime = itime + 1
      call read_ncwrfout_1datetime( ierr, fnamenc,
     &    itime, jyyyymmdd, jhhmmss )
      if (ierr .eq. -1) then
          write(*,9100) 'error reading time from met file', fnamenc
          stop
      else if (ierr .ne. 0) then
          write(*,9125) 'unable to locate date/time in met file', 
     &        'indj, itime =', indj, itime, fnamenc
          stop
      else 
          jul = juldate( jyyyymmdd, jhhmmss )
          duma = (jul-bdate)*86400.
          iduma = nint(duma)
          if (iduma .ne. wftime(indj)) goto 1100
      end if
      write(*,*) 
      write(*,*) 'readwind processing wrfout file ='
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
      lendim_max(1) = nwzmax+1
      ndims_exp = 2
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, dumarray_aa,
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of ZNW', fnamenc
          stop
      end if
      do k = 1, nwz
          if (abs(eta_w_wrf(k) - dumarray_aa(k)) 
     &            .gt. toler*abs(eta_w_wrf(k))) then
              write(*,9100) 'eta_w_wrf not consistent', fnamenc
              stop
          end if
      end do

      varname = 'ZNU'
      lendim_exp(1) = nwz-1
      lendim_max(1) = nwzmax+1
      ndims_exp = 2
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, dumarray_aa,
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of ZNU', fnamenc
          stop
      end if
      do k = 1, nwz-1
          if (abs(eta_u_wrf(k) - dumarray_aa(k)) 
     &            .gt. toler*abs(eta_u_wrf(k))) then
              write(*,9100) 'eta_u_wrf not consistent', fnamenc
              stop
          end if
      end do

      varname = 'P_TOP'
      lendim_exp(1) = 1
      lendim_max(1) = 1
! changed from ndims_exp =2 ; saji
      ndims_exp = 1
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, duma,
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of P_TOP', fnamenc
          stop
      end if
      if (abs(p_top_wrf - duma) .gt. toler*abs(p_top_wrf)) then
          write(*,9100) 'p_top_wrf not consistent', fnamenc
          stop
      end if

      varname = 'XLAT'
      lendim_exp(1) = nx
      lendim_max(1) = nxmax
      lendim_exp(2) = ny
      lendim_max(2) = nymax
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
      do j = 0, ny-1
      do i = 0, nx-1
          if (abs(ylat2d(i,j) - dumarray_pp(i,j,1)) .gt. 
     &                        toler*abs(ylat2d(i,j))) then
              write(*,9100) 'ylat2d not consistent', fnamenc
              write(*,'(a,2i5,2f16.6)') 'i,j,ylats =', i, j,
     &                ylat2d(i,j), dumarray_pp(i,j,1)
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
      do j = 0, ny-1
      do i = 0, nx-1
          if (abs(xlon2d(i,j) - dumarray_pp(i,j,1)) .gt. 
     &                        toler*abs(xlon2d(i,j))) then
              write(*,9100) 'xlon2d not consistent', fnamenc
              write(*,'(a,2i5,2f16.6)') 'i,j,xlons =', i, j,
     &                xlon2d(i,j), dumarray_pp(i,j,1)
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
      levdiff2=nlev_ec-nwz+1

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
      lendim_exp(1) = nx+1
      lendim_max(1) = nxmax
      lendim_exp(2) = ny
      lendim_max(2) = nymax
      lendim_exp(3) = nuvz-add_sfc_level
      lendim_max(3) = nuvzmax
      ndims_exp = 4
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, uuh(0,0,kbgn),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of U', fnamenc
          stop
      end if


c v wind velocity
c   the wrf output file contains (nuvz-add_sfc_level) levels
c   read the data into k=kbgn,nuvz
c   (interpolate it from "V-grid" to "T-grid" later)
      varname = 'V'
      lendim_exp(1) = nx
      lendim_max(1) = nxmax
      lendim_exp(2) = ny+1
      lendim_max(2) = nymax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, vvh(0,0,kbgn),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of V', fnamenc
          stop
      end if


c w wind velocity
c   this is on the "W-grid", and 
c   the wrf output file contains nwz levels, so no shifting needed
      varname = 'W'
      lendim_exp(1) = nx
      lendim_max(1) = nxmax
      lendim_exp(2) = ny
      lendim_max(2) = nymax
      lendim_exp(3) = nwz
      lendim_max(3) = nwzmax
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, wwh,
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of W', fnamenc
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
     &	  varname, pph(0,0,kbgn,n),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of PB', fnamenc
          stop
      end if

      varname = 'P'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, dumarray_pp(0,0,kbgn),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of P', fnamenc
          stop
      end if

      do k = kbgn, nuvz
      do j = 0, nymin1
      do i = 0, nxmin1
          pph(i,j,k,n) = pph(i,j,k,n) + dumarray_pp(i,j,k)
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
     &	  varname, zzh(0,0,kbgn,n),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of PB', fnamenc
          stop
      end if

      varname = 'PH'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, dumarray_pp(0,0,kbgn),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of P', fnamenc
          stop
      end if

      do k = kbgn, nwz+add_sfc_level
      do j = 0, nymin1
      do i = 0, nxmin1
          zzh(i,j,k,n) = 
     &            (zzh(i,j,k,n) + dumarray_pp(i,j,k))/9.81
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
          write(*,9100) 'error doing ncread of ALT', fnamenc
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
     &	  varname, tth(0,0,kbgn,n),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of T', fnamenc
          stop
      end if

      do k = kbgn, nuvz
      do j = 0, nymin1
      do i = 0, nxmin1
          tth(i,j,k,n) = (tth(i,j,k,n) + 300.0) *
     &            (pph(i,j,k,n)/1.0e5)**0.286
      end do
      end do
      end do


c specific humidity - read mixing ratio (kg-water-vapor/kg-dry-air),
c     then convert to (kg-water-vapor/kg-moist-air)
c   the wrf output file contains (nuvz-add_sfc_level) levels
c   read the data into k=kbgn,nuvz
      varname = 'QVAPOR'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, qvh(0,0,kbgn,n),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of QVAPOR', fnamenc
          stop
      end if

      do k = kbgn, nuvz
      do j = 0, nymin1
      do i = 0, nxmin1
          qvh(i,j,k,n) = max( qvh(i,j,k,n), 0.0 )
          qvh(i,j,k,n) = qvh(i,j,k,n)/(1.0 + qvh(i,j,k,n))
      end do
      end do
      end do


c surface pressure
      varname = 'PSFC'
      lendim_exp(3) = 0
      lendim_max(3) = 1
      ndims_exp = 3
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, ps(0,0,1,n),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of PSFC', fnamenc
          stop
      end if

c for the mexico city grid 3 simulation, the surface and
c   level 1 pressures are not as consistent as one would like,
c   with the problems occuring near the domain boundaries.
c so do the following
c   -- calculate surface pressure from lowest level pressure, temp, height
c   -- use wrf pressures (pph array) wherever possible
c      (avoid using surface pressure and the akz/bkz, akm/bkm)
      do j = 0, nymin1
      do i = 0, nxmin1
          duma = ps(i,j,1,n)
          dumdz = 0.5*(zzh(i,j,kbgn+1,n) - zzh(i,j,kbgn,n))
          tv = tth(i,j,kbgn,n)*(1.+0.61*qvh(i,j,kbgn,n))
          ps(i,j,1,n) = pph(i,j,kbgn,n)*exp( dumdz*ga/(r_air*tv) )
      end do
      end do


c 10 meter u velocity
c   note:  u10 is on the "T-grid" already
      varname = 'U10'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, u10(0,0,1,n),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of U10', fnamenc
          stop
      end if


c 10 meter v velocity
c   note:  v10 is on the "T-grid" already
      varname = 'V10'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, v10(0,0,1,n),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of V10', fnamenc
          stop
      end if


c 2 meter temperature
      varname = 'T2'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, tt2(0,0,1,n),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of T2', fnamenc
          stop
      end if


c 2 meter dew point - read 2 meter water vapor mixing ratio
c   then calculate the dew point
      varname = 'Q2'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, td2(0,0,1,n),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of Q2', fnamenc
          do j = 0, nymin1
          do i = 0, nxmin1
c 29-nov-2005 - changed qvh(i,j,1,n) to qvh(i,j,kbgn,n) here
              td2(i,j,1,n) = qvh(i,j,kbgn,n)
          end do
          end do
      end if

c calculate water vapor pressure in mb, from sfc pressure
c   and 2 m mixing ratio
      iduma = 0
      do j = 0, nymin1
      do i = 0, nxmin1
c 29-nov-2005 - added this to catch occasional tt2n=0.0 values
          duma = max( 100.0, tth(i,j,kbgn,n)-50.0 )
          if (tt2(i,j,1,n) .le. duma) then
              iduma = iduma + 1
              if (iduma .eq. 1) then
                  write(*,*) 'readwind - bad tt2 at'
                  write(*,*) 'i, j, tt2 =', i, j, tt2(i,j,1,n)
              end if
c             stop
              tt2(i,j,1,n) = tth(i,j,kbgn,n)
              td2(i,j,1,n) = qvh(i,j,kbgn,n)
          end if
          duma = td2(i,j,1,n)/0.622
          ewater_mb = 0.01*( 0.99976*ps(i,j,1,n)*duma/(1.0+duma) )
          esatwater_mb = 0.01*ew(tt2(i,j,1,n))
          ewater_mb = max( 1.0e-10, min( esatwater_mb, ewater_mb ) )
c then use the following, which is from an old 1970's report
c   (reference not available, but the formula works)
c   tdew(in C) = (4318.76/(19.5166 - ln(ewater(in mb)))) - 243.893
          td2(i,j,1,n) = 273.16 +
     +           (4318.76/(19.5166 - log(ewater_mb))) - 243.893
      end do
      end do
      if (iduma .gt. 0) write(*,*)
     &    'readwind - bad tt2 count =', iduma


c sea level pressure - calculate it from surface pressure and 
c    ground elevation using standard atmosphere relations
      do j = 0, nymin1
      do i = 0, nxmin1
          msl(i,j,1,n) = ps(i,j,1,n)/
     &            ((1.0 - 6.5e-3*oro(i,j)/288.0)**5.2553)
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
      do j = 0, nymin1
      do i = 0, nxmin1
          lsprec(i,j,1,n) = 0.0
          convprec(i,j,1,n) = 0.0
          tcc(i,j,1,n) = 0.0
      end do
      end do


c snow depth
      varname = 'SNOWH'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, sd(0,0,1,n),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of SNOWH', fnamenc
          do j = 0, nymin1
          do i = 0, nxmin1
              sd(i,j,1,n) = 0.0
          end do
          end do
      end if

      if(sfc_option .eq. sfc_option_wrf) then

c surface sensible heat flux (positive <--> upwards)
      varname = 'HFX'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, sshf(0,0,1,n),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of HFX', fnamenc
          do j = 0, nymin1
          do i = 0, nxmin1
              sshf(i,j,1,n) = 0.0
          end do
          end do
          hflswitch=.false.    ! Heat flux is not available
      else
          hflswitch=.true.     ! Heat flux is available
c limit to values to bounds originally used by flexpart?
c        do 1502 j=0,nymin1
c        do 1502 i=0,nxmin1
c           if(sshf(i,j,1,n).gt.200.) sshf(i,j,1,n)=200.
c           if(sshf(i,j,1,n).lt.-400.) sshf(i,j,1,n)=-400.
c1502    continue
      end if

c ustar
      varname = 'UST'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &      varname, ustar(0,0,1,n),
     &      itime,
     &      ndims, ndims_exp, ndims_max,
     &      lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
         write(*,9100) 'error doing ncread of UST', fnamenc
         do j = 0, nymin1
         do i = 0, nxmin1
             ustar(i,j,1,n) = 0.0
         end do
         end do
         strswitch=.false.    ! ustar is not available
      else
         strswitch=.true.     ! ustar is available
         do 1501 j=0,nymin1
         do 1501 i=0,nxmin1
            surfstr(i,j,1,n)=ustar(i,j,1,n)/dumarray_pp(i,j,kbgn)
1501     continue
      end if

c pblh
      varname = 'PBLH'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &    varname, hmix(0,0,1,n),
     &    itime,
     &    ndims, ndims_exp, ndims_max,
     &    lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of PBLH', fnamenc
          stop
      endif

      endif

c surface solar radiation flux (positive <--> downwards)
      varname = 'SWDOWN'
      call read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, ssr(0,0,1,n),
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
      if (ierr .ne. 0) then
          write(*,9100) 'error doing ncread of SWDOWN', fnamenc
          do j = 0, nymin1
          do i = 0, nxmin1
              ssr(i,j,1,n) = 0.0
          end do
          end do
      else
          do j = 0, nymin1
          do i = 0, nxmin1
              ssr(i,j,1,n) = max( ssr(i,j,1,n), 0.0 )
          end do
          end do
      end if


c ew & ns surface stress
c   Doesn't appear to be any 2-d cloud cover field in the
c       wrf output.
c   For now, set to zero
      do j = 0, nymin1
      do i = 0, nxmin1
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
      do j = 0, nymin1
      do i = 0, nxmin1
          uuh(i,j,k) = 0.5*(uuh(i,j,k) + uuh(i+1,j,k))
          vvh(i,j,k) = 0.5*(vvh(i,j,k) + vvh(i,j+1,k))
      end do
      end do
      end do



c for ecmwf flexpart, if nwz = nlev_ec+1, then wwh is set
c   to zero at the top level
c for wrf, nlev_ec==n_bottom_top and so nwz = nlev_ec+1.
c   however, it doesn't seem appropriate to zero wwh at 
c   the model top which might be ~100 hPa.
c so deactivate this for now
c      if(levdiff2.eq.0) then
c        iwmax=nlev_ec+1
c        do 60 i=0,nxmin1
c        do 60 j=0,nymin1
c60      wwh(i,j,nlev_ec+1)=0.
c      endif

C For global fields, assign the leftmost data column also to the rightmost
C data column; if required, shift whole grid by nxshift grid points
C
C FLEXPART_WRF - all "global" stuff is turned off
**************************************************************************

c     if (xglobal) then
c       call shift_field_0(ewss,nxfield,ny)
c       call shift_field_0(nsss,nxfield,ny)
c       call shift_field_0(oro,nxfield,ny)
c       call shift_field_0(excessoro,nxfield,ny)
c       call shift_field_0(lsm,nxfield,ny)
c       call shift_field(ps,nxfield,ny,1,1,2,n)
c       call shift_field(sd,nxfield,ny,1,1,2,n)
c       call shift_field(msl,nxfield,ny,1,1,2,n)
c       call shift_field(tcc,nxfield,ny,1,1,2,n)
c       call shift_field(u10,nxfield,ny,1,1,2,n)
c       call shift_field(v10,nxfield,ny,1,1,2,n)
c       call shift_field(tt2,nxfield,ny,1,1,2,n)
c       call shift_field(td2,nxfield,ny,1,1,2,n)
c       call shift_field(lsprec,nxfield,ny,1,1,2,n)
c       call shift_field(convprec,nxfield,ny,1,1,2,n)
c       call shift_field(sshf,nxfield,ny,1,1,2,n)
c       call shift_field(ssr,nxfield,ny,1,1,2,n)
c       call shift_field(tth,nxfield,ny,nuvzmax,nuvz,2,n)
c       call shift_field(qvh,nxfield,ny,nuvzmax,nuvz,2,n)
c       call shift_field(uuh,nxfield,ny,nuvzmax,nuvz,1,1)
c       call shift_field(vvh,nxfield,ny,nuvzmax,nuvz,1,1)
c       call shift_field(wwh,nxfield,ny,nwzmax,nwz,1,1)
c     endif

C CALCULATE SURFSTR
      if(sfc_option .eq. sfc_option_diagnosed) then
        do 65 i=0,nxmin1
        do 65 j=0,nymin1
65        surfstr(i,j,1,n)=sqrt(ewss(i,j)**2+nsss(i,j)**2)
        strswitch=.false.
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
        do 1500 i=0,nxmin1
          do 1500 j=0,nymin1
c           plev1=akz(3)+bkz(3)*ps(i,j,1,n)
            plev1=pph(i,j,k,n)
            pmean=0.5*(ps(i,j,1,n)+plev1)
            tv=tth(i,j,k,n)*(1.+0.61*qvh(i,j,k,n))
            fu=-r_air*tv/ga/pmean
            hlev1=fu*(plev1-ps(i,j,1,n))   ! HEIGTH OF FIRST MODEL LAYER
            ff10m= sqrt(u10(i,j,1,n)**2+v10(i,j,1,n)**2)
            fflev1=sqrt(uuh(i,j,k)**2+vvh(i,j,k)**2)
            call pbl_profile(ps(i,j,1,n),td2(i,j,1,n),hlev1,
     &                       tt2(i,j,1,n),tth(i,j,k,n),ff10m,fflev1,
     &                       surfstr(i,j,1,n),sshf(i,j,1,n))
            if(sshf(i,j,1,n).gt.200.) sshf(i,j,1,n)=200.
            if(sshf(i,j,1,n).lt.-400.) sshf(i,j,1,n)=-400.
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
      do 70 j = 0, nymin1
      do 70 i = 0, nxmin1
          uuh(i,j,1)   = u10(i,j,1,n)
          vvh(i,j,1)   = v10(i,j,1,n)
          tth(i,j,1,n) = tt2(i,j,1,n)
          qvh(i,j,1,n) = qvh(i,j,2,n)
c pressure at 2 m AGL
          pph(i,j,1,n) = 0.99976*ps(i,j,1,n)
c height (MSL) at ground level (shift it down)
          zzh(i,j,1,n) = zzh(i,j,2,n)
c height (MSL) at top of the added level
          zzh(i,j,2,n) = zzh(i,j,1,n) + 4.0
70    continue
      end if


      return    
      end

