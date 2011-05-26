! file map_proj_wrf_subaa.f - created 20-nov-2005
!
!   this file contains the "ps" (polar stereographic) 
!       and "lc" (lambert conformal) portions
!	file .../wrfsi2.1/src/mod/module_map_utils.F
!   the routines were converted from fortran90 to fortran77 as follows
!	the module variables and parameters were put into the 
!	    proj_info_cmn01 & _cmn02 common blocks
!	    (parameters were changed to variables)
!	the variables in the "proj_info" structure were renamed
!	    so that "rebydx" became "proj_rebydx", etc.
!	    the fortran90 usages of these variables
!	    changed from "proj%rebydx" to "proj_rebydx", etc.
!-----------------------------------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!dis   
!dis    open source license/disclaimer, forecast systems laboratory
!dis    noaa/oar/fsl, 325 broadway boulder, co 80305
!dis    
!dis    this software is distributed under the open source definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis    
!dis    in particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis    
!dis    - redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis    
!dis    - redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis    
!dis    - all modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis    
!dis    - if significant modifications or enhancements are made to this
!dis    software, the fsl software policy manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis    
!dis    this software and its documentation are in the public domain
!dis    and are furnished "as is."  the authors, the united states
!dis    government, its instrumentalities, officers, employees, and
!dis    agents make no warranty, express or implied, as to the usefulness
!dis    of the software and documentation for any purpose.  they assume
!dis    no responsibility (1) for the use of the software and
!dis    documentation; or (2) to provide technical support to users.
!dis   
!dis 

! module that defines constants, data structures, and
! routines used to convert grid indices to lat/lon
! and vice versa.   
!
! supported projections
! ---------------------
! cylindrical lat/lon (code = proj_latlon)
! mercator (code = proj_merc)
! lambert conformal (code = proj_lc)
! polar stereographic (code = proj_ps)
!
! remarks
! -------
! the routines contained within were adapted from routines
! obtained from the ncep w3 library.  the original ncep routines were less
! flexible (e.g., polar-stereo routines only supported truelat of 60n/60s)
! than what we needed, so modifications based on equations in hoke, hayes, and
! renninger (afgwc/tn/79-003) were added to improve the flexibility.  
! additionally, coding was improved to f90 standards and the routines were
! combined into this module.  
!
! assumptions
! -----------
!  grid definition:
!    for mercator, lambert conformal, and polar-stereographic projections,
!    the routines within assume the following:
!
!       1.  grid is dimensioned (i,j) where i is the east-west direction, 
!           positive toward the east, and j is the north-south direction, 
!           positive toward the north.  
!       2.  origin is at (1,1) and is located at the southwest corner,
!           regardless of hemispere.
!       3.  grid spacing (dx) is always positive.
!       4.  values of true latitudes must be positive for nh domains
!           and negative for sh domains.
!
!     for the latlon projection, the grid origin may be at any of the
!     corners, and the deltalat and deltalon values can be signed to 
!     account for this using the following convention:
!       origin location        deltalat sign      deltalon sign
!       ---------------        -------------      -------------
!        sw corner                  +                   +
!        ne corner                  -                   -
!        nw corner                  -                   +
!        se corner                  +                   -
!       
!  data definitions:
!       1. any arguments that are a latitude value are expressed in 
!          degrees north with a valid range of -90 -> 90
!       2. any arguments that are a longitude value are expressed in
!          degrees east with a valid range of -180 -> 180.
!       3. distances are in meters and are always positive.
!       4. the standard longitude (stdlon) is defined as the longitude
!          line which is parallel to the y-axis (j-direction), along
!          which latitude increases (not the absolute value of latitude, but
!          the actual latitude, such that latitude increases continuously
!          from the south pole to the north pole) as j increases.  
!       5. one true latitude value is required for polar-stereographic and
!          mercator projections, and defines at which latitude the 
!          grid spacing is true.  for lambert conformal, two true latitude
!          values must be specified, but may be set equal to each other to
!          specify a tangent projection instead of a secant projection.  
!       
! usage
! -----
! to use the routines in this module, the calling routines must have the 
! following statement at the beginning of its declaration block:
!   use map_utils
! 
! the use of the module not only provides access to the necessary routines,
! but also defines a structure of type (proj_info) that can be used
! to declare a variable of the same type to hold your map projection
! information.  it also defines some integer parameters that contain
! the projection codes so one only has to use those variable names rather
! than remembering the acutal code when using them.  the basic steps are
! as follows:
!  
!   1.  ensure the "use map_utils" is in your declarations.
!   2.  declare the projection information structure as type(proj_info):
!         type(proj_info) :: proj
!   3.  populate your structure by calling the map_set routine:
!         call map_set(code,lat1,lon1,dx,stdlon,truelat1,truelat2,nx,ny,proj)
!       where:
!         code (input) = one of proj_latlon, proj_merc, proj_lc, or proj_ps
!         lat1 (input) = latitude of grid origin point (i,j)=(1,1) 
!                         (see assumptions!)
!         lon1 (input) = longitude of grid origin 
!         dx (input) = grid spacing in meters (ignored for latlon projections)
!         stdlon (input) = standard longitude for proj_ps and proj_lc, 
!               deltalon (see assumptions) for proj_latlon, 
!               ignored for proj_merc
!         truelat1 (input) = 1st true latitude for proj_ps, proj_lc, and
!                proj_merc, deltalat (see assumptions) for proj_latlon
!         truelat2 (input) = 2nd true latitude for proj_lc, 
!                ignored for all others.
!         nx = number of points in east-west direction
!         ny = number of points in north-south direction
!         proj (output) = the structure of type (proj_info) that will be fully 
!                populated after this call
!
!   4.  now that the proj structure is populated, you may call any 
!       of the following routines:
!       
!       latlon_to_ij(proj, lat, lon, i, j)
!       ij_to_latlon(proj, i, j, lat, lon)
!       truewind_to_gridwind(lon, proj, ugrid, vgrid, utrue, vtrue)
!       gridwind_to_truewind(lon, proj, utrue, vtrue, ugrid, vgrid)
!       compare_projections(proj1, proj2, same_proj)
!
!       it is incumbent upon the calling routine to determine whether or
!       not the values returned are within your domain bounds.  all values
!       of i, j, lat, and lon are real values.
!
!
! references
! ----------
!  hoke, hayes, and renninger, "map preojections and grid systems for
!       meteorological applications." afgwc/tn-79/003(rev), air weather
!       service, 1985.
!
!  ncar mm5v3 modeling system, regridder program, module_first_guess_map.f
!  ncep routines w3fb06, w3fb07, w3fb08, w3fb09, w3fb11, w3fb12
!
! history
! -------
! 27 mar 2001 - original version
!               brent l. shaw, noaa/fsl (csu/cira)
! 02 apr 2001 - added routines to rotate winds from true to grid
!               and vice versa.
!               brent l. shaw, noaa/fsl (csu/cira)
! 09 apr 2001 - added compare_projections routine to compare two
!               sets of projection parameters.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c ! define data structures to define various projections
c
c type proj_info
c
c   logical        proj_init     ! flag to indicate if this struct is 
c                                ! ready for use
c   logical        proj_cyclic   ! flag indicating if this grid
c                                ! is cyclic in the longitudinal
c                                ! direction...happens with
c                                ! global lat/lon grids like gfs/avn
c   integer        proj_code     ! integer code for projection type
c   integer        proj_nx
c   integer        proj_ny
c   real           proj_lat1    ! sw latitude (1,1) in degrees (-90->90n)
c   real           proj_lon1    ! sw longitude (1,1) in degrees (-180->180e)
c   real           proj_dx       ! grid spacing in meters at truelats, used
c                                ! only for ps, lc, and merc projections
c   real           proj_dlat     ! lat increment for lat/lon grids
c   real           proj_dlon     ! lon increment for lat/lon grids
c   real           proj_clat     ! center latitude of grid
c   real           proj_clon     ! center longitude of grid
c   real           proj_stdlon   ! longitude parallel to y-axis (-180->180e)
c   real           proj_truelat1 ! first true latitude (all projections)
c   real           proj_truelat2 ! second true lat (lc only)
c   real           proj_hemi     ! 1 for nh, -1 for sh
c   real           proj_cone     ! cone factor for lc projections
c   real           proj_polei    ! computed i-location of pole point
c   real           proj_polej    ! computed j-location of pole point
c   real           proj_rsw      ! computed radius to sw corner
c   real           proj_rebydx   ! earth radius divided by dx
c
c end type proj_info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine map_init
    ! initializes the map projection structure to missing values

      implicit none

      include 'include_wrf_map_utils'

      pi = 3.1415927
      deg_per_rad = 180./pi
      rad_per_deg = pi / 180.

! mean earth radius in m.  the value below is consistent
! with nceps routines and grids.
      earth_radius_m = 6370000.

      proj_latlon = 0
      proj_merc = 1
      proj_lc = 3
      proj_ps = 5
      proj_rotlat = 203

      proj_lat1 =    -999.9
      proj_lon1 =    -999.9
      proj_dx    =    -999.9
      proj_stdlon =   -999.9
      proj_truelat1 = -999.9
      proj_truelat2 = -999.9
      proj_hemi     = 0.0
      proj_cone     = -999.9
      proj_polei    = -999.9
      proj_polej    = -999.9
      proj_rsw      = -999.9
      proj_init     = .false.
      proj_nx       = -99
      proj_ny       = -99 
      proj_cyclic   = .false.

      return
      end 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine map_set(
     &    proj_code_in, lat1, lon1, dx,
     &    stdlon, truelat1, truelat2,
     &    idim, jdim )
      ! given a partially filled proj_info structure, this routine computes
      ! polei, polej, rsw, and cone (if lc projection) to complete the 
      ! structure.  this allows us to eliminate redundant calculations when
      ! calling the coordinate conversion routines multiple times for the
      ! same map.
      ! this will generally be the first routine called when a user wants
      ! to be able to use the coordinate conversion routines, and it
      ! will call the appropriate routines based on the 
      ! proj_code which indicates which projection type  this is.

      implicit none
      
      ! declare arguments
      integer proj_code_in
      real lat1
      real lon1
      real dx
      real stdlon
      real truelat1
      real truelat2
      integer idim
      integer jdim

      ! local variables
      real center_i,center_j
      real center_lat, center_lon

      include 'include_wrf_map_utils'

      ! executable code

      proj_code = proj_code_in

      ! first, check for validity of mandatory variables in proj
      if ( abs(lat1) .gt. 90.001 ) then
        print '(a)', 'latitude of origin corner required as follows:'
        print '(a)', '    -90n <= lat1 < = 90.n'
        stop 'map_init'
      endif
      if ( abs(lon1) .gt. 180.) then
        print '(a)', 'longitude of origin required as follows:'
        print '(a)', '   -180e <= lon1 <= 180w'
        stop 'map_init'
      endif
      if ((dx .le. 0.).and.(proj_code .ne. proj_latlon)) then
        print '(a)', 'require grid spacing (dx) in meters be positive!'
        stop 'map_init'
      endif
      if ((abs(stdlon) .gt. 180.).and.(proj_code .ne. proj_merc)) then
        print '(a)', 'need orientation longitude (stdlon) as: '
        print '(a)', '   -180e <= lon1 <= 180w' 
        stop 'map_init'
      endif
      if (abs(truelat1).gt.90.) then
        print '(a)', 'set true latitude 1 for all projections!'
        stop 'map_init'
      endif
     
      call map_init       
      proj_code  = proj_code_in
      proj_lat1 = lat1
      proj_lon1 = lon1
      proj_dx    = dx
      proj_stdlon = stdlon
      proj_truelat1 = truelat1
      proj_truelat2 = truelat2
      proj_nx = idim
      proj_ny = jdim
      if (proj_code .ne. proj_latlon) then
        proj_dx = dx
        if (truelat1 .lt. 0.) then
          proj_hemi = -1.0 
        else
          proj_hemi = 1.0
        endif
        proj_rebydx = earth_radius_m / dx
      endif


      if (proj_code .eq. proj_ps) then
          !print '(a)', 'setting up polar stereographic map...'
          call set_ps

      else if (proj_code .eq. proj_lc) then
          !print '(a)', 'setting up lambert conformal map...'
          if (abs(proj_truelat2) .gt. 90.) then
            print '(a)', 
     &          'second true latitude not set, assuming a tangent'
            print '(a,f10.3)', 'projection at truelat1: ', proj_truelat1
            proj_truelat2=proj_truelat1
          else 
            ! ensure truelat1 < truelat2
            proj_truelat1 = min(truelat1,truelat2)
            proj_truelat2 = max(truelat1,truelat2)
          endif
          call set_lc
     
c     else if (proj_code .eq. proj_merc) then
c         !print '(a)', 'setting up mercator map...'
c         call set_merc
c    
c     else if (proj_code .eq. proj_latlon) then
c         !print '(a)', 'setting up cylindrical equidistant latlon map...'
c         ! convert lon1 to 0->360 notation
c         if (proj_lon1 .lt. 0.) proj_lon1 = proj_lon1 + 360.
c         proj_dlat = truelat1
c         proj_dlon = stdlon 
c         if (nint(proj_dlon*float(proj_nx)) .eq. 360) proj_cyclic = .true.

      else
          print '(a,i2)', 'unknown projection code: ', proj_code
          stop 'map_init'
      
      end if

      proj_init = .true.

      center_i = float(proj_nx+1)*0.5
      center_j = float(proj_ny+1)*0.5
      call ij_to_latlon(center_i,center_j,center_lat,center_lon)
      proj_clat = center_lat
      proj_clon = center_lon

      return
      end 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine latlon_to_ij(lat, lon, i, j)
      ! converts input lat/lon values to the cartesian (i,j) value
      ! for the given projection. 

      implicit none

      ! arguments
      real lat
      real lon
      real i
      real j

      include 'include_wrf_map_utils'

      if (.not.proj_init) then
        print '(a)', 'you have not called map_set for this projection!'
        stop 'latlon_to_ij'
      endif

      if (proj_code .eq. proj_ps) then
          call llij_ps(lat,lon,i,j)
        
      else if (proj_code .eq. proj_lc) then
          call llij_lc(lat,lon,i,j)
 
c     else if (proj_code .eq. proj_latlon) then
c         call llij_latlon(lat,lon,i,j)
c
c     else if (proj_code .eq. proj_merc) then
c         call llij_merc(lat,lon,i,j)
c
c     else if (proj_code .eq. proj_rotlat) then
c	write(6,*)   'doing nothing in latlon_to_ij'

      else 
          print '(a,i2)','unrecognized map projection code: ', proj_code
          stop 'latlon_to_ij'
   
      end if

      return
      end 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ij_to_latlon(i, j, lat, lon)
      ! computes geographical latitude and longitude for a given (i,j) point
      ! in a grid with a projection of proj

      implicit none

      ! arguments
      real i
      real j
      real lat
      real lon

      include 'include_wrf_map_utils'

      if (.not.proj_init) then
        print '(a)', 'you have not called map_set for this projection!'
        stop 'ij_to_latlon'
      endif

      if (proj_code .eq. proj_ps) then
          call ijll_ps(i, j, lat, lon)

      else if (proj_code .eq. proj_lc) then
          call ijll_lc(i, j, lat, lon)
 
c     else if (proj_code .eq. proj_latlon) then
c         call ijll_latlon(i, j, lat, lon)
c
c     else if (proj_code .eq. proj_merc) then
c         call ijll_merc(i, j, lat, lon)
c
c     else if (proj_code .eq. proj_rotlat) then
c	write(6,*)   'doing nothing in ij_to_latlon'

      else
          print '(a,i2)','unrecognized map projection code: ', proj_code
          stop 'ij_to_latlon'

      end if

      return
      end 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine set_ps
      ! initializes a polar-stereographic map projection from the partially
      ! filled proj structure. this routine computes the radius to the
      ! southwest corner and computes the i/j location of the pole for use
      ! in llij_ps and ijll_ps.

      implicit none
   
      ! local vars
      real ala1
      real alo1
      real reflon
      real scale_top
      real dumlat, dumlon

      include 'include_wrf_map_utils'

      ! executable code
      reflon = proj_stdlon + 90.
    
      ! cone factor
      proj_cone = 1.0

      ! compute numerator term of map scale factor
      scale_top = 1. + proj_hemi * sin(proj_truelat1 * rad_per_deg)

      ! compute radius to lower-left (sw) corner
      ala1 = proj_lat1 * rad_per_deg
      proj_rsw =proj_rebydx*cos(ala1)*scale_top/(1.+proj_hemi*sin(ala1))

      ! find the pole point
      alo1 = (proj_lon1 - reflon) * rad_per_deg
      proj_polei = 1. - proj_rsw * cos(alo1)
      proj_polej = 1. - proj_hemi * proj_rsw * sin(alo1)
!     print '(a,2f14.5)', 'set_ps - computed (i,j) of pole point: ',
!    &    proj_polei,proj_polej

      return
      end 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine llij_ps( lat, lon, i, j )
      ! given latitude (-90 to 90), longitude (-180 to 180), and the
      ! standard polar-stereographic projection information via the 
      ! public proj structure, this routine returns the i/j indices which
      ! if within the domain range from 1->nx and 1->ny, respectively.

      implicit none

      ! delcare input arguments
      real lat
      real lon

      ! declare output arguments     
      real i !(x-index)
      real j !(y-index)

      ! declare local variables
      real reflon
      real scale_top
      real ala
      real alo
      real rm

      include 'include_wrf_map_utils'

      ! begin code
    
      reflon = proj_stdlon + 90.
     
      ! compute numerator term of map scale factor

      scale_top = 1. + proj_hemi * sin(proj_truelat1 * rad_per_deg)

      ! find radius to desired point
      ala = lat * rad_per_deg
      rm = proj_rebydx * cos(ala) * scale_top/(1. + proj_hemi *sin(ala))
      alo = (lon - reflon) * rad_per_deg
      i = proj_polei + rm * cos(alo)
      j = proj_polej + proj_hemi * rm * sin(alo)
   
      return
      end 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ijll_ps( i, j, lat, lon )

      ! this is the inverse routine of llij_ps.  it returns the 
      ! latitude and longitude of an i/j point given the projection info 
      ! structure.  

      implicit none

      ! declare input arguments
      real i    ! column
      real j    ! row
      
      ! declare output arguments
      real lat     ! -90 -> 90 north
      real lon     ! -180 -> 180 east

      ! local variables
      real reflon
      real scale_top
      real xx,yy
      real gi2, r2
      real arccos

      include 'include_wrf_map_utils'

      ! begin code

      ! compute the reference longitude by rotating 90 degrees to the east
      ! to find the longitude line parallel to the positive x-axis.
      reflon = proj_stdlon + 90.
     
      ! compute numerator term of map scale factor
      scale_top = 1. + proj_hemi * sin(proj_truelat1 * rad_per_deg)

      ! compute radius to point of interest
      xx = i - proj_polei
      yy = (j - proj_polej) * proj_hemi
      r2 = xx**2 + yy**2

      ! now the magic code
      if (r2 .eq. 0.) then 
        lat = proj_hemi * 90.
        lon = reflon
      else
        gi2 = (proj_rebydx * scale_top)**2.
        lat = deg_per_rad * proj_hemi * asin((gi2-r2)/(gi2+r2))
        arccos = acos(xx/sqrt(r2))
        if (yy .gt. 0) then
          lon = reflon + deg_per_rad * arccos
        else
          lon = reflon - deg_per_rad * arccos
        endif
      endif
    
      ! convert to a -180 -> 180 east convention
      if (lon .gt. 180.) lon = lon - 360.
      if (lon .lt. -180.) lon = lon + 360.

      return
      end 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine set_lc      
      ! initialize the remaining items in the proj structure for a
      ! lambert conformal grid.

      implicit none

      include 'include_wrf_map_utils'
      
      real arg
      real deltalon1
      real tl1r
      real ctl1r

      ! compute cone factor
      call lc_cone(proj_truelat1, proj_truelat2, proj_cone)
      ! print '(a,f8.6)', 'computed cone factor: ', proj_cone
      ! compute longitude differences and ensure we stay out of the
      ! forbidden "cut zone"
      deltalon1 = proj_lon1 - proj_stdlon
      if (deltalon1 .gt. +180.) deltalon1 = deltalon1 - 360.
      if (deltalon1 .lt. -180.) deltalon1 = deltalon1 + 360.

      ! convert truelat1 to radian and compute cos for later use
      tl1r = proj_truelat1 * rad_per_deg
      ctl1r = cos(tl1r)

      ! compute the radius to our known lower-left (sw) corner
      proj_rsw = proj_rebydx * ctl1r/proj_cone * 
     &       (tan((90.*proj_hemi-proj_lat1)*rad_per_deg/2.) / 
     &        tan((90.*proj_hemi-proj_truelat1)*rad_per_deg/2.))
     &        **proj_cone

      ! find pole point
      arg = proj_cone*(deltalon1*rad_per_deg)
      proj_polei = 1. - proj_hemi * proj_rsw * sin(arg)
      proj_polej = 1. + proj_rsw * cos(arg)  
      !print '(a,2f10.3)', 'computed pole i/j = ', proj_polei, proj_polej

      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine lc_cone(truelat1, truelat2, cone)

    ! routine to compute the cone factor of a lambert conformal projection

      implicit none

      include 'include_wrf_map_utils'
      
      ! input args
      real truelat1  ! (-90 -> 90 degrees n)
      real truelat2  !   "   "  "   "     "

      ! output args
      real cone

      ! locals

      ! begin code

      ! first, see if this is a secant or tangent projection.  for tangent
      ! projections, truelat1 = truelat2 and the cone is tangent to the 
      ! earth surface at this latitude.  for secant projections, the cone
      ! intersects the earth surface at each of the distinctly different
      ! latitudes
      if (abs(truelat1-truelat2) .gt. 0.1) then

        ! compute cone factor following:
        cone=(alog(cos(truelat1*rad_per_deg))
     &       -alog(cos(truelat2*rad_per_deg))) / 
     &   (alog(tan((90.-abs(truelat1))*rad_per_deg*0.5 ))- 
     &    alog(tan((90.-abs(truelat2))*rad_per_deg*0.5 )) )
      else
         cone = sin(abs(truelat1)*rad_per_deg )  
      endif

      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ijll_lc( i, j, lat, lon)

    ! routine to convert from the (i,j) cartesian coordinate to the 
    ! geographical latitude and longitude for a lambert conformal projection.

    ! history:
    ! 25 jul 01: corrected by b. shaw, noaa/fsl
    ! 
      implicit none

      include 'include_wrf_map_utils'

      ! input args
      real i        ! cartesian x coordinate
      real j        ! cartesian y coordinate

      ! output args                 
      real lat      ! latitude (-90->90 deg n)
      real lon      ! longitude (-180->180 e)

      ! locals 
      real inew
      real jnew
      real r
      real chi,chi1,chi2
      real r2
      real xx
      real yy

      ! begin code

      chi1 = (90. - proj_hemi*proj_truelat1)*rad_per_deg
      chi2 = (90. - proj_hemi*proj_truelat2)*rad_per_deg

      ! see if we are in the southern hemispere and flip the indices
      ! if we are. 
      if (proj_hemi .eq. -1.) then 
        inew = -i + 2.
        jnew = -j + 2.
      else
        inew = i
        jnew = j
      endif

      ! compute radius**2 to i/j location
      xx = inew - proj_polei
      yy = proj_polej - jnew
      r2 = (xx*xx + yy*yy)
      r = sqrt(r2)/proj_rebydx
     
      ! convert to lat/lon
      if (r2 .eq. 0.) then
        lat = proj_hemi * 90.
        lon = proj_stdlon
      else
         
        ! longitude
        lon = proj_stdlon + deg_per_rad * atan2(proj_hemi*xx,yy)
     &        /proj_cone
        lon = amod(lon+360., 360.)

        ! latitude.  latitude determined by solving an equation adapted 
        ! from:
        !  maling, d.h., 1973: coordinate systems and map projections
        ! equations #20 in appendix i.  
          
        if (chi1 .eq. chi2) then
          chi = 2.0*atan( ( r/tan(chi1) )**(1./proj_cone) 
     &          * tan(chi1*0.5) )
        else
          chi = 2.0*atan( (r*proj_cone/sin(chi1))**(1./proj_cone) 
     &          * tan(chi1*0.5)) 
        endif
        lat = (90.0-chi*deg_per_rad)*proj_hemi

      endif

      if (lon .gt. +180.) lon = lon - 360.
      if (lon .lt. -180.) lon = lon + 360.

      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine llij_lc( lat, lon, i, j)

    ! routine to compute the geographical latitude and longitude values
    ! to the cartesian x/y on a lambert conformal projection.
      
      implicit none

      include 'include_wrf_map_utils'

      ! input args
      real lat      ! latitude (-90->90 deg n)
      real lon      ! longitude (-180->180 e)

      ! output args                 
      real i        ! cartesian x coordinate
      real j        ! cartesian y coordinate

      ! locals 
      real arg
      real deltalon
      real tl1r
      real rm
      real ctl1r
      

      ! begin code
      
      ! compute deltalon between known longitude and standard lon and ensure
      ! it is not in the cut zone
      deltalon = lon - proj_stdlon
      if (deltalon .gt. +180.) deltalon = deltalon - 360.
      if (deltalon .lt. -180.) deltalon = deltalon + 360.
      
      ! convert truelat1 to radian and compute cos for later use
      tl1r = proj_truelat1 * rad_per_deg
      ctl1r = cos(tl1r)     
     
      ! radius to desired point
      rm = proj_rebydx * ctl1r/proj_cone * 
     &    (tan((90.*proj_hemi-lat)*rad_per_deg/2.) / 
     &     tan((90.*proj_hemi-proj_truelat1)*rad_per_deg/2.))**proj_cone

      arg = proj_cone*(deltalon*rad_per_deg)
      i = proj_polei + proj_hemi * rm * sin(arg)
      j = proj_polej - rm * cos(arg)

      ! finally, if we are in the southern hemisphere, flip the i/j
      ! values to a coordinate system where (1,1) is the sw corner
      ! (what we assume) which is different than the original ncep
      ! algorithms which used the ne corner as the origin in the 
      ! southern hemisphere (left-hand vs. right-hand coordinate?)
      if (proj_hemi .eq. -1.) then
        i = 2. - i  
        j = 2. - j
      endif

      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

