***********************************************************************
* FLEXPART SOURCE FILE READ_NCWRFOUT - CONTAINS                       *
*                                                                     *
*    SUBROUTINE READ_NCWRFOUT_GRIDINFO                                *
*    SUBROUTINE READ_NCWRFOUT_1REALFIELD                              *
*    SUBROUTINE READ_NCWRFOUT_1DATETIME                               *
*                                                                     * 
***********************************************************************
*                                                                     * 
*             AUTHOR:      R. Easter & J. Fast, PNNL                  *
*             DATE:        2005-autumn-??                             *
*             LAST UPDATE: same                                       *
*                                                                     * 
***********************************************************************
*                                                                     *
* DESCRIPTION:                                                        *
*                                                                     *
* These routines read the netcdf wrf output files.                    *
*                                                                     *
***********************************************************************


c-----------------------------------------------------------------------
	subroutine read_ncwrfout_gridinfo( ierr, idiagaa, fnamenc,
     &	  n_west_east, n_south_north, n_bottom_top, 
     &	  dx_met, dy_met, 
     &	  m_grid_id, m_parent_grid_id, m_parent_grid_ratio, 
     &	  i_parent_start, j_parent_start, 
     &	  map_proj_id, map_stdlon, map_truelat1, map_truelat2 )
c
c   reads grid definition information from a netcdf wrf output file
c
c   arguments
c	ierr - output - if non-zero, an error occurred 
c		while opening or reading from the file
c	idiagaa - input - if positive, testing diagnostics are printed
c	fnamenc - input - path+filename of the wrf output file
c	n_west_east - output - east_west dimension of the "T-grid"
c	n_south_north - output - south_north dimension of the "T-grid"
c	n_bottom_top - output - bottom_top dimension of the "T-grid"
c	dx_met, dy_met - output - horizontal grid spacing (m)
c	m_grid_id - output - grid id number
c	m_parent_grid_id - output - grid id number of parent grid
c	m_parent_grid_ratio - output - ratio of parent grid dxy to current grid dxy
c	i_parent_start, j_parent_start - output - location of lower left corner
c		of current grid relative to the parent grid.
c	(if there is no parent grid, then the above 4 "...parent..." variables
c		area set to -987.)
c	map_proj_id - WRF map projection id (2=polar stereographic)
c	map_stdlon - map projection standard longitude (deg)
c	map_truelat1, truelat2 - map projection true latitudes (deg)
c
	implicit none

	include 'netcdf.inc'

c   arguments
	integer ierr, idiagaa,
     &  n_west_east, n_south_north, n_bottom_top, 
     &  m_grid_id, m_parent_grid_id, m_parent_grid_ratio, 
     &  i_parent_start, j_parent_start, map_proj_id

	real dx_met, dy_met, map_stdlon, map_truelat1, map_truelat2

	character*(*) fnamenc

c   local variables
	integer maxdim
	parameter (maxdim=20)	! max number of dimensions
	integer ibadaa
	parameter (ibadaa = -987)
	integer xbadaa
	parameter (xbadaa = -987.0)

	integer i, iatt, idimid_unlim, idum, iret, ivtype
	integer l, lenatt, lendim(maxdim)
	integer natts_tot, ncid, ndims_tot, nvars_tot
	integer n_west_east_stag, n_south_north_stag, n_bottom_top_stag 

	real duma
	real, allocatable, dimension(:) :: duma_alloc

	character*80 dimname(maxdim)
	character*80 attname
	character*1000 dumch1000

c   externals
c	integer nf_close
c	integer nf_inq
c	integer nf_inq_dim
c	integer nf_open



c   initialize with "missing values"
	n_west_east = ibadaa
	n_south_north = ibadaa
	n_bottom_top = ibadaa
	dx_met = xbadaa
	dy_met = xbadaa
	m_grid_id = ibadaa
	m_parent_grid_id = ibadaa
	m_parent_grid_ratio = ibadaa
	i_parent_start = ibadaa
	j_parent_start = ibadaa

c
c   open the netcdf file
c
	ncid = 10
	iret = nf_open( fnamenc, NF_NOWRITE, ncid )
	if (iret .ne. nf_noerr) then
	    write(*,9100) 'error doing open', fnamenc
	    ierr = -1
	    return
	end if

9100	format( / '*** read_ncwrfout_gridinfo -- ', a /
     &	'file = ', a )
9110	format( / '*** read_ncwrfout_gridinfo -- ', a, 1x, i8 /
     &	'file = ', a )
9120	format( / '*** read_ncwrfout_gridinfo -- ', a, 2(1x,i8) /
     &	'file = ', a )

90030	format( a, 2i6, 2(2x,a) )
	
c
c get information on dimensions
c
        iret = nf_inq( ncid, 
     &		ndims_tot, nvars_tot, natts_tot, idimid_unlim )
        if (iret .ne. nf_noerr) then
	    write(*,9100) 'error inquiring dimensions', fnamenc
	    ierr = -2
	    return
	end if

	n_west_east_stag = ibadaa
	n_south_north_stag = ibadaa
	n_bottom_top_stag = ibadaa

        do i = 1, min(ndims_tot,maxdim)
            iret = nf_inq_dim( ncid, i, dimname(i), lendim(i) )
            if (iret .ne. nf_noerr) then
		write(*,9110) 'error inquiring dimensions for dim#', i, fnamenc
		ierr = -2
		return
	    end if
	end do

        do i = 1, min(ndims_tot,maxdim)
            if (dimname(i) .eq. 'west_east')   
     &	                       n_west_east = lendim(i)
            if (dimname(i) .eq. 'south_north') 
     &	                       n_south_north = lendim(i)
            if (dimname(i) .eq. 'bottom_top')  
     &	                       n_bottom_top = lendim(i)
            if (dimname(i) .eq. 'west_east_stag')   
     &	                       n_west_east_stag = lendim(i)
            if (dimname(i) .eq. 'south_north_stag') 
     &	                       n_south_north_stag = lendim(i)
            if (dimname(i) .eq. 'bottom_top_stag')  
     &	                       n_bottom_top_stag = lendim(i)
	end do

	if (idiagaa .gt. 0) then
        write(*,9100) 'diagnostics', fnamenc
        do i = 1, min(ndims_tot,maxdim)
            write(*,90030) 'dim #, len, name =', 
     &		i, lendim(i), dimname(i)
	end do
	end if

	if ((n_west_east .le. 0) .or.
     &	    (n_west_east+1 .ne. n_west_east_stag)) then
	    write(*,9120) 'bad n_west_east, n_west_east_stag = ', 
     &		n_west_east, n_west_east_stag, fnamenc
	    ierr = -3
	    return
	end if

	if ((n_south_north .le. 0) .or.
     &	    (n_south_north+1 .ne. n_south_north_stag)) then
	    write(*,9120) 'bad n_south_north, n_south_north_stag = ', 
     &		n_south_north, n_south_north_stag, fnamenc
	    ierr = -3
	    return
	end if

	if ((n_bottom_top .le. 0) .or.
     &	    (n_bottom_top+1 .ne. n_bottom_top_stag)) then
	    write(*,9120) 'bad n_bottom_top, n_bottom_top_stag = ', 
     &		n_bottom_top, n_bottom_top_stag, fnamenc
	    ierr = -3
	    return
	end if

c
c   get information on global attributes
c

c   first just do diagnostics
	if (idiagaa .gt. 0) then
	    write(*,*)
	    write(*,*) 'attribute #, name, type, value'
	end if
	do 3400 iatt = 1, natts_tot
	    iret = nf_inq_attname( ncid, nf_global, iatt, attname)
            if (iret .ne. nf_noerr) goto 3600
	    iret = nf_inq_att( ncid, nf_global, attname, ivtype, lenatt )
            if (iret .ne. nf_noerr) goto 3600
	    if (ivtype .eq. 2) then
		iret = nf_get_att_text( ncid, nf_global, attname, dumch1000 )
        	if (iret .ne. nf_noerr) goto 3600
		i = max(1,min(1000,lenatt))
		if (idiagaa .gt. 0) write(*,91010)
     &			iatt, attname(1:40), ivtype, lenatt, dumch1000(1:i)
	    else if (ivtype .eq. 4) then
		iret = nf_get_att_int( ncid, nf_global, attname, idum )
        	if (iret .ne. nf_noerr) goto 3600
		if (idiagaa .gt. 0) write(*,91020)
     &			iatt, attname(1:40), ivtype, lenatt, idum
	    else if ((ivtype .eq. 5) .and. (lenatt .eq. 1)) then
		iret = nf_get_att_real( ncid, nf_global, attname, duma )
        	if (iret .ne. nf_noerr) goto 3600
		if (idiagaa .gt. 0) write(*,91030)
     &			iatt, attname(1:40), ivtype, lenatt, duma
	    else if ((ivtype .eq. 5) .and. (lenatt .gt. 1)) then
		allocate( duma_alloc(lenatt) )
		iret = nf_get_att_real( ncid, nf_global, attname, duma_alloc )
        	if (iret .ne. nf_noerr) goto 3600
		if (idiagaa .gt. 0) then
		    write(*,91010) iatt, attname(1:40), ivtype, lenatt
		    write(*,91040) (duma_alloc(i), i=1,lenatt)
		end if
		deallocate( duma_alloc )
	    else
		if (idiagaa .gt. 0) write(*,'(i4,1x,a,2(1x,i6))') 
     &			iatt, attname(1:40), ivtype, lenatt
		goto 3400
	    endif

	    if (attname .eq. 'GRID_ID') then
		m_grid_id = idum
	    else if (attname .eq. 'PARENT_ID') then
		m_parent_grid_id = idum
	    else if (attname .eq. 'PARENT_GRID_RATIO') then
		m_parent_grid_ratio = idum
	    else if (attname .eq. 'I_PARENT_START') then
		i_parent_start = idum
	    else if (attname .eq. 'J_PARENT_START') then
		j_parent_start = idum
	    else if (attname .eq. 'DX') then
		dx_met = duma
	    else if (attname .eq. 'DY') then
		dy_met = duma
	    else if (attname .eq. 'MAP_PROJ') then
		map_proj_id = idum
	    else if (attname .eq. 'STAND_LON') then
		map_stdlon = duma
	    else if (attname .eq. 'TRUELAT1') then
		map_truelat1 = duma
	    else if (attname .eq. 'TRUELAT2') then
		map_truelat2 = duma
	    end if

3400	continue
91010	format( i4, 1x, a, 2(1x,i6), 1x, a )
91020	format( i4, 1x, a, 2(1x,i6), 1x, i10 )
91030	format( i4, 1x, a, 2(1x,i6), 1x, 1pe12.4 )
91040	format(( 12x, 5(1pe12.4) ))

	goto 3900

3600	write(*,9110) 'error inquiring attribute', iatt, fnamenc
	stop

3900	continue


c
c   close and return
c
	iret = nf_close( ncid )
	ierr = 0

	return
	end



c-----------------------------------------------------------------------
	subroutine read_ncwrfout_1datetime( ierr, fnamenc,
     &	  itime, jyyyymmdd, jhhmmss )
c
c   a wrf output file may contain data at multiple time.  This routine returns
c	the date & time of the "itime" data group in the file.
c
c   arguments
c	ierr - output - if non-zero, an error occurred 
c		while opening or reading from the file,
c		or itime < 0, or itime > number of times in the file.
c	fnamenc - input - path+filename of the wrf output file
c	itime - input - specifies which data/time to return.  
c		1 for first, 2 for second, ...
c	jyyyymmdd - output - date as 8 decimal digits (yyyymmdd).
c		yyyy=year, mm=month, dd=day of month.
c	jhhmmss - output - time of day as 6 decimal digits (hhmmss).
c		hh=hour, mm=minute, ss=second
c	if (jyyyymmdd=jhhmmss=-1, then ierr is non-zero, and vice-versa)
c

	implicit none

	include 'netcdf.inc'

c   arguments
	integer ierr, itime, jyyyymmdd, jhhmmss
	character*(*) fnamenc

c   local variables
	integer ndims_maxbb
	parameter (ndims_maxbb=4)	! max number of dimensions for a variable

	integer i, id_var, iret, itype_var
	integer iduma, idumb, idumc
	integer id_dim(ndims_maxbb)
	integer istart(ndims_maxbb), icount(ndims_maxbb)
	integer lendim(ndims_maxbb)
	integer natts_tot, ncid, ndims

	character*32 timetext
	character*80 varname, varnamenc

c   externals
c	integer nf_close
c	integer nf_inq
c	integer nf_inq_dim
c	integer nf_open



	jyyyymmdd = -1
	jhhmmss = -1

c
c   open the netcdf file
c
	ncid = 10
	iret = nf_open( fnamenc, NF_NOWRITE, ncid )
	if (iret .ne. nf_noerr) then
	    write(*,9100) 'error doing open', fnamenc
	    ierr = -1
	    goto 8100
	end if

9100	format( / '*** read_ncwrfout_1datetime -- ', a /
     &	'file = ', a )
9110	format( / '*** read_ncwrfout_1datetime -- ', a, 1x, i8 /
     &	'file = ', a )
9120	format( / '*** read_ncwrfout_1datetime -- ', a, 2(1x,i8) /
     &	'file = ', a )
9130	format( / '*** read_ncwrfout_1datetime -- ', a, 3(1x,i8) /
     &	'file = ', a )
9115	format( / '*** read_ncwrfout_1datetime -- ', a / a, 1x, i8 /
     &	'file = ', a )
9125	format( / '*** read_ncwrfout_1datetime -- ', a / a, 2(1x,i8) /
     &	'file = ', a )
9135	format( / '*** read_ncwrfout_1datetime -- ', a / a, 3(1x,i8) /
     &	'file = ', a )

90030	format( a, 2i6, 2(2x,a) )
	
c
c get information on the variable
c
	varname = 'Times'
        iret = nf_inq_varid( ncid, varname, id_var )
        if (iret .ne. nf_noerr) then
c------modified 2011/03/29
	    write(*,9100) 
     &        'error inquiring var id for ' // varname, fnamenc
c------------------------
	    ierr = -1
	    goto 8100
	end if

        iret = nf_inq_var( ncid, id_var, 
     &		varnamenc, itype_var, ndims, id_dim, natts_tot )
        if (iret .ne. nf_noerr) then
c------modified 2011/03/29
	    write(*,9100)
     &        'error inquiring var info for ' // varname, fnamenc
c------------------------
	    ierr = -1
	    goto 8100
	end if

c   check variable type
	if (itype_var .ne. nf_char) then
	    write(*,9110) 'var type wrong for ' // varname, 
     &		itype_var, fnamenc
	    ierr = -1
	    goto 8100
	end if


c   check number of dimensions
	if (ndims .ne. 2) then
	    write(*,9115) 'var ndims is wrong for ' // varname, 
     &		'ndims =', ndims, fnamenc
	    ierr = -1
	    goto 8100
	end if

c   get sizes of dimensions
c   dimension 1 = # of characters in date/time string
c   dimension 2 = # of times in the file
	do i = 1, ndims
	    iret = nf_inq_dimlen( ncid, id_dim(i), lendim(i) )
	    if (iret .ne. nf_noerr) then
c------modified 2011/03/29
		write(*,9115)
     &            'error inquiring var dim len for ' // varname,
     &		  'idim =', i, fnamenc
c-------------------------
		ierr = -1
		goto 8100
	    end if
	end do

	if (itime .lt. 1) then
	    ierr = -11
	    goto 8100
	else if (itime .gt. lendim(2)) then
	    ierr = -12
	    goto 8100
        end if

c   get the data and extract the data & time
	do i = 1, ndims_maxbb
	    istart(i) = 1
	    icount(i) = 1
	end do
	istart(1) = 1
	icount(1) = lendim(1)
	istart(2) = itime
	icount(2) = 1
c------modified 2011/03/29
	iret = nf_get_vara_text(ncid, id_var, istart, icount, timetext)
        if (iret .ne. nf_noerr) then
	    write(*,9100) 'error reading var data for ' // varname, 
     &		fnamenc
	    ierr = -1
	    goto 8100
	end if

	read( timetext, '(i4,1x,i2,1x,i2)', iostat=iret ) 
     &		iduma, idumb, idumc
	if (iret .ne. 0) then
	    write(*,9125) 
     &		'error reading from timetext = "' // timetext // '"',
     &		'itime, lendim(1) =', itime, lendim(1), fnamenc
	    ierr = -1
	    goto 8100
	end if
	jyyyymmdd = iduma*10000 + idumb*100 + idumc

	read( timetext, '(11x,i2,1x,i2,1x,i2)', iostat=iret ) 
     &		iduma, idumb, idumc
	if (iret .ne. 0) then
	    write(*,9125) 
     &		'error reading from timetext = "' // timetext // '"',
     &		'itime, lendim(1) =', itime, lendim(1), fnamenc
	    ierr = -1
	    goto 8100
	end if
	jhhmmss = iduma*10000 + idumb*100 + idumc

c
c   success - close and return
c
	iret = nf_close( ncid )
	ierr = 0
	return

c
c   failure - close and return
c
8100	iret = nf_close( ncid )
	return

	end



c-----------------------------------------------------------------------
	subroutine read_ncwrfout_1realfield( ierr, idiagaa, fnamenc,
     &	  varname, vardata,
     &	  itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim, lendim_exp, lendim_max )
c
c   reads of real (single precision) field at one time from a netcdf wrf output file
c
c   arguments
c	ierr - output - if non-zero, an error occurred 
c		while opening or reading from the file
c		 -1 = error opening file
c		 -2 = requested variable is not in the file
c		 -3 = error while inquiring about the variable
c		 -4 = variable type is other than real (single precision)
c		... = check below, in the code, for explanation of other ierr values.
c	idiagaa - input - if positive, testing diagnostics are printed
c	fnamenc - input - path+filename of the wrf output file
c	varname - input - field name
c	vardata - output - the data for the field
c	itime - input - specifies which time to read.  
c		(1 for first time in the file, 2 for second, ...)
c	ndims - output - number of (netcdf) dimensions of the field.  
c		This includes the time dimension.
c	ndims_exp - input - expected number of dimensions of the field.  
c		An error occurs if ndims .ne. ndims_exp.
c	ndims_max - input - The dimension/size of the lendim_... arrays.
c	lendim - output - The size of each of the "ndims" dimensions.
c	lendim_exp - input - The expected size of each of the "ndims" dimensions.
c		If lendim_exp .gt. 0, then an error occurs if lendim .ne. lendim_exp.
c	lendim_max - input - The maximum size of each dimension.  These are
c		The dimensions of the vardata array.
c
	implicit none

	include 'netcdf.inc'

c   arguments
	integer ierr, idiagaa, itime,
     &	  ndims, ndims_exp, ndims_max,
     &	  lendim(ndims_max), lendim_exp(ndims_max),  
     &	  lendim_max(ndims_max)

	real vardata( lendim_max(1), lendim_max(2), lendim_max(3) )

	character*(*) fnamenc, varname

c   local variables
	integer ndims_maxbb
	parameter (ndims_maxbb=4)	! max number of dimensions for a variable
	integer ibadaa
	parameter (ibadaa = -987)
	integer xbadaa
	parameter (xbadaa = -987.0)
	integer ijktestmax
	parameter (ijktestmax=500)

	integer i, iatt, id_var, iret, itype_var
	integer itot, ijktot, ii
	integer id_dim(ndims_maxbb)
	integer istart(ndims_maxbb), icount(ndims_maxbb)
	integer j, jtot, jj
	integer k, ktot, kk
	integer l, lenatt
	integer lendim_use(ndims_maxbb)
	integer m
	integer natts_tot, ncid, ndiffa, ndiffb

	real duma, dumb
	real testavg(ijktestmax,3,2), 
     &	    testmin(ijktestmax,3,2), testmax(ijktestmax,3,2)

	character*80 varnamenc
	character*80 dimname(ndims_maxbb)

c   externals
c	integer nf_close
c	integer nf_inq
c	integer nf_inq_dim
c	integer nf_open



c
c   open the netcdf file
c
	ncid = 10
	iret = nf_open( fnamenc, NF_NOWRITE, ncid )
	if (iret .ne. nf_noerr) then
	    write(*,9100) 'error doing open', fnamenc
	    ierr = -1
	    return
	end if

9100	format( / '*** read_ncwrfout_1realfield -- ', a /
     &	'file = ', a )
9110	format( / '*** read_ncwrfout_1realfield -- ', a, 1x, i8 /
     &	'file = ', a )
9120	format( / '*** read_ncwrfout_1realfield -- ', a, 2(1x,i8) /
     &	'file = ', a )
9130	format( / '*** read_ncwrfout_1realfield -- ', a, 3(1x,i8) /
     &	'file = ', a )
9115	format( / '*** read_ncwrfout_1realfield -- ', a / a, 1x, i8 /
     &	'file = ', a )
9125	format( / '*** read_ncwrfout_1realfield -- ', a / a, 2(1x,i8) /
     &	'file = ', a )
9135	format( / '*** read_ncwrfout_1realfield -- ', a / a, 3(1x,i8) /
     &	'file = ', a )

90030	format( a, 2i6, 2(2x,a) )
	
c
c get information on the variable
c
        iret = nf_inq_varid( ncid, varname, id_var )
        if (iret .ne. nf_noerr) then
c------modified 2011/03/29
	    write(*,9100)
     &        'error inquiring var id for ' // varname, fnamenc
c-------------------------
	    ierr = -2
	    goto 8100
	end if

        iret = nf_inq_var( ncid, id_var, 
     &		varnamenc, itype_var, ndims, id_dim, natts_tot )
        if (iret .ne. nf_noerr) then
c------modified 2011/03/29
	    write(*,9100)
     &        'error inquiring var info for ' // varname, fnamenc
c------------------------
	    ierr = -3
	    goto 8100
	end if

c   check variable type
	if (itype_var .ne. nf_float) then
	    write(*,9110) 'var type wrong for ' // varname, 
     &		itype_var, fnamenc
	    ierr = -4
	    goto 8100
	end if


c   check number of dimensions
	if (ndims_exp .le. 0) then
	    write(*,9115) 
     &		'bad ndims_exp for ' // varname, 
     &		'ndims_exp =', ndims_exp, fnamenc
	    ierr = -11
	    goto 8100
	end if
	if (ndims .ne. ndims_exp) then
	    write(*,9125) 'var ndims mismatch for ' // varname, 
     &		'ndims_exp, ndims =', ndims_exp, ndims, fnamenc
	    ierr = -12
	    goto 8100
	end if
	if (ndims .gt. ndims_max) then
	    write(*,9125) 'var ndims > ndims_max for ' // varname, 
     &		'ndims, ndims_max =', ndims, ndims_max, fnamenc
	    ierr = -13
	    goto 8100
	end if
	if (ndims .gt. ndims_maxbb) then
	    write(*,9125) 'var ndims > ndims_maxbb for ' // varname, 
     &		'ndims, ndims_maxbb =', ndims, ndims_maxbb, fnamenc
	    ierr = -14
	    goto 8100
	end if

c   check size of each dimension
	do i = 1, ndims_exp
            iret = nf_inq_dimlen( ncid, id_dim(i), lendim(i) )
            if (iret .ne. nf_noerr) then
c------modified 2011/03/29
		write(*,9110)
     &            'error inquiring var dim len for ' // varname,
c-------------------------
     &			i, fnamenc
		ierr = -15
		goto 8100
	    end if
	    if ((i .lt. ndims_exp) .and. (lendim_exp(i) .gt. 0) .and. 
     &	        (lendim(i) .ne. lendim_exp(i))) then
		write(*,9130) 'var lendim mismatch for ' // varname, 
     &		    i, lendim_exp(i), lendim(i), fnamenc
		ierr = -16
		goto 8100
	    end if
	    if ((i .lt. ndims_exp) .and. 
     &		(lendim(i) .gt. lendim_max(i))) then
		write(*,9130) 'var lendim too big for ' // varname, 
     &		    i, lendim_max(i), lendim(i), fnamenc
		ierr = -17
		goto 8100
	    end if
	    if ((i .eq. ndims_exp) .and. (lendim(i) .lt. itime)) then
		write(*,9130) 'var itime < ntimes for ' // varname, 
     &		    i, itime, lendim(i), fnamenc
		ierr = -18
		goto 8100
	    end if
	end do

c   do diagnostics on the dimensions
	if (idiagaa .gt. 0) then
	  write(*,'(/a)') 
     &		'read_ncwrfout_1realfield - dim info for var = ' // 
     &		varname(1:20)
	  do i = 1, ndims
            iret = nf_inq_dim( ncid, id_dim(i), dimname(i), lendim(i) )
            if (iret .ne. nf_noerr) then
c------modified 2011/03/29
		write(*,9115)
     &            'error inquiring var dim info for ' // varname,
     &		  'idim =', i, fnamenc
c------------------------
		ierr = -19
		goto 8100
	    end if
	    write(*,'(a,3i5,2x,a)') '     i,id,len,name =', 
     &		i, id_dim(i), lendim(i), dimname(i)(1:32)
	  end do
	end if

c
c   get the data
c
	do i = 1, ndims_maxbb
	    istart(i) = 1
	    icount(i) = 1
	end do
	do i = 1, ndims_exp - 1
	    istart(i) = 1
	    icount(i) = lendim(i)
	end do
c   in wrfout files, the last dimension should always be time
	istart(ndims_exp) = itime
	icount(ndims_exp) = 1

	iret = nf_get_vara_real( ncid, id_var, istart, icount, vardata )
        if (iret .ne. nf_noerr) then
	    write(*,9120) 'error reading var data for ' // varname, 
     &		fnamenc
	    ierr = -21
	    goto 8100
	end if

c
c   reorder the data
c
	lendim_use(1) = lendim(1)
	lendim_use(2) = 1
	if (ndims_exp .ge. 3) lendim_use(2) = lendim(2)
	lendim_use(3) = 1
	if (ndims_exp .ge. 4) lendim_use(3) = lendim(3)

	call reorder_ncwrfout_1realfield( ierr, idiagaa, 
     &	  varname, vardata, vardata,
     &	  lendim_use, lendim_max )

	if (ierr .ne. 0) then
	    write(*,9120) 'error re-ordering var data for ' // varname, 
     &		fnamenc
	    ierr = -22
	    goto 8100
	end if

c
c   success - close and return
c
	iret = nf_close( ncid )
	ierr = 0
	return

c
c   error - close and return
c
8100	iret = nf_close( ncid )
	return

	end



c-----------------------------------------------------------------------
	subroutine reorder_ncwrfout_1realfield( ierr, idiagaa, 
     &	  varname, vardata_in, vardata_out,
     &	  lendim_use, lendim_max )
c
c   reorders a real (single precision) field 
c	the nf_get_vara_real loads the data for a field into 
c	    a contiguous block of memory starting at vardata(1,1,1)
c	it does not know if perhaps lendim() < lendim_max()
c	this routine corrects for that, so that the data are
c	    loaded into non-contiguous blocks when lendim() < lendim_max()
c
c   arguments
c	ierr - output - if non-zero, an error occurred while re-ordering
c	idiagaa - input - if positive, testing diagnostics are printed
c	varname - input - field name
c
c	vardata_in  - input  - the data for the field
c	vardata_out - output - the data for the field
c		In the calling program, vardata_in & vardata_out are usually
c		the same array.  This routine "pretends" that they are 
c		different, and specifies their dimensions differently,
c		to facilitate the reordering.
c
c	lendim_use - input - The actual size of the spatial dimensions of the field.
c	lendim_max - input - The actual spatial dimensions of the vardata array.
c		Most wrf fields are spatially either 1d (z), 2d (xy), or 3d (xyz).
c		For a 1d spatial field (e.g., z only), set
c		    lendim_use(1) = nz,  lendim_max(1) = nz_max
c		    lendim_use(2) = 1,   lendim_max(2) = 1
c		    lendim_use(3) = 1,   lendim_max(3) = 1
c		For a 2d spatial field (e.g., xy only), set
c		    lendim_use(1) = nx,  lendim_max(1) = nx_max
c		    lendim_use(2) = ny   lendim_max(2) = ny_max
c		    lendim_use(3) = 1,   lendim_max(3) = 1
c		For a 3d spatial field (xyz), set
c		    lendim_use(1) = nx,  lendim_max(1) = nx_max
c		    lendim_use(2) = ny   lendim_max(2) = ny_max
c		    lendim_use(3) = nz   lendim_max(3) = nz_max
c

	implicit none

	include 'netcdf.inc'

c   arguments
	integer ierr, idiagaa, 
     &	  lendim_use(3), lendim_max(3)

	real vardata_in(  lendim_use(1), lendim_use(2), lendim_use(3) )
	real vardata_out( lendim_max(1), lendim_max(2), lendim_max(3) )

	character*(*) varname

c   local variables
	integer ijktestmax
	parameter (ijktestmax=500)
	integer check_reordering
	parameter (check_reordering = 1)

	integer i, j, k, m, n
	integer itestend, jtestend, ktestend
	integer ijk, ijktestend
	integer ndiffa, ndiffb

	real duma, dumb
	real testavg(ijktestmax,3,2), 
     &	    testmin(ijktestmax,3,2), testmax(ijktestmax,3,2)

c
c   the testavg/min/max are avg, min, and max values for
c	a i (or j or k) fixed and j,k (or i,k or i,j) varying
c   they are computed before and after the data has been reordered
c	then compared at the end
c   an error occurs if they do not match
c
	if (check_reordering .gt. 0) then

	do n = 1, 2
	do m = 1, 3
	do i = 1, ijktestmax
	    testavg(i,m,n) = 0.0
	    testmin(i,m,n) = +1.0e37
	    testmax(i,m,n) = -1.0e37
	end do
	end do
	end do

	ktestend = min( ijktestmax, lendim_use(3) )
	jtestend = min( ijktestmax, lendim_use(2) )
	itestend = min( ijktestmax, lendim_use(1) )

c   pass 1 -- compute the test---(:,:,1) from vardata_in
	do k = 1, ktestend
	do j = 1, jtestend
	do i = 1, itestend
	    duma = vardata_in(i,j,k)
	    testavg(i,1,1) =      testavg(i,1,1) + duma
	    testmin(i,1,1) = min( testmin(i,1,1),  duma )
	    testmax(i,1,1) = max( testmax(i,1,1),  duma )
	    testavg(j,2,1) =      testavg(j,2,1) + duma
	    testmin(j,2,1) = min( testmin(j,2,1),  duma )
	    testmax(j,2,1) = max( testmax(j,2,1),  duma )
	    testavg(k,3,1) =      testavg(k,3,1) + duma
	    testmin(k,3,1) = min( testmin(k,3,1),  duma )
	    testmax(k,3,1) = max( testmax(k,3,1),  duma )
	end do
	end do
	end do

	end if    ! if (check_reordering .gt. 0) then

c   pass 2 -- shift the data values
	do k = lendim_use(3), 1, -1
	do j = lendim_use(2), 1, -1
	do i = lendim_use(1), 1, -1
	    vardata_out(i,j,k) = vardata_in(i,j,k)
	end do
	end do
	end do

c   pass 3 -- compute the test---(:,:,2) from vardata_out
	if (check_reordering .gt. 0) then

	do k = 1, ktestend
	do j = 1, jtestend
	do i = 1, itestend
	    duma = vardata_out(i,j,k)
	    testavg(i,1,2) =      testavg(i,1,2) + duma
	    testmin(i,1,2) = min( testmin(i,1,2),  duma )
	    testmax(i,1,2) = max( testmax(i,1,2),  duma )
	    testavg(j,2,2) =      testavg(j,2,2) + duma
	    testmin(j,2,2) = min( testmin(j,2,2),  duma )
	    testmax(j,2,2) = max( testmax(j,2,2),  duma )
	    testavg(k,3,2) =      testavg(k,3,2) + duma
	    testmin(k,3,2) = min( testmin(k,3,2),  duma )
	    testmax(k,3,2) = max( testmax(k,3,2),  duma )
	end do
	end do
	end do

c   now compare the test---(:,:,1) & test---(:,:,2)
	ndiffb = 0
	do m = 1, 3
	    if (m .eq. 1) then
		ijktestend = itestend
		duma = 1.0/(jtestend*ktestend)
		if (idiagaa .gt. 0) write(*,'(a,a)') varname(1:20),
     &			'i, testavg(i,1), testmin(i,1), testmax(i,1)'
	    else if (m .eq. 2) then
		ijktestend = jtestend
		duma = 1.0/(itestend*ktestend)
		if (idiagaa .gt. 0) write(*,'(a,a)') varname(1:20),
     &			'j, testavg(j,2), testmin(j,2), testmax(j,2)'
	    else
		ijktestend = ktestend
		duma = 1.0/(itestend*jtestend)
		if (idiagaa .gt. 0) write(*,'(a,a)') varname(1:20),
     &			'k, testavg(k,3), testmin(k,3), testmax(k,3)'
	    end if

	    ndiffa = 0
	    do ijk = 1, ijktestend
		i = ijk
c------modified 2011/03/29
		dumb = 
     &           max( abs(testavg(i,m,1)), abs(testavg(i,m,2)) )*2.0e-7
		if (abs(testavg(i,m,1)-testavg(i,m,2)) .gt. dumb)
     &              ndiffa = ndiffa + 1
		dumb =
     &           max( abs(testmin(i,m,1)), abs(testmin(i,m,2)) )*2.0e-7
		if (abs(testmin(i,m,1)-testmin(i,m,2)) .gt. dumb)
     &              ndiffa = ndiffa + 1
		dumb =
     &           max( abs(testmax(i,m,1)), abs(testmax(i,m,2)) )*2.0e-7
		if (abs(testmax(i,m,1)-testmax(i,m,2)) .gt. dumb)
     &              ndiffa = ndiffa + 1
	    end do

	    if (ndiffa .le. 0) then
		if (idiagaa .gt. 0) write(*,*) '     *** no differences'
	    else
	      do ijk = 1, ijktestend
		i = ijk
		if (idiagaa .gt. 0) write(*,'(i3,1p,3(2x,2e11.3))') i,
     &		    testavg(i,m,1)*duma, 
     &		    (testavg(i,m,1)-testavg(i,m,2))*duma, 
     &		    testmin(i,m,1), (testmin(i,m,1)-testmin(i,m,2)), 
     &		    testmax(i,m,1), (testmax(i,m,1)-testmax(i,m,2))
	      end do
	    end if

	    ndiffb = ndiffb + ndiffa
	end do

	if (ndiffb .gt. 0) then
	    ierr = -12
	    goto 8100
	end if

	end if    ! if (check_reordering .gt. 0) then

c
c   success
c
	ierr = 0
	return

c
c   error
c
8100	return

	end



