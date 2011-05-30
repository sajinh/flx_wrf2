undef("only_season_or_month")
function only_season_or_month()
begin
  return("can use only one of options, season or months")
end

function afs_obsdata_reader(data_dir,var_name,mon,years,opt)
begin
  infil = data_dir+"/"+var_name+".nc"

  if .not. isfilepresent(infil)
    if isatt(opt,"terminate_if_no_data") \
             .and. (opt@terminate_if_no_data .eq. False)
      return(AFS_FillValue)
    else
      afs_abort(infil+" not present")
    end if
  end if

  fin = addfile(infil,"r")
  lat_rng=afs_dim_range(fin,"lat",opt)
  lon_rng=afs_dim_range(fin,"lon",opt)


  time=fin&time
  nyr=dimsizes(years)
  tunits=time@units
  t0 = time(0)
  tdat = get_yr_mo_dy(t0)
  smon=ut_inv_calendar(years(0),mon,tdat@day,0,0,0,tunits,0)
  emon=ut_inv_calendar(years(nyr-1),mon,tdat@day,0,0,0,tunits,0)
  istart=ind(time.eq.smon)
  iend=ind(time.eq.emon)
  return(fin->$var_name$(istart:iend:12,\
                    {lat_rng(0):lat_rng(1)},{lon_rng(0):lon_rng(1)}))
end

function get_obs_data(data_dir,var_name,start_mon,nmon,years,opt)
begin
  ; default case
  months = ispan(start_mon,start_mon+nmon-1,1)
  nmon2 = nmon
 
  ; if options "months" or "season" are passed change the rule
  if isatt(opt,"months") .and. all(opt@months .ne. "NONE")
    if isatt(opt,"season") .and. opt@season .ne. "NONE"
      terminate_program(only_season_or_month)
    end if
    delete(months)
    months=month_to_int(opt@months)
  end if

  if isatt(opt,"season") .and. opt@season .ne. "NONE"
    if isatt(opt,"months") .and. any(opt@months .ne. "NONE")
      terminate_program(only_season_or_month)
    end if
    delete(months)
    months=months_for_season(opt@season)  + 1
  end if

  nmon2=dimsizes(months)
  smonths=month_to_string(months)
  needed_month=get_obs_mon_year(start_mon,months(0),years)
  needed_years=needed_month@years
  var0=afs_obsdata_reader(data_dir,var_name,needed_month,needed_years,opt)
  dims=dimsizes(var0)
  nlat=dims(1)
  nlon=dims(2)
  var=new( (/dimsizes(years),nlat,nlon,nmon2/), typeof(var0) )
  var(:,:,:,0)=var0
  
  do imon = 1,nmon2-1
    needed_month=get_obs_mon_year(start_mon,months(imon),years)
    needed_years=needed_month@years
    var(:,:,:,imon)=afs_obsdata_reader(data_dir,var_name,needed_month,needed_years,opt)
  end do
  var!0="year"
  var&year=years
  var@variable=var_name
  var!1="lat"
  var!2="lon"
  var!3="month"
  if isatt(opt,"season") .and. opt@season .ne. "NONE"
    vopt=True
    vopt@preserve_coordinates=True
    vseas = cut_from_last_dim(ensure_last_dim(var,"month"),0,vopt)
    vseas = reorder_dims_as_in_template(vseas,average_along_given_dim(var,"month"))
    new_var=reorder_dims_as_in_template(var,vseas)
    new_var!3="season"
    new_var&season=opt@season
    return(new_var(year|:,season|:,lat|:,lon|:))
  end if
  var&month=smonths
  return(var(year|:,month|:,lat|:,lon|:))
end

function get_obs_anom(data_dir,var_name,start_mon,nmon,years,opt)
begin
  var=get_obs_data(data_dir,var_name,start_mon,nmon,years,opt)
  vanom=anomaly_along_dim(var,"year")
  if isatt(opt,"season") .and. opt@season .ne. "NONE"
    return(reorder_dims_like_this((/"year","season","lat","lon"/),vanom))
  else
    return(reorder_dims_like_this((/"year","month","lat","lon"/),vanom))
  end if
end
