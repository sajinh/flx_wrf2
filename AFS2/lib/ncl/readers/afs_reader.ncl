function only_season_or_month()
begin
  return("can use only one of options, season or months")
end

function mme_in_path(db_root,store,cabinet,model,mon,lead,yr,opt)
begin
  month=integer_to_month(mon)
  if opt@read_mme_out .eq. True
    if isatt(opt,"season") .and. str_upper(opt@season) .ne. "NONE"
      month=month+"/"+str_upper(opt@season)
    end if
  end if
  indir=lead+"-MON"+"/"+store+"/"+cabinet+"/"+model+"/"+month+"/"+yr
  return(db_root+"/"+str_upper(indir))
end

function afs_dataout_reader(var,opt)
begin
  return(var)
end

function afs_data_monthly_value(var,months)
begin
  new_var=ensure_last_dim(var,"time")
  opt=True
  opt@preserve_coordinates=True
  mon_val=cut_from_last_dim(new_var,months,opt)
  return(reorder_dims_as_in_template(var,mon_val))
end

function afs_data_reader(data_dir,data_period,var_name,opt)
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

  fin=addfile(infil,"r")
  lat_rng=afs_dim_range(fin,"lat",opt)
  lon_rng=afs_dim_range(fin,"lon",opt)
  if isatt(opt,"nlev")
    nlev=opt@nlev
    var=fin->$var_name$(:,0:nlev-1,{lat_rng(0):lat_rng(1)},{lon_rng(0):lon_rng(1)})
  else
    var=fin->$var_name$(:,:,{lat_rng(0):lat_rng(1)},{lon_rng(0):lon_rng(1)})
  end if

  if opt@read_mme_out .eq. True
    return(afs_dataout_reader(var,opt))
  end if

  ;------------code only for mme_in ; start -------------
  if isatt(opt,"season")
    season=str_upper(opt@season)
    if season.ne."NONE"
      if isatt(opt,"months") .and. any(opt@months .ne. "NONE")
        terminate_program(only_season_or_month)
      end if
      time_opts = data_period
      t1 = ind(_seasons .eq. season) + 1 - time_opts(0)
      if t1 .lt. 0
        t1 = t1 + 12
      end if
      vseas = var(lat|:,lon|:,level|:,time|t1:t1)
      vseas= dim_avg(var(lat|:,lon|:,level|:,time|t1:t1+2))
      delete(var)
      var = vseas(time|:,level|:,lat|:,lon|:)
    end if
  end if

  if (isatt(opt,"months") .and. (all(opt@months .ne. "NONE")))
    if isatt(opt,"season") .and. opt@season .ne. "NONE"
      terminate_program(only_season_or_month)
    end if
    time_opts = data_period
    req_month = month_to_int(opt@months)
    req_month = req_month - time_opts(0)
    req_month = where(req_month .lt. 0, req_month + 12, req_month)
    if any(req_month .lt. 0)
      msg1="Please provide valid month to retrieve"
      msg2="  Start Month     == "+time_opts(0)
      msg3="  Requested Month == "+opt@months
      terminate_program(msg1+cr()+msg2+cr()+msg3)
    end if
    vmon = afs_data_monthly_value(var,req_month)
    delete(var)
    var=vmon
  end if


  if isatt(opt,"ensemble_average")
    if opt@ensemble_average
      return(dim_avg_Wrap(var(time|:,lat|:,lon|:,level|:)))
    end if
  end if

  ;------------code only for mme_in ; end -------------

  return(var)
end
