load "$AFS2/lib/ncl/readers/afs_reader.ncl"

function mme_in_reader(db_root,data_type,model,var_name,data_period,opt)
begin
  ; mme_in_reader opens multiple year data
  ; and return a data object
  ;(options = DataType, DataPeriod, EnsembleAverage)
  ;(options = mme_out to read_mme out data)
  
  opt@read_mme_out = False
  time_opts        = data_period
  start_mon        = time_opts(0)
  end_mon          = time_opts(1)
  years            = time_opts(2:)
  nyr              = dimsizes(years)
  first_year       = years(0)
  lead_time=end_mon-start_mon+1
  data_dir=mme_in_path(db_root,"MME_IN",data_type,model,start_mon,lead_time,\
                           first_year,opt)

  var0=afs_data_reader(data_dir,data_period,var_name,opt)
  dims=dimsizes(var0)
  nmon=dims(0)

  named_dims=getvardims(var0)
  if named_dims(1).eq."level"
    rank=4
    nlev=dims(1)
    nlat=dims(2)
    nlon=dims(3)
    if isatt(opt,"constrain_nlev")
      if opt@constrain_nlev
        opt@nlev=nlev
      end if
    end if
    var=new( (/nmon,nlev,nlat,nlon,nyr/), typeof(var0) )
    var!4="year"
    var(:,:,:,:,0)=var0
  else
    rank=3
    nlat=dims(1)
    nlon=dims(2)
    var=new( (/nmon,nlat,nlon,nyr/), typeof(var0) )
    var(:,:,:,0)=var0
    var!3="year"
  end if
  delete(var0)
 
  do iyr=1,nyr-1
    data_dir=mme_in_path(db_root,"MME_IN",data_type,model,start_mon,\
                           lead_time,years(iyr),opt)
    var0=afs_data_reader(data_dir,data_period,var_name,opt)
    if rank .eq. 4
      var(:,:,:,:,iyr)=var0
    else
      var(:,:,:,iyr)=var0
    end if
  end do
  var&year=years
  var@variable=var_name

  delete(var&time)
  if isatt(opt,"season") .and. opt@season .ne. "NONE"
    var!0="season"
    var&season=str_upper(opt@season)
    if rank .eq. 4
      return(var(year|:,season|:,level|:,lat|:,lon|:))
    else
      return(var(year|:,season|:,lat|:,lon|:))
    end if
  else
    var!0="month"
    if isatt(opt,"months") .and. any(opt@months .ne. "NONE")
      var&month=str_upper(opt@months)
    else
      var&month=month_to_string(ispan(start_mon,end_mon,1))
    end if
    if rank .eq. 4
      return(var(year|:,month|:,level|:,lat|:,lon|:))
    else
      return(var(year|:,month|:,lat|:,lon|:))
    end if
  end if
end

function multi_year_predictions(db_root,data_type,models,\
                     var_name,start_mon,lead_time,start_years,opt)
begin
  data_period=new( dimsizes(start_years)+2, integer)
  data_period(0) = start_mon
  data_period(1) = start_mon+lead_time-1
  data_period(2:) = start_years
  opt@ensemble_average=True
  var0=mme_in_reader(db_root,data_type,models(0),var_name,data_period,opt)
  dims=dimsizes(var0)
  nyr=dims(0)
  nmon=dims(1)
  nlat=dims(2)
  nlon=dims(3)
  nmodel=dimsizes(models)
  var=new( (/nyr,nmon,nlat,nlon,nmodel/), typeof(var0) )
  var!4="models"
  var&models=models
  var(:,:,:,:,0) = var0
  do imodel=1,nmodel-1
    var(:,:,:,:,imodel) = mme_in_reader(db_root,data_type,\
                    models(imodel),var_name,data_period,opt)
  end do
  return(var)
end

function single_year_predictions(db_root,data_type,models,\
var_name,start_mon,lead_time,start_year,opt)
begin
  return(multi_year_predictions(db_root,data_type,models,\
  var_name,start_mon,lead_time,start_year,opt))
end

function hindcast_clim(db_root,models,start_mon,lead_time, \
                      start_years, var_name, opt)
begin
  data_type="Hindcast"
  var= multi_year_predictions(db_root,data_type,models,\
                     var_name,start_mon,lead_time,start_years,opt)
  return(average_along_given_dim(var,"year"))
end

function single_year_predictions_as_anom(db_root,data_type,models,\
                     var_name,start_mon,lead_time,start_year,opt)
begin
  var=single_year_predictions(db_root,data_type,models,\
                     var_name,start_mon,lead_time,start_year,opt)

  ; calculate climatology based on the hindcast period
  hindcast_years = find_correct_hindcast_years(db_root,data_type,var_name,models,\
                   start_mon,lead_time, start_year, opt)
  
  var_clim= hindcast_clim(db_root,models,start_mon,lead_time, \
                              hindcast_years,var_name,opt) 
  vclim = conform_dims( dimsizes(var), var_clim, (/1,2,3,4/) )
  vanom = var
  vanom = var - vclim
  vanom@hindcast_years=hindcast_years
  return(vanom)
end

function multi_year_hindcasts_as_anom(db_root,models,\
                     var_name,start_mon,lead_time,start_years,opt)
begin
  data_type="Hindcast"
  var=multi_year_predictions(db_root,data_type,models,\
                     var_name,start_mon,lead_time,start_years,opt)
  var_clim = average_along_given_dim(var,"year")
  vclim = conform_dims( dimsizes(var), var_clim, (/1,2,3,4/) )
  vanom = var
  vanom = var - vclim
  return(vanom)
end

;function get_model_hindcast_mme_clim(lead_time,models,start_mon, \
;                      start_years, var_name, opt)
;begin
;  data_type="Hindcast"
;  var=get_model_data(db_root,data_type,lead_time,models,\
;                     var_name,start_mon,start_years,opt)
;  return (wgt_areaave_Wrap(var,1,1,0))
;end

