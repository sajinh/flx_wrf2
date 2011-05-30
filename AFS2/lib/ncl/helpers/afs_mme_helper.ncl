function afs_dim_range(fin,dim,opt)
; provide desired range of desired dimension (say lat)
; separated by ":" example opt@lat="30:120"
begin
  if isatt(opt,dim)
    dim_var=str_split(opt@$dim$,":")
    dim1=dim_var(0)
    dim2=dim_var(1)
  else
    dim_var=fin->$dim$
    dim1=dim_var(0)
    dim2=dim_var(dimsizes(dim_var)-1)
  end if
  return((/to_float(dim1),to_float(dim2)/))
end

function find_hindcast_years(db_root,models,var_name,start_mon,lead_time,opt)

begin
  nmodels = dimsizes(models)
  minyears = new(nmodels,integer)
  maxyears = new(nmodels,integer)
  mon=integer_to_month(start_mon)
  TopDir="/"+lead_time+"-MON"
  RootDir=db_root

  Store = "MME_IN"
  if isatt(models,"Store")
    Store = models@Store
  end if

  do imod = 0, nmodels-1
    OutDir=TopDir+"/"+Store+"/Hindcast/"+models(imod)+"/"+mon+"/"
    OutDir=RootDir+"/"+str_upper(OutDir)
    ls_cmd="ls "+OutDir+"/*/"+var_name+".nc"
    tr_cmd=" | tr -s / "
    awk_cmd=" | awk  -F"+models(imod)+"/"+mon+"/ '{print $2}'"
    awk_spl=" | awk  -F/ '{print $1}'"
    years=stringtoint(systemfunc(ls_cmd+tr_cmd+awk_cmd+awk_spl))
    minyears(imod)=min(years)
    maxyears(imod)=max(years)
    delete(years)
  end do

  start_year=max(minyears)
  end_year=min(maxyears)
  start_end_years=(/start_year,end_year/)
  return(start_end_years)
end

function common_training_period(db_root,models[*]:string,start_mon,lead_time,opt)
local nmodels, minyears, maxyears, imod,  years, \
       start_year, end_year, start_end_years
begin

  nmodels = dimsizes(models)
  minyears = new(nmodels,integer)
  maxyears = new(nmodels,integer)
  mon=integer_to_month(start_mon)
  TopDir="/"+lead_time+"-MON"
  RootDir=db_root

  Store = "MME_IN"
  if isatt(models,"Store")
    Store = models@Store
  end if

  do imod = 0, nmodels-1
    OutDir=TopDir+"/"+Store+"/Hindcast/"+models(imod)+"/"+mon+"/"
    OutDir=RootDir+"/"+str_upper(OutDir)
    years=stringtoint(systemfunc("ls "+OutDir))
    minyears(imod)=min(years)
    maxyears(imod)=max(years)
    delete(years)
  end do

  start_year=max(minyears)
  end_year=min(maxyears)
  start_end_years=(/start_year,end_year/)
  return(start_end_years)
end

function find_correct_hindcast_years(db_root,data_type, var_name, models,\
               i_start_mon,lead_time,start_year, opt)
begin
  hindcast_years = find_hindcast_years(db_root,models,\
                     var_name,i_start_mon,lead_time,opt)
  start_years = ispan(hindcast_years(0),hindcast_years(1),1)
  if str_upper(data_type) .eq. "FORECAST"
    return(start_years)
  end if
  wanted_yr_indices=ind(start_years.ne.start_year)
  correct_hindcast_years=start_years(wanted_yr_indices)
  return(correct_hindcast_years)
end
