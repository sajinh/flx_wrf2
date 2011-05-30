; subroutines to add data to the CMS Database
; use custom functions to calculate climatology
; these will calculate climatology and write
; appropriate information and also calculate
; std dev or  other appropriate measures
load "$AFS2/lib/ncl/writers/writem.ncl"

; + - + 
;procedure remove_unnecessary_attribute_of(var,att)
;begin
; if isatt(var,att)
;   delete(var@$att$)
; end if
;end

procedure check_rank_of(var)
local var
begin
  ; Check if we can process the variable using this module
  nrank=dimsizes(dimsizes(var))
  if nrank.lt.3 .or. nrank.gt.4
    print("This module cannot process your order : please check manual")
    printVarSummary(var)
    print("eXit")
    exit
  end if
  var@nrank=nrank
end

function create_outdir(darray)
local outdir
begin
  darray@join="/"
return(array_to_string(darray))
end

procedure remove_xtra_slashes_of(aString)
local aString
begin
  tmp=aString
  aString=systemfunc("echo "+tmp+" | tr -s '//' '/'")
end

function name_of_outfil(outdir,vname)
begin
  outfil=outdir+"/"+vname+".nc"
return(outfil)
end
; ++++
function create_file_if_not_exists(outdir,vname)
local outdir,vname,outfil
begin
  outfil = name_of_outfil(outdir,vname)
  print("Checking existence of :"+outfil)
  if isfilepresent(outfil)
    print("Output file exists : eXit")
    TerminateProgram
  end if
return(outfil)
end
; ++++
procedure check_time(var)
begin
  dims=getvardims(var)
  ndims=dimsizes(dims)
  do i = 0,ndims-1
    if ismissing(dims(i))
      continue
    end if
    if dims(i).eq."time"
      var@_tindx=i
      return(i)
    end if
  end do
terminate_program("sorry, your data has no time dimension")
end
; ++++
procedure check_lev(var)
begin
  dims=getvardims(var)
  ndims=dimsizes(dims)
  do i = 0,ndims-1
    if ismissing(dims(i))
      continue
    end if
    if dims(i).eq."level"
      var@_lindx=i
      return(i)
    end if
  end do
terminate_program("sorry, your data has no level dimension")
end
; ++++
procedure check_time_object_of(var)
local time
begin
  check_time(var)
  time=var&time
  if .not. isatt(time,"units")
    terminate_program("time object does not have time units")
  end if
  if .not. isatt(time,"calendar")
    terminate_program("please specify a calendar type for time, e.g: standard")
  end if
end
; ++++
procedure make_standard_time_axis_for(var)
local time
begin
  time=var&time
  ntim=dimsizes(time)
  D1=get_yr_mo_dy(time(0))
  D2=get_yr_mo_dy(time(ntim-1))
;  tim2=to_std_time(time)
; add code to write time depending on data frequency
; all monthly time to written on mid-month
  tunits="days since 1900-01-01"

  y1=D1@year
  y2=D2@year
  m1=D1@month
  m2=D2@month
  d1=D1@day
  d2=D2@day

  if var@DataFreq.eq."Monthly"
    tim2=make_monthly_time(y1,m1,15,y2,m2,15,tunits,ntim)
  end if

  if var@DataFreq.eq."Daily"
    tim2=make_daily_time(y1,m1,d1,y2,m2,d2,tunits,ntim)
  end if


  delete(var&time)
  tindx=var@_tindx
  var!tindx="time"
  var&time=tim2
  delete(var@_tindx)
  print(var@time1)
  print(D1)
  if isatt(var,"time1")
    delete(var@time1)
    delete(var@time2)
  end if
  var@time1=D1
  var@time2=D2
end
; + +
procedure update_or_new(var,vname)
begin
  vname@_update=False
  if isatt(var,"update")
    vname@_update=var@update
  end if
end
; + +
procedure var_update(outdir,var_name,vname,var1)
begin
; 1. open the existing file (outdir/vname.nc)
; 2. find last time of var1 in existing file
; 3. if current time is consistent with last time, append new data
  ;printVarSummary(var1)
  remove_unnecessary_attribute_of(vname,"_update")
  outfil=name_of_outfil(outdir,var_name)
  fin1  = addfile(outfil,"r")
  tim   = fin1->time
  ntim0 = dimsizes(tim)
  tim1  = tim(ntim0-1)
  if var1@DataFreq .eq. "Monthly"
    next_time = next_month_of(tim1)
  end if
  if var1@DataFreq .eq. "Daily"
    next_time = next_day_of(tim1)
  end if
  delete(tim)
  delete(tim1)
  tim   = var1&time
  ntim1 = dimsizes(tim)
  tim1  = tim(0)
  delete(tim)
  if (tim1.eq.next_time) then
    print(tim1+"::"+next_time+" : End times match")
    print("Cat the updated data with the original data")
    var0 = fin1->$vname$
    nlat = dimsizes(var0&lat)
    nlon = dimsizes(var0&lon)
    var2 = new( (/ntim0+ntim1,nlat,nlon/), typeof(var0) )
    var2!0 = "time"
    newtime           = new( ntim0+ntim1, typeof(var0&time) )
    newtime(:ntim0-1)  = var0&time
    newtime(ntim0:)    = var1&time
    var2&time         = newtime
    var2(:ntim0-1,:,:) = var0
    var2(ntim0:,:,:)   = var1
    print("Deleting old file")
    system("rm "+outfil)
    print("Adding new file with updated data")
    fout      = addfile(outfil,"c")
    write_TLL(fout,vname,var2)
  else
    print("What is wrong ?")
    terminate_program(tim1+"::"+next_time+" : End times do not match")
  end if
  delete(tim1)
  delete(next_time)
end
; +++
procedure write_4D_data(var,outdir,tim1,tim2)
begin
  vname     = var@vnam
  levs      = var&level
  nlev     = dimsizes(levs)
  remove_unnecessary_attribute_of(var,"_lindx")
  do ilev   = 0, nlev-1
    var1=var({time|tim1:tim2},level|ilev,lat|:,lon|:)
    var_name  = vname+"_"+levs(ilev)
    update_or_new(var1,vname)   ; check if output file is to be updated
    remove_unnecessary_attribute_of(var1,"update")
    remove_unnecessary_attribute_of(var1,"vnam")
    if vname@_update
      var_update(outdir,var_name,vname,var1)
    else                       ; start a new file
      remove_unnecessary_attribute_of(vname,"_update")
      outfil    = create_file_if_not_exists(outdir,var_name)
      fout      = addfile(outfil,"c")
      write_TLL(fout,vname,var1)
    end if
    delete(var1)
  end do
end

procedure write_3D_data(var,outdir,tim1,tim2)
begin
  var1=var({time|tim1:tim2},lat|:,lon|:)
  vname     = var@vnam
  remove_unnecessary_attribute_of(var1,"_lindx")
  remove_unnecessary_attribute_of(var1,"vnam")
  if isatt(var,"filNam")
    var_name=var@filNam
    remove_unnecessary_attribute_of(var1,"filNam")
  else
    var_name=vname
  end if
  update_or_new(var1,vname)   ; check if output file is to be updated
    remove_unnecessary_attribute_of(var1,"update")
    remove_unnecessary_attribute_of(var1,"vnam")
  if vname@_update
    var_update(outdir,var_name,vname,var1)
    remove_unnecessary_attribute_of(var,"update")
    remove_unnecessary_attribute_of(vname,"_update")
  else                       ; start a new file
    remove_unnecessary_attribute_of(vname,"_update")
 
    outfil    = create_file_if_not_exists(outdir,var_name)
    fout      = addfile(outfil,"c")
    write_TLL(fout,vname,var1)
  end if
end

procedure write_data_for_each_year(var,var_name,yr1,yr2,outdir)
local nrank, nyr, outdir0, tunits, year
begin
  nrank     = var@nrank
  nyr       = yr2-yr1+1
  outdir0   = outdir
  tunits    = var&time@units
  do i = 1,nyr
    var@vnam  = var_name
    year=yr1+i-1
    outdir    = create_outdir((/outdir0,numeral_to_string(year)/))
    make_the_dir(outdir)
    start_time=ut_inv_calendar(year,1,\
                       1,0,0,0,tunits,0)
    end_time =ut_inv_calendar(year,12,\
                       31,0,0,0,tunits,0)
    remove_unnecessary_attribute_of(var,"nrank")
    remove_unnecessary_attribute_of(var,"time1")
    remove_unnecessary_attribute_of(var,"time2")
    print("Output the standardized data for "+year)
    if nrank.eq.4
      write_4D_data(var,outdir,start_time,end_time)
    else
      write_3D_data(var,outdir,start_time,end_time)
    end if
  end do
end

procedure CMS_Write(var,var_name,dataRoot,
           dataSrc,dataFrq,dataTyp)
begin
;  init_write_library()
  check_rank_of(var) 
  check_time_object_of(var)
  if var@nrank.eq.4  
    check_lev(var)
  end if
  dataFrq   = capitalize(dataFrq)
  var@DataFreq = dataFrq
  make_standard_time_axis_for(var)
  ;printVarSummary(var)

  dataSrc   = changeCase(dataSrc,"up")
  dataTyp   = capitalize(dataTyp)
  outdir    = create_outdir((/dataRoot,dataSrc,dataFrq,dataTyp/))
  remove_xtra_slashes_of(outdir)

  D1        = var@time1
  D2        = var@time2
  y1        = floattointeger(D1(0,0))
  y2        = floattointeger(D2(0,0))
  pp("Writing data to database")
  pp(var_name)
  write_data_for_each_year(var,var_name,y1,y2,outdir)
end

function CMS_FillGaps(var,fill_req)
local ndim, idim
begin

  ndim=dimsizes(fill_req)
  do idim = 0,ndim-1
    print("Filling data gaps in "+ fill_req(idim))
    fillup_along_time(var)
  end do
return(var)
end

function CMS_Anom(var,vclim,y1,y2)
local var1, iyr, tim1, tim2, tindx
begin
  vanom=var
  tunits=var&time@units
  t0=0
  do iyr=y1,y2  
    tim1=ut_inv_calendar(iyr,1,1,0,0,0,tunits,0)
    tim2 =ut_inv_calendar(iyr,12,31,0,0,0,tunits,0)
    var1 = var({tim1:tim2},:,:)

  ; Find day of year for the daily data set
  ; From time variable calculate year, month and day
    time=var1&time
    tarray=ut_calendar(time,0)

    if var@DataFreq.eq."Monthly"
      tindx=floattoint(tarray(:,1))
    end if
    if var@DataFreq.eq."Daily"
      tindx=floattoint(tarray(:,2))
    end if
    tindx=tindx-1

    time!0="time"
    time&time=time
    vtims=time({tim1:tim2})
    nvtims=dimsizes(vtims)
    if dimsizes(tindx).eq.1
      vanom(t0:t0+nvtims-1,:,:) = var1- (/ vclim(tindx:tindx,:,:) /)
    else
      vanom(t0:t0+nvtims-1,:,:) = var1- (/ vclim(tindx,:,:) /)
    end if
    t0=t0+nvtims

    delete(var1)
    delete(time)
    delete(tarray)
    delete(vtims)
    delete(tindx)
  end do
return(vanom)
end

procedure Make_CMS_Anomaly(opt)
local opt
begin
  allowed_freq=(/"Monthly","Daily"/)
  if .not.any(allowed_freq.eq.opt@DataFreq)
    terminate_program(opt@DataFreq+" anomaly not supported")
  end if

  ;cr=integertochar(10) ; already defined
  if opt@DataFreq.eq."Daily"
    if .not. isatt(opt,"NumHarmonics")
      msg="You need to specify the number of harmonics to smooth the data"
      msg=msg+cr+"        Please use opt@NumHarmonics=(1, 2 or 3)"
      terminate_program(msg)
    end if
  end if

  print("Retrieving requested data for anomaly construction")
  opt@DataType="Mean"
  var=slice_data(opt)
  ;printVarSummary(var)
  print("Retrieving climatology data")
  opt@DataType="Climatology"
  vclim=slice_data(opt)
  ;printVarSummary(vclim)

  time=var&time
  ntim=dimsizes(time)
  D1=get_yr_mo_dy(time(0))
  D2=get_yr_mo_dy(time(ntim-1))
  y1=D1@year
  y2=D2@year
  m1=D1@month
  m2=D2@month
  d1=D1@day
  d2=D2@day

  dataFrq   = capitalize(opt@DataFreq)
  opt@DataFreq=dataFrq ; for some reason after capitalization, this option
                       ; is downcased
  var@DataFreq=dataFrq
  if dataFrq.eq."Monthly"
    ;vanom=calcMonAnomTLL (var,vclim)
    vanom=CMS_Anom (var,vclim,y1,y2)
  end if
  if dataFrq.eq."Daily"
    ddd      = day_of_year(y1, m1, d1)
    yyyyddd  = y1*1000 + ddd
    vclim_sm = smthClmDayTLL(vclim, opt@NumHarmonics) 
    ;vanom    = calcDayAnomTLL (var, yyyyddd, vclim_sm)
    vanom=CMS_Anom (var,vclim_sm,y1,y2)
  end if

  dataSrc   = changeCase(opt@DataSource,"up")
  dataTyp   = "Anomaly"
  outdir    = create_outdir((/opt@DataRoot,dataSrc,dataFrq,dataTyp/))
  remove_xtra_slashes_of(outdir)

  delete(var)
  delete(vclim)
  var=vanom
  delete(vanom)
  var_name=opt@vname
  if isatt(opt,"Lev")
   var_fil=var_name+"_"+opt@Lev
   var@filNam=var_fil
  end if
  
  check_rank_of(var) 
  if isatt(opt,"update")
    var@update=opt@update
  end if
  if isatt(opt,"FillGaps")
    var=CMS_FillGaps(var,opt@FillGaps)
  end if
  write_data_for_each_year(var,var_name,y1,y2,outdir)
end

procedure Make_CMS_Climatology(opt)
begin
  allowed_freq=(/"Monthly","Daily"/)

  if .not.any(allowed_freq.eq.opt@DataFreq)
    terminate_program(opt@DataFreq+" climatology not supported")
  end if

  ;cr=integertochar(10)
  if .not. isatt(opt,"ClimPeriod")
    msg="You need to specify the period for climatology calculation"
    msg=msg+cr+"        Please use the opt@ClimPeriod=(/period1,period2/)"
    terminate_program(msg)
  end if

  opt@Time=opt@ClimPeriod
  print("Retrieving requested data for climate construction")
  var=slice_data(opt)
  tunits="days since 1900-01-01"

  print("Make Climatology")
  if opt@DataFreq.eq."Monthly"
    clim = clmMonTLL(var)
    std  = stdMonTLL(var)
    ntim=12
    time=make_monthly_time(1904,1,15,1904,12,15,tunits,ntim)
  end if

  if opt@DataFreq.eq."Daily" 
    time    = var&time
    TIME    = ut_calendar(time, 0)          ; type float 
    delete(time)
    year    = floattointeger( TIME(:,0) )
    month   = floattointeger( TIME(:,1) )
    day     = floattointeger( TIME(:,2) ) 
    ddd     = day_of_year(year, month, day) 
    yyyyddd = year*1000 + ddd
    clim=clmDayTLL(var,yyyyddd) 
    ntim=366
    time=make_daily_time(1904,1,15,1904,12,15,tunits,ntim)
  end if
  CPeriod=opt@ClimPeriod
  CPeriod@join=" "
  clim@ClimPeriod=array_to_string(CPeriod)

  dataSrc   = changeCase(opt@DataSource,"up")
  dataFrq   = capitalize(opt@DataFreq)
  dataTyp   = "Climatology"
  outdir    = create_outdir((/opt@DataRoot,dataSrc,dataFrq,dataTyp/))
  remove_xtra_slashes_of(outdir)
  print(outdir)
  dims=dimsizes(clim)
  ntim=dims(0)
  nlat=dims(1)
  nlon=dims(2)
  nlev=1
  clim2=new( (/ntim,nlev,nlat,nlon/), typeof(clim) )
  clim2(:,0,:,:)=clim
  ;clim2(:,1,:,:)=clim ; till we develop code to write std deviation
  delete(clim2&$clim2!0$)
  ;printVarSummary(clim2)
  clim2!0="time"
  clim2&time=time
  clim2!1="level"
  clim2&level=ispan(1,nlev,1)
  clim2&level@description="climatology"
  vname=opt@vname
  var_name=vname
  if isatt(opt,"Lev")
   var_name=var_name+"_"+opt@Lev
  end if
  make_the_dir(outdir)
  outfil    = create_file_if_not_exists(outdir,var_name)
  fout      = addfile(outfil,"c")
  remove_unnecessary_attribute_of(clim2,"nrank")
  remove_unnecessary_attribute_of(clim2,"time1")
  remove_unnecessary_attribute_of(clim2,"time2")
  remove_unnecessary_attribute_of(clim2,"vnam")
  remove_unnecessary_attribute_of(clim2,"info")
  remove_unnecessary_attribute_of(clim2,"time_op_ncl")
  remove_unnecessary_attribute_of(clim2,"valid_range")
  remove_unnecessary_attribute_of(clim2,"actual_range")
  remove_unnecessary_attribute_of(clim2,"precision")
  if isatt(opt,"FillGaps")
    clim2=CMS_FillGaps(clim2,opt@FillGaps)
  end if
  print("Write Climatology to DataBase")
  write_TLLL(fout,vname,clim2)
end
