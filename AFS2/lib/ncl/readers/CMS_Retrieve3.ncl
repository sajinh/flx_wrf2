TimeUnits="days since 1900-01-01"

function findCurrentTime()
begin
  cyr=stringtoint(systemfunc("date +%Y"))
  cmo=stringtoint(systemfunc("date +%m"))
  cdy=stringtoint(systemfunc("date +%d"))
return( (/cyr,cmo,cdy/) )
end


function cut_data(DataDir,vname,yrs,tim,lat,lon)
local ntim,nlat,nlon,lat1,lat2,tim1,tim2,lon1,lon2
begin
  ntim=dimsizes(tim)
  nlat=dimsizes(lat)
  nlon=dimsizes(lon)
  var2=new((/ntim,nlat,nlon/),"float")

  lat1=lat(0)
  lon1=lon(0)
  tim1=tim(0)

  lat2=lat(nlat-1)
  lon2=lon(nlon-1)
  tim2=tim(ntim-1)

  nyrs=dimsizes(yrs)
    t0=0
    do iyr=0,nyrs-1
    DataDir_yr=DataDir+yrs(iyr)+"/"+vname+".nc"
      fin=addfile(DataDir_yr,"r")
      var1=fin->$vname$({tim1:tim2},{lat1:lat2},{lon1:lon2})
      tims=var1&time
      tims!0="tims"
      tims&tims=tims
      vtims=tims({tim1:tim2})
      nvtims=dimsizes(vtims)
      var2(t0:t0+nvtims-1,:,:)=var1
      delete(var1)
      delete(tims)
      delete(vtims)
      t0=t0+nvtims
    end do

  var_dims=(/ntim,nlat,nlon/)
  dim_names=(/"time","lat","lon"/)
  new_indcs=ind(var_dims.gt.1) ; chose indices of the
			     ; the input array for
                             ; which dimsize > 1
  newdims=var_dims(new_indcs)
  new_dimnames=(dim_names(new_indcs))
  ndims=dimsizes(newdims)
  newvar=new(newdims,"float")
  newvar=var2
  do i=0,ndims-1
    newvar!i=new_dimnames(i)
  end do
return(newvar)
end

function cut_level_data(DataDir,vname,yrs,tim,lev,lat,lon)
local ntim,nlev,nlat,nlon,lat1,lat2,tim1,tim2,lon1,lon2,nstat
begin
  ntim=dimsizes(tim)
  nlev=dimsizes(lev)
  nlat=dimsizes(lat)
  nlon=dimsizes(lon)
  var2=new((/ntim,nlev,nlat,nlon/),"float")

  lat1=lat(0)
  lon1=lon(0)
  tim1=tim(0)

  lat2=lat(nlat-1)
  lon2=lon(nlon-1)
  tim2=tim(ntim-1)

  nyrs=dimsizes(yrs)
    t0=0
    do iyr=0,nyrs-1
      if typeof(yrs).eq."string"
        DataDir_yr=DataDir+"/"+vname+"_"
        nstat=vname@nstat
      else
        DataDir_yr=DataDir+yrs(iyr)+"/"+vname+"_"
       print("Now processing "+yrs(iyr)+" from "+DataDir_yr)
      end if
      do ilev=0,nlev-1
      fin=addfile(DataDir_yr+lev(ilev)+".nc","r")
      if typeof(yrs).eq."string"
        var1=fin->$vname$({tim1:tim2},nstat,{lat1:lat2},{lon1:lon2})
      else
        var1=fin->$vname$({tim1:tim2},{lat1:lat2},{lon1:lon2})
      end if
      tims=var1&time
      tims!0="tims"
      tims&tims=tims
      vtims=tims({tim1:tim2})
      nvtims=dimsizes(vtims)
      var2(t0:t0+nvtims-1,ilev,:,:)=var1
      delete(var1)
      delete(tims)
      delete(vtims)
      end do
      t0=t0+nvtims
    end do

  var_dims=(/ntim,nlev,nlat,nlon/)
  dim_names=(/"time","level","lat","lon"/)
  new_indcs=ind(var_dims.gt.1) ; chose indices of the
			     ; the input array for
                             ; which dimsize > 1
  newdims=var_dims(new_indcs)
  new_dimnames=(dim_names(new_indcs))
  ndims=dimsizes(newdims)
  newvar=new(newdims,"float")
  newvar=var2
  do i=0,ndims-1
    newvar!i=new_dimnames(i)
  end do
return(newvar)
end

function retrieve_lats(fin,lat1,lat2)
begin
  lat=fin&lat
  lat!0="lat"
  lat&lat=lat
  tmp_lat=lat({lat1:lat2})
return(tmp_lat)
end

function retrieve_lons(fin,lon1,lon2)
begin
  lon=fin&lon
  lon!0="lon"
  lon&lon=lon
  tmp_lon=lon({lon1:lon2})
return(tmp_lon)
end

function retrieve_tims(fin)
begin
  return(fin&time)
end

function find_datalevs(DataDir,vname,lev1,lev2)
begin

if vname@Access.eq."Remote"
  if DataDir@url.eq."210.98.49.102"
    CMD="rsync "+DataDir@url+"::cms/Data_Org/"+DataDir+vname+"_*.nc"
  else
    CMD="rsync "+DataDir@url+"::cms/"+DataDir+vname+"_*.nc"
  end if
  junk=systemfunc(CMD)
  njnk=dimsizes(junk)

  ; Below is an algorithm to detect meaningful
  ; info, in this case, the years
  
  jind=0 ; by default
  do ijnk=0,njnk-1
  jchar=stringtochar(junk(ijnk))
  njchar=dimsizes(jchar)
    if njchar.gt.1 
      if (jchar(njchar-3:njchar-2).eq."nc")
        jind=ijnk
        break
      end if
    end if
  delete(jchar)
  end do

  file_levs1=junk(jind:)
  lflevs=dimsizes(file_levs1)
  file_levs=new(lflevs,float) ;  i am using float to prevent 
                            ; mistaking it as index
  do ilev=0,lflevs-1
    CMD="echo "+file_levs1(ilev)+" | awk  '{print $5}'"
    CMD=CMD+" | awk -F"+vname+"_ '{print $2}'"
    file_levs(ilev)=stringtofloat(systemfunc(CMD))
  end do

else
  CMD="ls "+DataDir+vname+"_*.nc "
  CMD=CMD+" | awk -F"+DataDir+vname+"_ '{print $2}'"
  CMD=CMD+" | awk -F. '{print $1}'"
  file_levs=stringtofloat(systemfunc(CMD))
end if

  qsort(file_levs)

  file_levs!0="levs"
  file_levs&levs=file_levs
  my_levs=file_levs({lev1:lev2})
  return(my_levs)
end

procedure check_timespecs(fin1,fin2,start_time,end_time)

begin
  time_of_first_data=fin1&time(0)
  tmp_time=fin2&time
  time_of_last_data=tmp_time(dimsizes(tmp_time)-1)

    if time_of_last_data.lt.end_time
      print("Data request exceeds data collection :")
      print("Last date available is  "+ut_calendar(time_of_last_data,-2))
      print("You requested data upto "+ut_calendar(end_time,-2))
      print("Please revise request and try again : eXiT")
      exit
    end if

    if time_of_first_data.gt.start_time
      print("Data request exceeds data collection :")
      print("First date available is  "+ut_calendar(time_of_first_data,-2))
      print("You requested data from  "+ut_calendar(start_time,-2))
      print("Please revise request and try again : eXiT")
      exit
    end if
end

;function make_newvar(var, var_dims )
;begin
;  dim_names=(/"time","level","lat","lon"/)
;  new_indcs=ind(var_dims.gt.1) ; chose indices of the
;			     ; the input array for
;                             ; which dimsize > 1
;  newdims=var_dims(new_indcs)
;  new_dimnames=(dim_names(new_indcs))
;  ndims=dimsizes(newdims)
;  newvar=new(newdims,"float")
;  newvar=var
;  do i=0,ndims-1
;    newvar!i=new_dimnames(i)
;  end do
;return(newvar)
;end

undef("integer_value_of_month")
;------------------------------------------------------------------------------
function integer_value_of_month(month[1]:string)
;------------------------------------------------------------------------------
begin

  months=(/"JAN","FEB","MAR","APR","MAY","JUN", \
              "JUL", "AUG","SEP","OCT","NOV","DEC"/)
  return(ind(months.eq.month)+1)

end




function CMS_RetrieveData(org_DataDir,vname,domain)
begin
DataDir=org_DataDir

; The default settings.
; We discover the ranges of these 4 dimensions
; by parsing the user requests.

nlat=0
nlon=0
nlev=0
nyrs=0

if isatt(domain,"lat1")
  lat1=domain@lat1
  lat2=domain@lat2
 if lat1.eq.lat2
 nlat=1 ; if nlat=0 we should exit, at least one latitude should be specified
 else
 nlat=2
 end if
end if

if isatt(domain,"lon1")
  lon1=domain@lon1 ; nlon should not be zero, if zero we exit
  lon2=domain@lon2
 if lon1.eq.lon2
 nlon=1
 else
 nlon=2
 end if
end if

nlev=0 ; Default
  vnam2=vname
if isatt(domain,"lev1")
  lev1=stringtofloat(domain@lev1)
  lev2=stringtofloat(domain@lev2)
  nlev=1 ; nlev can be zero, because some data do not have level info
  vnam2=vname+"_"+lev1
end if

DataDir=DataDir+"/"+domain@DataFreq+"/"+domain@DataType+"/"

if vname@Access.eq."Local" .and. .not.isfilepresent(DataDir)
  print("The directory "+DataDir+" is not available")
  TerminateProgram
end if

if vname@Access.eq."Remote"
  DataUrl=DataDir@url
  if DataDir@url.eq."210.98.49.102"
    CMD="rsync "+DataDir@url+"::cms/Data_Org/"+DataDir
  else
    CMD="rsync "+DataDir@url+"::cms/"+DataDir
  end if
  junk=systemfunc(CMD)
  njnk=dimsizes(junk)

  ; Below is an algorithm to detect meaningful 
  ; info, in this case, the years
  
  jind=0
  do ijnk=0,njnk-1
  jchar=stringtochar(junk(ijnk))
  njchar=dimsizes(jchar)
    if njchar.gt.1 
      if (jchar(njchar-2).eq.".")
        jind=ijnk
        break
      end if
    end if
  delete(jchar)
  end do

  ; Here we are *NOT* checking for a specific file
  ; hence the first entry is likely to be a
  ; file pointing to current directory

  jind=jind+1

  file_yrs1=junk(jind:)
  lfyrs=dimsizes(file_yrs1)
  file_yrs=new(lfyrs,integer)
  do iyr=0,lfyrs-1
    file_yrs(iyr)=stringtoint(systemfunc("echo "+file_yrs1(iyr)+" | awk  '{print $5}'"))
  end do
else
  if domain@DataType.eq."Climatology"
    file_yrs="0"
    nstat=domain@nstat
  else
    file_yrs=stringtoint(systemfunc("ls "+DataDir))
  end if
end if


if typeof(file_yrs).eq."integer"
  tim1=stringtochar(domain@tim1)
  tim2=stringtochar(domain@tim2)
  nsiz=sizeof(tim1)-1

  if nsiz.lt.8 .or. nsiz.gt.9
    print("Error in Time specification")
    print("Time to be specified as ddMONyyyy, e.g. 20Jan2006")
    print("However current time is provided as "+domain@tim1)
    print("eXiT")
    exit
  end if

  yr1=stringtoint(chartostring(tim1(nsiz-4:nsiz-1)))
  smo1=changeCase(chartostring(tim1(nsiz-7:nsiz-5)),"up")
  dy1=stringtoint(chartostring(tim1(:nsiz-8)))


  nsiz=sizeof(tim2)-1
  yr2=stringtoint(chartostring(tim2(nsiz-4:nsiz-1)))
  smo2=changeCase(chartostring(tim2(nsiz-7:nsiz-5)),"up")
  dy2=stringtoint(chartostring(tim2(:nsiz-8)))

  mo1=integer_value_of_month(smo1)
  mo2=integer_value_of_month(smo2)

  tunits=TimeUnits ; retrieved from top-level Defaults file

  if yr1.lt. min(file_yrs) .or. yr2.gt.max(file_yrs)
  print("Prescribed Year not to be found in Data Base")
  print("This could be either "+yr1+" or "+yr2)
  exit
  else
  nyrs=yr2-yr1+1
  my_yrs=yr1+nyrs-1
  end if

; Check for consistent specification

  if nlat.eq.0 .or. nlon.eq.0 .or. nyrs.eq.0
  print("At least one latitude, longitude and time point to be specified")
  print("Current specs are given below: ")
  printVarSummary(domain)
  exit
  end if
end if

  if nlat.eq.0 .or. nlon.eq.0 
  print("At least one latitude and longitude point to be specified")
  print("Current specs are given below: ")
  printVarSummary(domain)
  exit
  end if

_datFreq=domain@DataFreq
allowed_data_freqs=(/"Daily","Monthly","Yearly"/)

; we define monthly data to be located at mid-month
if _datFreq.eq."Monthly"
  dy1=15
  dy2=15
end if


if typeof(file_yrs).eq."integer"
  start_time=ut_inv_calendar(yr1,mo1,dy1,0,0,0,tunits,0)
  end_time=ut_inv_calendar(yr2,mo2,dy2,0,0,0,tunits,0)

  if _datFreq.eq."Yearly"
    ntim=doubletointeger(fabs(yr2-yr1+1))
    data_tim=make_yearly_time(yr1,mo1,dy1,yr2,mo2,dy2,tunits,ntim)
  end if
  if _datFreq.eq."Monthly"
    ntim=doubletointeger(fabs( (yr2-yr1-1)*12+(12-mo1+1)+mo2))
    data_tim=make_monthly_time(yr1,mo1,dy1,yr2,mo2,dy2,tunits,ntim)
  end if
  if _datFreq.eq."Daily"
    ntim=doubletointeger(end_time-start_time+1)
    data_tim=make_daily_time(yr1,mo1,dy1,yr2,mo2,dy2,tunits,ntim)
  end if


  if .not.any(allowed_data_freqs.eq._datFreq)
    terminate_program("DataFreq of "+_datFreq+" not allowed")
  end if

  data_tim!0="time"
  data_tim&time=data_tim
end if
; Parse the time and find starting and ending time
 
if typeof(file_yrs).eq."integer"
  DataDir_y1=DataDir+yr1+"/"
else
  DataDir_y1=DataDir+"/"
end if

if vname@Access.eq."Remote"
  DataDir_y1@url=DataDir@url
  DataDir="http://"+DataDir@url+":8080/opendap/"+DataDir
end if

if nlev.ge.1  		; we are looking at Level Data
  ; Find if data for specified level is available
  data_lev=find_datalevs(DataDir_y1,vname,lev1,lev2)
  nlev=dimsizes(data_lev)
  
  if typeof(file_yrs).eq."integer"
    fin1=addfile(DataDir+yr1+"/"+vname+"_"+data_lev(0)+".nc","r")
    fin2=addfile(DataDir+yr2+"/"+vname+"_"+data_lev(0)+".nc","r")
    check_timespecs(fin1,fin2,start_time,end_time)
  else
    fin1=addfile(DataDir+"/"+vname+"_"+data_lev(0)+".nc","r")
  end if

  data_lat=retrieve_lats(fin1,lat1,lat2)
  data_lon=retrieve_lons(fin1,lon1,lon2)
  nlat=dimsizes(data_lat)
  nlon=dimsizes(data_lon)


  if typeof(file_yrs).eq."integer"
    years=ispan(yr1,yr2,1)
    var=cut_level_data(DataDir,vname,years,data_tim,data_lev,data_lat,data_lon)
  else
    data_tim=retrieve_tims(fin1)
    years="null"
    vname@nstat=nstat
    var=cut_level_data(DataDir,vname,years,data_tim,data_lev,data_lat,data_lon)
  end if
else
  if typeof(file_yrs).eq."integer"
    fin1=addfile(DataDir+yr1+"/"+vname+".nc","r")
    fin2=addfile(DataDir+yr2+"/"+vname+".nc","r")
    check_timespecs(fin1,fin2,start_time,end_time)
  else
    fin1=addfile(DataDir+"/"+vname+".nc","r")
  end if
  data_lat=retrieve_lats(fin1,lat1,lat2)
  data_lon=retrieve_lons(fin1,lon1,lon2)
  nlat=dimsizes(data_lat)
  nlon=dimsizes(data_lon)

  if typeof(file_yrs).eq."integer"
    years=ispan(yr1,yr2,1)
    var=cut_data(DataDir,vname,years,data_tim,data_lat,data_lon)
  else
    var=fin1->$vname$(:,nstat,{lat1:lat2},{lon1:lon2})
  end if
end if

return(var)
end

function CMS_FindDataDir(DataSource,vname)
; Locate DataDirectories and returns
; address of data in the DataBase
; return value is an array with string elements
begin

BaseDir="/apcc01/CIS_Data/"
if vname@Access.eq."Remote"
BaseDir=""
end if

TMI=BaseDir+"TMI/"
QUK=BaseDir+"QuikScat/"
AVS=BaseDir+"AVISO/"
JRA=BaseDir+"JRA/"
JRS=BaseDir+"JRA_SFC/"
NCP=BaseDir+"NCEP/"
OLR=BaseDir+"NOAA_OLR/"
NCS=BaseDir+"NCEP_SFC/"

if DataSource.eq."TMI"
DataDir=TMI
end if

if DataSource.eq."QUIKSCAT"
DataDir=QUK
end if

if DataSource.eq."AVISO"
DataDir=AVS
end if

if DataSource.eq."JRA"
DataDir=JRA
end if

if DataSource.eq."JRA_SFC"
DataDir=JRS
end if

if DataSource.eq."NOAA_OLR"
DataDir=OLR
end if

if DataSource.eq."NCEP"
DataDir=NCP
 if any( (/"air","slp"/).eq.vname)
   if .not.isatt(vname,"Lev")
     DataDir=NCS
   end if
 end if
end if

if DataSource.eq."NCEP_SFC"
 DataDir=NCS
end if


if vname@Access.eq."Remote"
;(Remote servers are at 190.1.1.24 and 210.98.49.102)
;DataDir@url="190.1.1.24" 
DataDir@url="210.98.49.102"
end if

return(DataDir)
end
 
procedure space_remover(val)
; remove empty spaces from a string

begin
  tmp=val
  val=systemfunc("echo "+tmp+" | tr -d [:blank:]")
end

procedure CMS_subsetter(val1,val2,vals)
; subroutine that will map a 2-element
; domain array specification into 
; the two domain boundaries
begin

  nvals=dimsizes(vals)
  val1=vals(0)
  if nvals.eq.1
  val2=val1
  else
  val2=vals(1)
  end if

  if typeof(vals).eq."string"
  space_remover(val1)
  space_remover(val2)
  end if

  if  vals@DEBUG
  print(val1+" -- "+val2)
  end if
  

end

procedure debug_checker(opt,vals)
begin

if isatt(opt,"DEBUG")
vals@DEBUG=True
else
vals@DEBUG=False
end if

end



function CMS_getdata(opt)

; Plotting routine requests data through this
; interface.
; The routine can specify
; Data Source, Data period, Data type , Data domain etc. as
; attributes to a logical variable.

; Sample usage:
; opt@Lat=(/-30,30/) ; we need domain specs because we may plot
; opt@Lon=(/40,120/) ; hovmoeller diagrams and it will be efficient
; opt@Lev=(/300,500/); to retrieve subsets instead of huge domains
; opt@Time=(/19Jan2006,20Mar2007/)
; opt@DataSource="TMI"
; opt@Variable="sst"
; opt@DataType="interannual anomalies" ; band-passed anoms at a particular day
; opt@DataType="intraseasonal anomalies"; band-passed anoms at a particular day
; opt@DataType="monthly anomalies"; raw anomalies averaged for the requested period
; opt@DataType="seasonal anomalies"; raw anomalies averaged for requested period
; opt@DataType="weekly anomalies"; raw anomalies averaged for requested period
; opt@DataType="climatology"; climatology averaged for requested period


; Our maps are 7-day, monthly and 90-day
begin

   domain=True
   ;anomaly_types=(/"monthly","seasonal","weekly"/)+"anomalies"

; Determine which DataType to use
; Account for inconsistencies of user input related to case
; and spaces.

domain@DataType="Anomaly" ; default data type 
                          ; defaulted to if nothing is specified
                          ; or opt@DataType is any of
                          ; monthly weekly or seasonal anomalies
domain@DataFreq=opt@DataFreq

if isatt(opt,"DataType")
  opt@DataType=changeCase(opt@DataType,"low")
  space_remover(opt@DataType)

  if opt@DataType.eq."mean"
    domain@DataType="Mean"
  end if

  if opt@DataType.eq."intraseasonalanomalies"
    domain@DataType="ISV"
  end if

  if opt@DataType.eq."climatology"
    domain@DataType="Climatology"
    ostat=opt@Statistic
    ostat=changeCase(ostat,"low")
    if ostat.eq."mean"
      domain@nstat=0
    end if
    if ostat.eq."stddev"
      domain@nstat=1
    end if
  end if

  if opt@DataType.eq."lanom"
    domain@DataType="LAnom"
  end if
end if

  if (opt.and.isatt(opt,"Lat")) then
  lats=opt@Lat
  debug_checker(opt,lats)
  domain@lat1=1.0
  domain@lat2=1.0
  CMS_subsetter(domain@lat1,domain@lat2,lats)
  end if

  if (opt.and.isatt(opt,"Lon")) then
  lons=opt@Lon
  domain@lon1=1.0
  domain@lon2=1.0
  debug_checker(opt,lons)
  CMS_subsetter(domain@lon1,domain@lon2,lons)
  end if

  if (opt.and.isatt(opt,"Time")) then
  tims=opt@Time
  domain@tim1="1"
  domain@tim2="1"
  debug_checker(opt,tims)
  CMS_subsetter(domain@tim1,domain@tim2,tims)
  end if

  if (opt.and.isatt(opt,"Lev")) then
  levs=opt@Lev
  domain@lev1="1"
  domain@lev2="1"
  debug_checker(opt,levs)
  CMS_subsetter(domain@lev1,domain@lev2,levs)
  end if

  opt@DataSource=changeCase(opt@DataSource,"up")
  vname=changeCase(opt@Variable,"low")
  vname@Access="Local"
  if isatt(opt,"Lev")
    vname@Lev=opt@Lev
  end if

  if (opt.and.isatt(opt,"Access")) then
    vname@Access=opt@Access
  end if

  if isatt(opt,"DataBaseDir")
    DataDir=opt@DataBaseDir+"/"+opt@DataSource+"/"
  else
    DataDir=CMS_FindDataDir(opt@DataSource,vname); returns list of directories from which
  end if
                                         ; to get data
  var=CMS_RetrieveData(DataDir,vname,domain) ; opens all the directories
                                                ; gets all the data
  return(var)
end


function slice_data(sopt)
begin

  opt=True
  opt@Variable=sopt@vname
  opt@DataSource=sopt@DataSource
  if isatt(sopt,"Time")
    opt@Time=sopt@Time
  end if
  if isatt(sopt,"Lat")
    opt@Lat=sopt@Lat
  else
    opt@Lat=(/-90,90/)
  end if
  if isatt(sopt,"Lon")
    opt@Lon=sopt@Lon
  else
    opt@Lon=(/0,360/)
  end if

; Optional Parameters
  if isatt(sopt,"Lev")
    opt@Lev=sopt@Lev
  end if

  if isatt(sopt,"DataType")
    opt@DataType=sopt@DataType
    if opt@DataType.eq."Climatology"
      if isatt(sopt,"Statistic")
        opt@Statistic=sopt@Statistic
      else
        opt@Statistic="Mean"
      end if
    end if
  else
    opt@DataType="Anomalies"
  end if

  if isatt(sopt,"Access")
    opt@Access=sopt@Access
  end if

  if isatt(sopt,"DataFreq")
    opt@DataFreq=sopt@DataFreq
  else
    opt@DataFreq="Daily"
  end if

  if isatt(sopt,"DataRoot")
    opt@DataBaseDir=sopt@DataRoot
  end if

; Data is retrieved using a complex algorithm
  vdat=CMS_getdata(opt)

; Find no of dimensions to average
nave=0
dstring=""

if isatt(sopt,"TimAve")
  nave=nave+1
  dstring=dstring+"tim"
end if

if isatt(sopt,"LevAve")
  nave=nave+1
  dstring=dstring+"lev"
end if

if isatt(sopt,"LatAve")
  nave=nave+1
  dstring=dstring+"lat"
end if

if isatt(sopt,"LonAve")
  nave=nave+1
  dstring=dstring+"lon"
end if


; Strategy 1: Assume that the request is correct 
;             a - Find the dimensions to average
;             b - Find the number of dimensions of retrieved data
;             c - create a string array to hold the re-ordered dimensions
;             d - at the end of the new array, put the dimensions we have to 
;                 average
;             e - for the other positions, do an iterative search


;Find the dimensions we have to average

; dstring consists of dimension names
; specified as 3-letter words, such as lev, lat, lon and tim

if nave.gt.0
  cdstring=stringtochar(dstring)
  ndstring=nave ; no of dimensions to average


  my_dstring=new(ndstring,"string")

  ; create an array my_dstring to sort out the dimensions
  ; in the right order
  ik=0
  do is=0,ndstring-1
  my_dstring(is)=chartostring(cdstring(ik:ik+2))
  ik=ik+3
  end do

  Std_dimensions=(/"time","level","lat","lon"/)
  Abr_dimensions=(/"tim","lev","lat","lon"/)


  DimsToAve=new(ndstring,string) ; Map into real dimension names

  do is=0,ndstring-1
    do ks=0,3
       DimsToAve(is)=Std_dimensions(ks)
       if Abr_dimensions(ks).eq.my_dstring(is)
         break
       end if
    end do
  end do

;  print(DimsToAve)


; Find the dimensions of the retrieved data
  data_dims=getvardims(vdat)
  nrank=dimsizes(data_dims)

ndiff=nrank-nave
if ndiff.le.0
print("You are asking for the ridiculous : eXiT")
exit
end if

ReorderedDims=new(nrank,"string")
ReorderedDims(nrank-ndstring:)=DimsToAve


; Now fill up the remaining of the ReorderedDims with
; the appropriate dimensions

ncount=ndiff
do it=0,ndiff-1
  do is=0,nrank-1
       if any( data_dims(is).eq.ReorderedDims(ncount:) )
       continue
       end if
       ReorderedDims(ncount-1)=data_dims(is)
  end do
ncount=ncount-1
end do


;if ndiff.lt.1 .or. ndiff.gt.2 
;  print("We need either 2-D or 1-D data")
;  print("The resulting dimension will be "+ndiff)
;  print("eXiT")
;  exit
;end if


; How to re-order co-ordinate variables
; and then to use the right averaging routine
; if nave=1, use dim_avg
; if nave=2, use wgt_areaave
; if nave=3, use wgt_volave
; dims=(/"time","lat","lon"/)
; var=vdat($dims(1)$|:,$dims(2)$|:,$dims(0)$|:)
; printVarSummary(var)

  dims=ReorderedDims
  if nrank.eq.4
  var1=vdat($dims(0)$|:,$dims(1)$|:,$dims(2)$|:,$dims(3)$|:)
  end if

  if nrank.eq.3
  var1=vdat($dims(0)$|:,$dims(1)$|:,$dims(2)$|:)
  end if

  if nrank.eq.2
  var1=vdat($dims(0)$|:,$dims(1)$|:)
  end if

  if nrank.eq.1
  var=vdat
  end if

    if nave.eq.3
    var=wgt_volave(var1,1,1,1,0)
    end if

    if nave.eq.2
    var=wgt_areaave_Wrap(var1,1,1,0)
    end if

    if nave.eq.1
    var=dim_avg_Wrap(var1)
    end if
 
else
var=vdat
end if

return(var)
end

function findFirstTimeOfCollection(DataDir)
begin
  years=systemfunc("ls "+DataDir)
  year=stringtoint(years(0))
  ;print("The year "+year+" is the first year")

  MeanDir=DataDir+"/"+year
  files=systemfunc("ls "+MeanDir+"/*.nc")
  fin=addfile(files(0),"r")
  torig=fin->time
  ntim0=dimsizes(torig)
  tfirst=torig(0)
  tnarr=ut_calendar(tfirst,0)
  y0=floattoint ( tnarr(0,0) )
  m0=floattoint ( tnarr(0,1) )
  d0=floattoint ( tnarr(0,2) )
  firstTime=(/y0,m0,d0/)
  firstTime@join=","
  return(array_to_string(firstTime))
end
function findLastTimeOfCollection(DataDir)
begin
  years=systemfunc("ls "+DataDir)
  nyrs = dimsizes(years)
  year=stringtoint(years(nyrs-1))
  ;print("The year "+year+" is the last year")

  MeanDir=DataDir+"/"+year
  files=systemfunc("ls "+MeanDir+"/*.nc")
  fin=addfile(files(0),"r")
  torig=fin->time
  ntim0=dimsizes(torig)
  tlast=torig(ntim0-1)

  tnarr=ut_calendar(tlast,0)
  y0=floattoint ( tnarr(0,0) )
  m0=floattoint ( tnarr(0,1) )
  d0=floattoint ( tnarr(0,2) )
  lastTime=(/y0,m0,d0/)
  lastTime@join=","
  return(array_to_string(lastTime))
end

function getTime(opt)
	;opt@vname
	;opt@DataSource
	;opt@DataRoot
	;opt@DataFreq="Monthly"
	;opt@DataType="Mean"
begin
  if .not. isatt(opt,"DataSource")
    terminate_program("You need to specify the source of the data")
  end if
  if .not. isatt(opt,"WhichTime")
    terminate_program("You need to specify which Time to return")
  end if
  if isatt(opt,"DataRoot")
    opt@DataBaseDir=opt@DataRoot
    remove_unnecessary_attribute_of(opt,"DataRoot")
  end if
  if isatt(opt,"DataBaseDir")
    DataDir=opt@DataBaseDir+"/"+opt@DataSource+"/"
  else
    DataDir=CMS_FindDataDir(opt@DataSource,vname); 
  end if
  if .not.isatt(opt,"DataFreq")
    opt@DataFreq="Daily"
  end if
  if .not.isatt(opt,"DataType")
    opt@DataType="Mean"
  end if
  DataDir=DataDir+"/"+opt@DataFreq+"/"+opt@DataType
  if changeCase(opt@WhichTime,"low").eq."first"
   return(findFirstTimeOfCollection(DataDir))
  end if
  if  changeCase(opt@WhichTime,"low").eq."last"
    return(findLastTimeOfCollection(DataDir))
  end if
end
