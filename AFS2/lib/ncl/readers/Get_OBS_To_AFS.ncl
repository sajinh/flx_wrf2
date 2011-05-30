load "$AFS/CODE/COMMON/NCL/writem.ncl"
load "$AFS/CODE/COMMON/NCL/AFS.Definitions"
load "$AFS/CODE/COMMON/NCL/CMS_Retrieve.ncl"

function ODir()
begin
;return("http://"+210.98.49.102:8080/opendap/")
return("/apcc01/OPER/AFS/DATA/OBS/")
end

procedure AFS_WriteOBS(var)
begin
ObsDir="/apcc01/OPER/AFS/DATA/OBS/"
outfil=ObsDir+var@vname+".nc"
print(outfil)
if isfilepresent(outfil)
  system("rm "+outfil)
end if
fout=addfile(outfil,"c")
  write_TLL(fout,var@vname,var)
end

procedure AFS_CatOBS(var)
begin
ObsDir="/apcc01/OPER/AFS/DATA/OBS/"

; Our DataBase File or Target File
infil=ObsDir+var@vname+".nc"

; This is our newly retrieved data
; we will cat it to the Target File
; For now we store it in the temporary
; directory and will delete after
; we finish ncrcat

tmpfil="./"+var@vname+".nc"
ftmp=addfile(tmpfil,"c")
write_TLL(ftmp,var@vname,var)

system("ncrcat -h -O "+infil+" "+tmpfil+" "+infil)
system("rm "+tmpfil)
end

function lookUp(vname)
begin
  if vname.eq."z500"
    datasource="NCEP"
    datasource@var="hgt"
    datasource@lev="500"
  end if
  if vname.eq."u200"
    datasource="NCEP"
    datasource@var="uwnd"
    datasource@lev="200"
  end if
  if vname.eq."v200"
    datasource="NCEP"
    datasource@var="vwnd"
    datasource@lev="200"
  end if
  if vname.eq."u850"
    datasource="NCEP"
    datasource@var="uwnd"
    datasource@lev="850"
  end if
  if vname.eq."v850"
    datasource="NCEP"
    datasource@var="vwnd"
    datasource@lev="850"
  end if
  if vname.eq."t850"
    datasource="NCEP"
    datasource@var="air"
    datasource@lev="850"
  end if
  if vname.eq."ts"
    datasource="NCEP_SFC"
    datasource@var="air"
  end if
  if vname.eq."slp"
    datasource="NCEP_SFC"
    datasource@var="slp"
  end if
  if vname.eq."prec"
    datasource="GPCP"
    datasource@var="prec"
  end if

  if .not.isdefined("datasource")
  print("This data is not available for download yet")
  exit
  end if
  datasource@url="210.98.49.102"
  return(datasource)
end

function getLogFile(vname)
begin
  printVarSummary(vname)
  if vname@Access.eq."Local" 
    LogFile="http://"+DataSource@url+":8080/opendap/"
    LogFile=LogFile+DataSource+"/Daily/Logger/"+vnam2+"_first_time.nc"
  else
    LogFile="http://"+DataSource@url+":8080/opendap/"
    LogFile=LogFile+DataSource+"/Daily/Logger/"+vnam2+"_first_time.nc"
  end if

end

procedure OBS_Seeder(vname)
; Finds first time from CMS data base
; Loads first month data corresponding to first time
; from CMS
begin

  ;writeGraphVizLog("lookUp")
  DataSource=lookUp(vname)
  ; Now get firsttime.nc
  ; do a opendap call to CMS data server

  vnam2=DataSource@var
  if isatt(DataSource,"lev")
    vnam2=DataSource@var+"_"+DataSource@lev
  end if

  LogFile=getLogFile(vname)
  fin=addfile(LogFile,"r")
 
  time=fin->time
  date=ut_calendar(time,0) 
  year=doubletoint(date(:,0))
  mo=doubletoint(date(:,1))
  dy=doubletoint(date(:,2))
  if dy.ne.1
   print("We do not have enough data: eXit")
   exit
  end if

  month=integer_to_month(mo)
  dys=days_in_month(year,mo)

  opt=True
  opt@vname=DataSource@var
  opt@DataSource=DataSource
  opt@Access=vname@Access
  
  tim1="01"+month+year
  tim2=dys+month+year

  opt@Time=(/tim1,tim2/)
  opt@TimAve=True
  opt@Lat=(/-90,90/)
  opt@Lon=(/0,360/)

  if isatt(DataSource,"lev")
  opt@Lev=DataSource@lev
  end if


opt@DEBUG=True
opt@DataType="Mean"

;writeGraphVizLog("slice_data")
var1=slice_data(opt)
nlat=dimsizes(var1&lat)
nlon=dimsizes(var1&lon)
var=new( (/1,nlat,nlon/),typeof(var1))
add_attributes(var,(/"time","lat","lon"/))
tunits=AFS_TimeUnits()
time=ut_inv_calendar(year,mo,15,0,0,0,tunits,0)
var&time=time
var&lat=var1&lat
var&lon=var1&lon
var(0,:,:)=var1
var@vname=vname
;writeGraphVizLog("AFS_WriteOBS")
AFS_WriteOBS(var)
end

procedure OBS_Update(vname)
; Finds first time from CMS data base
; Loads first month data corresponding to first time
; from CMS
begin

  DataSource=lookUp(vname)
  ; Now get firsttime.nc
  ; do a opendap call to CMS data server

  vnam2=DataSource@var
  if isatt(DataSource,"lev")
    vnam2=DataSource@var+"_"+DataSource@lev
  end if

; Check last available time in AFS DataBase
  ObsDir=ODir()
  fin=addfile(ObsDir+vname+".nc","r")
  time=fin->time
  ntim=dimsizes(time)
  time0=time(ntim-1) ; the last available time
  delete(time)
  date=ut_calendar(time0,0) 
; Last Year, Month and Day in AFS DataBase

  yrT=doubletoint(date(:,0))
  moT=doubletoint(date(:,1))
  dyT=doubletoint(date(:,2))
  delete(date)

; Check last available time in Source Repository
  LogFile="http://"+DataSource@url+":8080/opendap/"
;  LogFile=getenv("AFS")
  LogFile=LogFile+"/"+DataSource+"/Daily/Logger/"+vnam2+"_last_time.nc"
  fin=addfile(LogFile,"r")
print(fin)
 
exit
  time=fin->time
  date=ut_calendar(time,0) 
; Last Year, Month and Day in Source Repository

  yrS=doubletoint(date(:,0))
  moS=doubletoint(date(:,1))
  dyS=doubletoint(date(:,2))

  nyr=yrS-yrT                             

  nmons=(nyr-1)*12 + (12-moT) + moS ; No of months to retrieve

; if the last month does not have enough daily data to make monthly
; means ignore that month

  requiredLastDay=days_in_month(yrS,moS)
  if dyS.lt.requiredLastDay
  nmons=nmons-1
  end if

; Now do a loop over 'nmons', retrieve data and append it to
; target (AFS DataBase for OBS)

  mo=moT
  year=yrT

do imon=0,nmons-1
  mo=mo+1
  if mo.gt.12
    mo=1
    year=year+1
  end if
  month=integer_to_month(mo)
  print((imon+1)+") Getting Data for "+month+" "+year)

  dys=days_in_month(year,mo)

  opt=True
  opt@vname=DataSource@var
  opt@DataSource=DataSource
  
  tim1="01"+month+year
  tim2=dys+month+year

  opt@Time=(/tim1,tim2/)
  opt@TimAve=True
  opt@Lat=(/-90,90/)
  opt@Lon=(/0,360/)

  if isatt(DataSource,"lev")
  opt@Lev=DataSource@lev
  end if


  opt@DEBUG=True
  opt@DataType="Mean"
  opt@Access=vname@Access

  var1=slice_data(opt)
  nlat=dimsizes(var1&lat)
  nlon=dimsizes(var1&lon)
  var=new( (/1,nlat,nlon/),typeof(var1))
  add_attributes(var,(/"time","lat","lon"/))
  tunits=AFS_TimeUnits()
  time=ut_inv_calendar(year,mo,15,0,0,0,tunits,0)
  var&time=time
  var&lat=var1&lat
  var&lon=var1&lon
  var(0,:,:)=var1
  var@vname=vname
  AFS_CatOBS(var)

end do

end

doWhat=getenv("doWhat")
vname=getenv("var")
vname@Access="Remote" ; use the OpenDAP Server
;vname@Access="Local" ; Default gets from /apcc02

if doWhat.eq."seed"
  OBS_Seeder(vname)
else
  OBS_Update(vname)
end if

