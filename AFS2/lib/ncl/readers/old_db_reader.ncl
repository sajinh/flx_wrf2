procedure add_units_and_long_name(var,var_name)
begin
  if (var_name.eq."prec") then
    var@units="mm/day"
    var@long_name="Precipitation"
  end if
  if (var_name.eq."slp") then
        var@units="mb"
        var@long_name="Sea Level Pressure"
  end if
  if (var_name.eq."sst") then
    var@units="deg K"
    var@long_name="Sea Surface Temperature"
  end if
  if (var_name.eq."t2m") then
    var@units="deg K"
    var@long_name="Temperature at 2m"
  end if
  if (var_name.eq."t850") then
    var@units="deg K"
    var@long_name="Air temperature at 850mb"
  end if
  if (var_name.eq."ts") then
    var@units="deg K"
    var@long_name="Temperature at the surface layer/level"
  end if
  if (var_name.eq."u200") then
    var@units="m/s"
    var@long_name="Zonal wind at 200mb"
  end if
  if (var_name.eq."u850") then
    var@units="m/s"
    var@long_name="Zonal wind at 850mb"
  end if
  if (var_name.eq."v200") then
    var@units="m/s"
    var@long_name="Meridional wind at 200mb"
  end if
  if (var_name.eq."v850") then
    var@units="m/s"
    var@long_name="Meridional wind at 850mb"
  end if
  if (var_name.eq."z500") then
    var@units="m"
    var@long_name="Geopotential height at 500 mb"
  end if
end

function binary_reader(dir,yr,smon,lead,nlat,nlon,model,vname)
begin
  file_prefx=dir+"/"+model+"/total/"+vname+"/"+vname
  infil=file_prefx+yr+sprinti("%0.2i",smon)
  if isfilepresent(infil+"_1m.gdat") then
    nensembles=dimsizes(systemfunc("ls "+infil+"_*m.gdat"))-1
  else
    afs_abort(infil+"_1m.gdat not present")
  end if

  var=new((/lead,nensembles,nlat,nlon/), "float",AFS_FillValue() )
  y1=yr
  ilead=0
  emon=smon+lead-1
  do imon=smon,emon
    if imon.gt.12
      imon2=imon-12
      y1=y1+1
    else
      imon2=imon
    end if
    do iens=1,nensembles
      infil=file_prefx+y1+sprinti("%0.2i",imon2)+"_"+iens+"m.gdat"
      if isfilepresent(infil) then
        var(ilead,iens-1,:,:)=fbindirread(infil,0,(/nlat,nlon/),"float")
      else
        afs_abort(infil+" not present")
      end if
    end do
  ilead=ilead+1
  end do
  var@ensembles=ispan(1,nensembles,1)
  add_units_and_long_name(var,vname)
  var@missing_value = var@_FillValue
  return(var)
end
