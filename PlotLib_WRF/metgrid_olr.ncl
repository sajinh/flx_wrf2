load "$AFS2/lib/ncl/helper_libs.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

function regrid_wrf_var(fin,vname)
begin
  var = wrf_user_getvar(fin,vname,-1)
  lat2d  = fin->XLAT_M(0,:,:)
  lon2d  = fin->XLONG_M(0,:,:)
  dimll  = dimsizes(lat2d)                ; get size of dimensions
  nlat   = dimll(0)
  mlon   = dimll(1)
  lat1=(lat2d(0,0))
  lon1=(lon2d(0,0))
  lat2=(lat2d(nlat-1,mlon-1))
  lon2=(lon2d(nlat-1,mlon-1))
  lat = fspan(lat1,lat2,nlat)
  lon = fspan(lon1,lon2,mlon)
  var_r  = rcm2rgrid(lat2d,lon2d,var,lat,lon,0)
  lat!0="lat"
  lon!0="lon"
  lat&lat=lat
  lon&lon=lon
  lat@units="degrees_north"
  lon@units="degrees_east"
  var_r!0="time"
  var_r!1="lat"
  var_r!2="lon"
  var_r&lat=lat
  var_r&lon=lon
  return(var_r)
end

begin
  root="/mnt/disk1/work/wrf-demo/"
  dir="pre_process/"
  fil="met_em.d01.2010-08-07_18:00:00.nc"
  type = "x11"
  fin = fopen_read(root+dir+fil)
  wks = openwks(type,"olr","rainbow")

  mopt=True
  mopt@gsnDebugWriteFileName="wrf_map.debug"
  plot = wrf_map(wks,fin,mopt)
  draw(plot)
  frame(wks)
exit

  times = wrf_user_list_times(fin)
  ntimes = dimsizes(times)
  olr = regrid_wrf_var(fin,"PMSL")
  info(olr)

  res                       = True        ; plot mods desired
  res@gsnMaximize           = True        ; uncomment to maximize size
  res@gsnAddCyclic          = False       ; regional data: not cyclic
  res@mpOutlineDrawOrder    = "PostDraw"  ; draw continental outline last
  res@mpFillDrawOrder       = "Draw"
  res@mpFillOn              = True
  nodraw(res)
  noframe(res)
;************************************************
; establish title
;************************************************
  opt=True
  lat=olr&lat
  lon=olr&lon
  nlat=dimsizes(lat)
  nlon=dimsizes(lon)
  opt@latlon=(/lat(0),lat(nlat-1),lon(0),lon(nlon-1)/)
  mymap(res,opt)
  ;opt@clevs=  (/ 100., 340., 20./)
  myfill(res,opt)

  res@gsnSpreadColors  = True
  res@lbLabelBarOn    = True
  ;res@gsnSpreadColorStart = 51 
  res@mpDataSetName = "Earth..4"
  res@tiMainString          = "Dianmu 05W"


  do it = 0,ntimes-1                             ; start at time 2

    res@gsnLeftString = times(it)                   ; set Valid time on the plot
    plot = gsn_csm_contour_map(wks, olr(it,:,:), res)
    draw(plot)
    frame(wks)
  end do

end
