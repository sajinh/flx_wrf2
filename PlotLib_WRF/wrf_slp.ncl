load "$AFS2/lib/ncl/helper_libs.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "WRFUserARW_add.ncl"

begin
  root="/mnt/disk1/work/WRF/wrf-typhoon/"
  ;dir="aug6-8_out/2010/8/6/"
  ;fil="wrfout_d01_2010-08-05_18:00:00.nc"
  dir="run/"
  FILES=systemfunc("ls -1 "+root+dir+"wrfout_d01_2010-08*")+".nc"
  fil=FILES(0)
  type = "x11"
  fin = fopen_read(fil)

  wks = openwks(type,"olr","rainbow")
  times = wrf_user_getvar_from_files(FILES,"times",-1,True)
  ntimes = dimsizes(times)
  olr = wrf_user_getvar_from_files(FILES,"slp",-1,True)

  lat = fin->XLAT(0,:,:)
  lon = fin->XLONG(0,:,:)
  dimll  = dimsizes(lat)                ; get size of dimensions
  nlat   = dimll(0)
  nlon   = dimll(1)

  ;olr@lat2d = lat
  ;olr@lon2d = lon
    res                     = True              ; Plot modes desired.
    res@gsnMaximize         = True              ; Maximize plot
    res@gsnSpreadColors     = True              ; use full range of colormap
    res@cnFillOn            = True              ; color plot desired
    res@cnLinesOn           = False
  res@lbLabelBarOn    = True
  res@mpLimitMode                 = "Corners"           ; method to zoom
  res@mpLeftCornerLatF            = lat(0,0)
  res@mpLeftCornerLonF            = lon(0,0)
  res@mpRightCornerLatF           = lat(nlat-1,nlon-1)
  res@mpRightCornerLonF           = lon(nlat-1,nlon-1)
  res@mpProjection        = "LambertConformal"
  res@mpLambertParallel1F = 33.028
  res@mpLambertParallel2F = 33.028
  res@mpLambertMeridianF  = 131.009
  lat2d=lat
  lon2d=lon
  true_lat1 = (max(lat2d)+min(lat2d))*0.5
  true_lat2 = (max(lat2d)+min(lat2d))*0.5
  true_lon = (max(lon2d)+min(lon2d))*0.5
  ;res@mpLambertParallel1F = true_lat1
  ;res@mpLambertParallel2F = true_lat2
  ;res@mpLambertMeridianF  = true_lon
  pp(true_lat1+"::"+true_lat2+"::"+true_lon)


  
; usually, when data is placed onto a map, it is TRANSFORMED to the specified
; projection. Since this model is already on a native lambert conformal grid,
; we want to turn OFF the tranformation.
  
  res@tfDoNDCOverlay = True


  info(olr)

  res@gsnAddCyclic          = False       ; regional data: not cyclic
  res@mpFillOn              = True
  res@pmTickMarkDisplayMode = "Always"
  res@mpFillOn              = False           ; turn off map fill
  res@mpOutlineDrawOrder    = "PostDraw"      ;
;************************************************
; establish title
;************************************************
  opt=True
  ;opt@clevs=  (/ 100., 340., 20./)
  myfill(res,opt)

  res@gsnSpreadColors  = True
  res@lbLabelBarOn    = True
  ;res@gsnSpreadColorStart = 51 
  res@mpDataSetName = "Earth..4"
  res@tiMainString          = "Dianmu 05W"
  noframe(res)

  gres                  = True                ; polyline mods desired
  gres@gsLineThicknessF = 3.0                 ; line thickness




  do it = 0,ntimes-1                             ; start at time 2

    res@gsnLeftString = times(it)                   ; set Valid time on the plot
    plot = gsn_csm_contour_map(wks, olr(it,:,:), res)
    gres@gsLineColor      = "Red"               ; line color 
    gsn_polyline (wks, plot, lon2d(nlat/2,:), lat2d(nlat/2,:), gres)
  
    gres@gsLineColor      = "Blue" 
    gsn_polyline (wks, plot, lon2d(:,nlon/2), lat2d(:,nlon/2), gres)
  
    gres@gsLineColor      = "Black" 
    gsn_polyline (wks, plot, lon2d(:,0)     , lat2d(:,0), gres)
    gsn_polyline (wks, plot, lon2d(:,nlon-1), lat2d(:,nlon-1), gres)
    gsn_polyline (wks, plot, lon2d(0,:)     , lat2d(0,:), gres)
    gsn_polyline (wks, plot, lon2d(nlat-1,:), lat2d(nlat-1,:), gres)
    frame(wks)
  end do

end
