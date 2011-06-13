load "$AFS2/lib/ncl/helper_libs.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "WRFUserARW_add.ncl"

begin
  root="/mnt/disk1/work/WRF/wrf-typhoon/"
  dir="out/2010/8/9/"
  ;fil="wrfout_d01_2010-08-05_18:00:00.nc"
  ;dir="run/"
  FILES=systemfunc("ls -1 "+root+dir+"wrfout_d02_2010-08*")+".nc"
  fil=FILES(0)
  type = "x11"
  fin = fopen_read(fil)
  wks = openwks(type,"olr","wxpEnIR")
  ;wks = openwks(type,"olr","testcmap")

  ;times = wrf_user_list_times(fin)
  times = wrf_user_getvar_from_files(FILES,"times",-1,True)
  ntimes = dimsizes(times)
  ;olr = regrid_wrf_var(fin,"OLR")
  olr = wrf_user_getvar_from_files(FILES,"OLR",-1,True)
  lat = fin->XLAT(0,:,:)
  lon = fin->XLONG(0,:,:)
  dimll  = dimsizes(lat)                ; get size of dimensions
  nlat   = dimll(0)
  nlon   = dimll(1)

  res                     = True              ; Plot modes desired.
  res@gsnMaximize         = True              ; Maximize plot
  res@gsnSpreadColors     = True              ; use full range of colormap
  res@cnFillOn            = True              ; color plot desired
  res@cnLinesOn           = False
  res@lbLabelBarOn    = True
  res@gsnSpreadColorStart = 51 
  res@mpLimitMode                 = "Corners"           ; method to zoom
  res@mpLeftCornerLatF            = lat(0,0)
  res@mpLeftCornerLonF            = lon(0,0)
  res@mpRightCornerLatF           = lat(nlat-1,nlon-1)
  res@mpRightCornerLonF           = lon(nlat-1,nlon-1)
  res@mpProjection        = "LambertConformal"
  res@mpLambertParallel1F = 33.028
  res@mpLambertParallel2F = 33.028
  res@mpLambertMeridianF  = 131.009
  
; usually, when data is placed onto a map, it is TRANSFORMED to the specified
; projection. Since this model is already on a native lambert conformal grid,
; we want to turn OFF the tranformation.
  
  res@tfDoNDCOverlay = True


  info(olr)

  res@gsnAddCyclic          = False       ; regional data: not cyclic
  ;res@mpOutlineDrawOrder    = "PostDraw"  ; draw continental outline last
  ;res@mpFillDrawOrder       = "Draw"
  res@mpFillOn              = True
;************************************************
; establish title
;************************************************
  opt=True
  opt@clevs=  ispan(100, 340, 20)
  myfill(res,opt)

  res@gsnSpreadColors  = True
  res@lbLabelBarOn    = True
  res@gsnSpreadColorStart = 51 
  res@mpDataSetName = "Earth..4"
  res@tiMainString          = "Dianmu 05W"


  do it = 1,ntimes-1                             ; start at time 2

    res@gsnLeftString = times(it)                   ; set Valid time on the plot
    plot = gsn_csm_contour_map(wks, olr(it,:,:), res)
  end do

end
