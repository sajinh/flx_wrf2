load "$AFS2/lib/ncl/helper_libs.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

wrf_dir="../wrf-oper/test_out/2010/7/20/"
wrf_dir="../wrf-demo/run/"
do idom = 1,1
  f = addfile(wrf_dir+"wrfout_d0"+idom+"_2010-07-20_12:00:00.nc","r")            
  ofil="wrf_dom_"+idom+".temp"
  type = "x11"
  wks = gsn_open_wks(type,ofil)
  gsn_define_colormap(wks,"temp_19lev")   ; overwrite the .hluresfile color map
  times  = wrf_user_list_times(f)  ; get times in the file
  ntimes = dimsizes(times)         ; number of times in the file
  var = wrf_user_getvar(f,"T2",ntimes-1) 
  var = var - 273.16
print(max(var))
print(min(var))
  
  lat2d  = f->XLAT(0,:,:)     
  lon2d  = f->XLONG(0,:,:) 
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
  var_r!0="lat"
  var_r!1="lon"
  var_r&lat=lat
  var_r&lon=lon  
                                                         
  res                       = True        ; plot mods desired
  res@gsnMaximize           = True        ; uncomment to maximize size
  res@cnLineThicknessF      = 1.0         ; make twice as thick as default
  res@cnLineLabelBackgroundColor = "background" ; no lines through labels
  res@gsnAddCyclic          = False       ; regional data: not cyclic 
  res@mpOutlineDrawOrder    = "PostDraw"  ; draw continental outline last
  res@mpFillDrawOrder       = "Draw"
  res@mpFillOn              = True
;************************************************
; establish title 
;************************************************
  res@tiMainString          = ""
  nodraw(res)
  noframe(res)
  opt=True
  opt@latlon=(/lat1,lat2,lon1,lon2/)
  mymap(res,opt)
  ;opt@clevs=fspan(25,39,19)
  myfill(res,opt)
  res@gsnSpreadColors  = True
  res@lbLabelBarOn    = True
  res@cnSmoothingOn = False
  res@cnFillMode = "RasterFill"
  res@lbLabelAutoStride = True
  res@cnFillMode = "CellFill"
  res@mpProjection          = "LambertConformal"
  res@gsnMaskLambertConformal = True            ; turn on lc masking
  if ntimes.gt.1
  res@cnCellFillEdgeColor           = -1
  else
  res@cnCellFillEdgeColor           = 1
  end if
  res@cnCellFillMissingValEdgeColor = "-1"
  res@cnFillDrawOrder = "Predraw"

;************************************************
; Turn off the fill boundaries, since we are going
; to explicitly set which boundaries we want to
; fill.
;************************************************
fill_specs = (/"water"/)
mask_specs = (/"Japan"/)
 ;res@mpGridMaskMode       = "MaskLand"  ; 
  res@mpFillBoundarySets    = "NoBoundaries"
  res@mpDataBaseVersion     = "Ncarg4_1"
  res@mpFillAreaSpecifiers  = fill_specs
  res@mpAreaMaskingOn       = 1
  res@mpMaskAreaSpecifiers  = mask_specs
  res@mpSpecifiedFillColors = (/7,0/)

;************************************************
; Turn on masking, and set resource indicating
; areas to mask.
;************************************************
  res@mpDataSetName = "Earth..4"

  plot                      = gsn_csm_contour_map(wks, var_r, res)

     ; attributes of the polylines 
     resp = True 
     resp@gsLineColor      = "blue" ; polylines color 
     resp@gsLineThicknessF = 1.0    ; polylines thickness 


  draw(plot)
  frame(wks)
  delete(wks)
  delete(times)
  delete(var)
  delete(var_r)
  delete(lat2d)
  delete(lon2d)
  delete(lon)
  delete(lat)

end do
end

