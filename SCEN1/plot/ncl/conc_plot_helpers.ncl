procedure read_flxout_dims(nx,ny,nz,nt)
begin
  nx=288
  ny=285
  nz=18
  nt=11
end

procedure mark_point(wks,plot,lat1, lon1)
begin
  pmres = True
  pmres@gsMarkerColor = "Red"
  pmres@gsMarkerIndex = 13
  pmres@gsMarkerSizeF = 0.01
  gsn_polymarker(wks,plot,lon1,lat1,pmres)
end

procedure set_plot_resources(res,true_lat, true_lon, lat1,lon1,lat2,lon2)
begin
  res@cnFillOn               = True               ; color fill  
  res@cnLinesOn              = False              ; no contour lines
  res@cnLineLabelsOn         = False              ; no contour labels
  res@cnInfoLabelOn          = False              ; no contour info label
  
  res@mpGridLineDashPattern  = 2                  ; lat/lon lines as dashed
  res@pmTickMarkDisplayMode  = "Always"           ; turn on tickmarks
  
  res@tiMainFontHeightF      = 0.020              ; smaller title
  
  res@gsnAddCyclic           = False              ; regional data
; !!!!! any plot of data that is on a native grid, must use the "corners"
; method of zooming in on map.
  
  res@mpLimitMode       = "Corners"            ; choose range of map
  res@mpLeftCornerLatF  = lat1
  res@mpLeftCornerLonF  = lon1
  res@mpRightCornerLatF = lat2
  res@mpRightCornerLonF = lon2
  
; The following 4 pieces of information are REQUIRED to properly display
; data on a native lambert conformal grid.  This data should be specified
; somewhere in the model itself.
  
  res@mpProjection        = "LambertConformal"
  res@mpLambertParallel1F = true_lat
  res@mpLambertParallel2F = true_lat
  res@mpLambertMeridianF  = true_lon
  
; usually, when data is placed onto a map, it is TRANSFORMED to the specified
; projection. Since this model is already on a native lambert conformal grid,
; we want to turn OFF the tranformation.
  
  res@tfDoNDCOverlay = True
  res@lbLabelAutoStride = True
  res@lbBoxMinorExtentF        = 0.13
  res@lbLabelFontHeightF       = 0.012

  
   ; Add map resources
    res@mpDataBaseVersion     = "HighRes"          ; Default is LowRes
    res@mpOutlineDrawOrder    = "PostDraw"           ; Draw map outlines last
    res@mpGridAndLimbOn       = False                ; Turn off lat/lon lines
    res@pmTickMarkDisplayMode = "Always"             ; Turn on map tickmarks

  res@gsnMaximize = True
  noframe(res);
end

function flx_levs(ilev)
local levs
begin
 levs=(/15,37,66,110,267,395,579,774,974,1178,1386,1600,1818,2041,2270,\
        2504, 2744, 2989, 3242, 3501, 3767, 4040, 4322, 4612, 4911, \
        5219, 5537, 5867, 6208, 6563, 6932, 7317, 7719, 8141, 8583, \
        9541, 10064, 10682, 11242, 11877, 12600, 13409 /)
  return(levs(ilev)) 
end

function read_conc(vname,path,nx1,ny1,ngx1,ngx2,ngy1,ngy2,ilev)
begin
  nx=36
  ny=68
  nz=43
  nt=17
  read_flxout_dims(nx,ny,nz,nt)
  if vname .eq. "conc"
    conc0 = cbinread( path,  (/nt,nz,ny,nx/), "float")
    conc = conc0(:,ilev,ngy1:(ny1-ngy2),ngx1:(nx1-ngx2))
  else
    conc0 = cbinread( path,  (/nt,ny,nx/), "float")
    conc = conc0(:,ngy1:(ny1-ngy2),ngx1:(nx1-ngx2))
  end if
  mconc=asciiread("maxconc.dat",1,"float")

; Scale Concentrations
  conc=100*conc/mconc
  conc@nt=nt
  conc@nx=nx
  conc@ny=ny
  conc@nz=nz
  return(conc)
end
