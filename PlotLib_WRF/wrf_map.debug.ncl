load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
;
; Read in data.
;
  cdf_file = addfile("wrf_map.debug.nc","r")

;
; Open X11 file.
;
  wks = gsn_open_wks("x11", "wrf_map.debug")
;  gsn_define_colormap(wks, )
res                             = True
  res@mpFillOn                    = False        ; turn off gray fill
  res@mpOutlineBoundarySets       = "National"   ; turn on country boundaries
  res@mpGeophysicalLineColor      = "Navy"       ; color of cont. outlines
  res@mpGeophysicalLineThicknessF = 1.5          ; thickness of outlines
  res@mpRightCornerLonF = 167.366
  res@mpRightCornerLatF = 39.7568
  res@mpLeftCornerLonF = 114.723
  res@mpLeftCornerLatF = 5.57973
  res@mpLimitMode = "Corners"
  ;res@mpCenterRotF = 65.06
  res@mpCenterRotF = 24.94
  res@mpCenterLonF = 180

  plot = gsn_csm_map_ce(wks,res)                 ; draw map

exit

  res = True
  res@mpProjection = "CylindricalEquidistant"
  ;res@gsnMaximize = True
;  res@gsnFrame = False
;  res@gsnDraw = False
  ;res@mpRightCornerLonF = 167.366
  ;res@mpRightCornerLatF = 39.7568
  ;res@mpLeftCornerLonF = 114.723
  ;res@mpLeftCornerLatF = 5.57973
  ;res@mpLimitMode = "Corners"
  ;res@mpGridSpacingF = 45
  ;res@mpCenterRotF = 24.94
  ;res@mpCenterRotF = 65.06
  ;res@mpCenterLonF = 180
  res@mpFillOn = True
res@mpOutlineBoundarySets       = "National"   ; turn on country boundaries
  res@mpGeophysicalLineColor      = "Navy"       ; color of cont. outlines
  res@mpGeophysicalLineThicknessF = 1.5     

  plot = gsn_csm_map_ce(wks, res)
end
