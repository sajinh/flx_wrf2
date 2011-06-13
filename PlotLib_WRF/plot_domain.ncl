;   Example script display location of model domains
;   Only works for ARW domains
;   Beta code used
;   November 2008

load "$AFS2/lib/ncl/helper_libs.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "./draw_pref.ncl"
begin
;

; approx acc. 10 m with 4 digits
np1_lat= 37.4230
np1_lon=141.0329
aizu_lat = 37.498
aizu_lon = 139.926
ref_lat   = 37.827
ref_lon   = 140.253
max_dom=4

; We generate plots, but what kind do we prefer?
  type = "x11"
 type = "x11"
; type = "ps"
; type = "ncgm"
  wks = gsn_open_wks(type,"wrf_flex")

  colors = (/"white","black","White","ForestGreen","DeepSkyBlue"/)
  gsn_define_colormap(wks, colors)  

  mpres = True
  ;mpres@mpFillOn = True
  ;mpres@mpFillColors  = (/"background","DeepSkyBlue","ForestGreen","DeepSkyBlue", "transparent"/)
  mpres@mpGeophysicalLineColor      = "Black"
  mpres@mpGridLineColor             = -1
  mpres@mpLimbLineColor             = "Black"
  mpres@mpNationalLineColor         = "Black"
  mpres@mpPerimLineColor            = "Black"
  mpres@mpUSStateLineColor          = "Black"

  lnres = True 
  lnres@gsLineThicknessF = 2.5

  txres = True
  txres@txFont = "helvetica-bold"
  txres@txJust = "BottomLeft"
  txres@txPerimOn = False
  txres@txFontHeightF = 0.015
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  mpres@max_dom           = max_dom
  mpres@parent_id         = (/   1,   1,   2,   3/)
  mpres@parent_grid_ratio = (/   1,   3,   3,   3/)
  mpres@e_we              = (/  70, 64, 160, 340/)
  mpres@e_sn              = (/  60, 55, 130, 241/)
  mpres@i_parent_start    = (/   1,  24, 6, 15/)
  mpres@j_parent_start    = (/   1,  18, 7, 15/)
  mpres@dx                = 15000
  mpres@dy                = 15000
  mpres@map_proj          = "lambert"
  mpres@ref_lat           =  ref_lat
  mpres@ref_lon           =  ref_lon
  mpres@truelat1          =  ref_lat
  mpres@truelat2          =  ref_lat
  mpres@stand_lon         =  ref_lon

  mpres@mpDataBaseVersion = "HighRes"

  
  mp = wrf_wps_dom (wks,mpres,lnres,txres)


; Add a mark where DC is located
  pmres = True
  pmres@gsMarkerColor = "Blue"
  pmres@gsMarkerIndex = 16
  pmres@gsMarkerSizeF = 0.01
  gsn_polymarker(wks,mp,aizu_lon,aizu_lat,pmres)
  gsn_polyline(wks,mp,(/aizu_lon-10,aizu_lon+10/),(/aizu_lat,aizu_lat/),pmres)
  gsn_polyline(wks,mp,(/aizu_lon,aizu_lon/),(/aizu_lat-5,aizu_lat+5/),pmres)

  pmres@gsMarkerColor = "Red"
  pmres@gsMarkerIndex = 18
  gsn_polymarker(wks,mp,np1_lon,np1_lat,pmres)

  draw_pref_on_plot_v2(wks,mp)
  frame(wks)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

end

