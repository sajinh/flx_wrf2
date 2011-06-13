
load "$AFS2/lib/ncl/helper_libs.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$WRF/draw_pref.ncl"
load "$WRF/wrf_helpers.ncl"
load "$WRF/gsn_color.ncl"

procedure set_col_style(res,clevs,wks)
begin
  grid_min=min(clevs)
  grid_max=max(clevs)
  nlevs=dimsizes(clevs)
  grid_step=(grid_max-grid_min)/(1.0*nlevs)
  grid_center = 0.0
  color_map = gsn_retrieve_colormap(wks)
  ncolors = dimsizes(color_map)
  color_end = ncolors(0)
  trans_color = new(1,"integer")
  if isatt(wks,"trans_color")
    trans_color = wks@trans_color
  else
    trans_color = ((color_end - 2)/2)
  end if
  gsn_range = gsnColorRange(grid_min, grid_max, grid_step, grid_center, color_end, trans_color)
end


procedure set_lbar(res)
begin
  res@lbLabelBarOn = True
  res@lbLabelAutoStride = True
end

procedure set_con_style(res,clevs)
begin
  opt=True
  opt@clevs=clevs
  myfill(res,opt)
end

procedure set_vec_style(res)
begin
  res@vcRefMagnitudeF         = 10.0              ; make vectors larger
  res@vcRefLengthF            = 0.050            ; ref vec length
  res@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors
  res@vcMinDistanceF          = 0.017            ; thin out vectors
end

procedure set_map_style(res,fin)
begin
  res@mpDataBaseVersion     = "HighRes"          ; Default is LowRes
  res@mpOutlineDrawOrder    = "PostDraw"           ; Draw map outlines last
  res@mpGridAndLimbOn       = False                ; Turn off lat/lon lines
  res@pmTickMarkDisplayMode = "Always"             ; Turn on map tickmarks
  res = set_mp_wrf_map_resources(fin,res)
  res@mpFillOn = True
  res@mpOutlineDrawOrder    = "PostDraw"           ; Draw map outlines last
  res@mpOceanFillColor        = 15
  res@mpInlandWaterFillColor = 2
end

procedure set_vwport_from_latlon(res,lat2d,lon2d)
begin
  dom_dims = dimsizes(lat2d)
  dom_rank = dimsizes(dom_dims)
  nx1 = dom_dims(dom_rank - 1) - 1
  ny1 = dom_dims(dom_rank - 2) - 1

  res@mpLimitMode        = "Corners"               ; Portion of map to zoom
  res@mpLeftCornerLatF   = lat2d(0,0)
  res@mpLeftCornerLonF   = lon2d(0,0)
  res@mpRightCornerLatF  = lat2d(ny1,nx1)
  res@mpRightCornerLonF  = lon2d(ny1,nx1)
end

wrf_root = "/home/saji/WRF-FKS/out/2011/3/12"
fname=wrf_root+"/wrfout_d03_2011-03-14_18:00:20"
fin=fopen_read(fname+".nc")
topt=True
times = wrf_user_getvar(fin,"times",-1)

do ptime=1,23
;ptime=3
  time=reformat_wrf_time(times(ptime),True)
u = wrf_user_getvar(fin,"U10",ptime)
v = wrf_user_getvar(fin,"V10",ptime)
lat2d = wrf_user_getvar(fin,"XLAT",0)
lon2d = wrf_user_getvar(fin,"XLONG",0)
u@lat2d=lat2d
u@lon2d=lon2d

wks = gsn_open_wks("x11", "test")  ; Open graphics file
color_map="BlWhRe"
gsn_define_colormap(wks, color_map)


res=True
nodraw(res)
noframe(res)
res2=res

clevs = ispan(-10,10,1)
set_con_style(res2,clevs)
set_col_style(res2,clevs,wks)
set_lbar(res2)
set_vwport_from_latlon(res2,lat2d,lon2d)
set_map_style(res2,fin)
setstring(res2,0.01,"","")

con=gsn_csm_contour_map(wks,u,res2)

setstring(res,0.01,time,"wind")
set_vec_style(res)
vec=gsn_csm_vector(wks,u,v,res)

overlay(con,vec)

draw(con)

con@resolution="prefectures"
draw_pref_on_plot_v2(wks,con)
frame(wks)
delete(wks)
delete(res)
end do
