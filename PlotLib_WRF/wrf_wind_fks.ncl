; Overlay information from 2 domains
; November 2009

load "$AFS2/lib/ncl/helper_libs.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$WRF/draw_pref.ncl"
load "$WRF/wrf_helpers.ncl"


procedure make_plot(idom,yr,mo,dy)
begin
; plot specs
  ;--- by default wind at 10 meters

  least_lead=0
  uwnd="ua" ; variable to plot
  vwnd="va" ; variable to plot
  voffset=0.0         ; offset to apply before plotting
  dom="d"+sprinti("%02i",idom)
  plot_type="x11"
; ---

  mon=sprinti("%02i",mo)
  day=sprinti("%02i",dy)
  hr="18:00:20"

  root="/home/saji/WRF-FKS/"
  dir="out/"+yr+"/"+mo+"/"+dy
 
  plot_root=getenv("WRF")+"/WRF-FKS-PLOT-Test/"

  day="14"
  d1 = addfile(root+dir+"/wrfout_"+dom+"_"+yr+"-"+mon+"-"+day+"_"+hr+".nc", "r")
;wrfout_d03_2011-03-14_18:00:20

  tim_opt=True
  times = wrf_user_getvar(d1,"times",-1)
  ntimes = dimsizes(times)

  do ilead=0,22,2
    print(ilead)
    req_time=  ilead
    u = wrf_user_getvar(d1,uwnd,req_time) ; - voffset
    v = wrf_user_getvar(d1,vwnd,req_time) ; - voffset
    u1 = u(0,:,:)
    v1 = v(0,:,:)
    time=reformat_wrf_time(times(req_time),tim_opt)

    hora=str_split(hr,":")
    plot_dir = plot_root+"init."+yr+"-"+mo+"-"+dy+"-"+hora(0)+"/"+time@yr+"/"+time@mo+"/"
    plot_dir = plot_dir+time@dy+"/"+time@hr+"/"
    if plot_type .eq. "pdf"
      system("mkdir -p "+plot_dir)
    end if
    plot_file = plot_dir+"wind_"+dom

    lat1 = wrf_user_getvar(d1,"XLAT",0)
    lon1 = wrf_user_getvar(d1,"XLONG",0)

    u1@lat2d = lat1
    u1@lon2d = lon1
    v1@lat2d = lat1
    v1@lon2d = lon1

    dom_dims = dimsizes(u1)
    dom_rank = dimsizes(dom_dims)
    nx1 = dom_dims(dom_rank - 1) - 1
    ny1 = dom_dims(dom_rank - 2) - 1

    speed1    = u1;
    speed1    = sqrt(u1^2+v1^2)
    speed1@long_name = "Wind Speed"

    wks = gsn_open_wks(plot_type, plot_file)  ; Open graphics file
    colors = (/ "white", "black", \
                "white", "yellow", "PowderBlue", "lightseagreen",\
                "PaleGreen","wheat", "Brown", "Pink" /)
  clevs = (/0,0.5,1,1.5,2,3,4,5/)
    gsn_define_colormap(wks,"rainbow")
    res                 = True


    res@vcRefMagnitudeF         = 10.0              ; make vectors larger
    res@vcRefLengthF            = 0.040            ; ref vec length
    res@vcGlyphStyle            = "CurlyVector"    ; turn on curly vectors
    ;res@vcGlyphStyle            = "FillArrow"    ; turn on curly vectors
    res@vcMonoFillArrowFillColor = True
    ;res@vcFillArrowEdgeColor = -1
    res@vcMinDistanceF          = 0.017            ; thin out vectors
    res@vcVectorDrawOrder       = "PostDraw"
    res@vcLineArrowHeadMinSizeF = 0.008
    res@vcLineArrowHeadMaxSizeF = 0.008

    res@tiMainString            = "Wind vectors colored by wind speed"
    res@gsnLeftString = ""
    res@gsnRightString = ""

    res@gsnDraw         = False
    res@gsnFrame        = False


  ; Add map resources
    res@mpDataBaseVersion     = "HighRes"          ; Default is LowRes
    res@mpOutlineDrawOrder    = "PostDraw"           ; Draw map outlines last
    res@mpGridAndLimbOn       = False                ; Turn off lat/lon lines
    res@pmTickMarkDisplayMode = "Always"             ; Turn on map tickmarks
    res = set_mp_wrf_map_resources(d1,res)
    res@mpLimitMode        = "Corners"               ; Portion of map to zoom
    res@mpLeftCornerLatF   = lat1(0,0)
    res@mpLeftCornerLonF   = lon1(0,0)
    res@mpRightCornerLatF  = lat1(ny1,nx1)
    res@mpRightCornerLonF  = lon1(ny1,nx1)
    res@mpFillOn = True
    res@mpOutlineDrawOrder    = "PostDraw"           ; Draw map outlines last
    res@mpOceanFillColor        = 15
    res@mpInlandWaterFillColor = 2


  ; Add label bar resources
    ;res@gsnMaximize     = True    ; Maximize plot in frame
    res@gsnLeftString = time



  ; we need these to later draw boxes for the location of the nest domain
    xbox_out = new(5,float)
    ybox_out = new(5,float)
    lnres = True
    lnres@gsLineThicknessF  = 1.5

  ; make images
    nodraw(res)
    noframe(res)
    res@vpYF = 0.8
    ;map = gsn_csm_vector_map(wks, u1, v1, res)
    myfill(res,True)
    res@lbLabelBarOn = True; 
    res@lbLabelAutoStride = True
    map = gsn_csm_contour_map(wks, u1, res)

    res3 = True
    res3@gsnSpreadColors     = True      ; use full colormap
;    res3@gsnSpreadColorStart = 6         ; first color to use
;    res3@gsnSpreadColorEnd   = 196       ; last color to use
    opt = True
    opt@clevs = clevs
    myfill(res3,opt)
    nodraw(res3)
    noframe(res3)
    res3@lbLabelBarOn = True; 
    res3@lbLabelAutoStride = True
    ;res3@gsnMaximize  = True  ; Use maximization from original plot
    res3@gsnLeftString = ""
    wsp = gsn_csm_contour(wks,u1,res3)
    ;overlay(map,wsp)

  ; PLOT 3
    draw(map)   ; domain 2 already overlaid here - so just draw again

  ; Add a mark where AizuWakamatsu is located
    aizu_lat = 37.498
    aizu_lon = 139.926 

    pmres = True
    pmres@gsMarkerColor = "Red"
    pmres@gsMarkerIndex = 16
    pmres@gsMarkerSizeF = 0.005
    gsn_polymarker(wks,map,aizu_lon,aizu_lat,pmres)

    if idom .le. 2
      map@resolution = "prefectures"
    end if
    draw_pref_on_plot_v2(wks,map)
    map@resolution = "prefectures"
    draw_pref_on_plot_v2(wks,map)

    frame(wks)

    plot_file=plot_file+".pdf"
    pp(plot_file)
    if isfilepresent(plot_file)
      pp("SUCCESS")
    end if

    delete(wks)
  end do
end

do idom=3,3
  print("Here")
  yr= 2011; stringtointeger(getenv("WFK_YR"))
  mo=    3; stringtointeger(getenv("WFK_MO"))
  dy=   12; stringtointeger(getenv("WFK_DY"))
  make_plot(idom,yr,mo,dy)
end do
exit

