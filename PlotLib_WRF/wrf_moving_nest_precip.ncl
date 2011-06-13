load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "./WRFUserARW_add.ncl"

begin

  wks = gsn_open_wks("x11", "moving_nest_precip")  ; Open graphics file

  root="/mnt/disk1/work/WRF/wrf-typhoon/"
  dir="out/2010/8/9/"
  FILES=systemfunc("ls -1 "+root+dir+"wrfout_d01_2010-08*")+".nc"

  numFILES = dimsizes(FILES)
  print("numFILES = " + numFILES)
  print(FILES)

  f = addfiles(FILES,"r")

  var_in = wrf_user_getvar_from_files(FILES,"RAINNC",-1,True)
  lat = wrf_user_getvar_from_files(FILES,"XLAT",-1,True)
  lon = wrf_user_getvar_from_files(FILES,"XLONG",-1,True)
  times = wrf_user_getvar_from_files(FILES,"times",-1,True)

  dom_dims = dimsizes(var_in)
  dom_rank = dimsizes(dom_dims)
  nx = dom_dims(dom_rank - 1) - 1
  ny = dom_dims(dom_rank - 2) - 1

  res                 = True

; Set some contouring resources.
  res@cnFillOn        = True
  res@cnLinesOn       = False
  res@cnLineLabelsOn  = False
  res@cnInfoLabelOn   = False
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = (/ .1, .2, .4, .8, 1.6, 3.2, 6.4, \
                               12.8, 25.6, 51.2, 102.4/)
  res@cnFillColors         = (/"White","White","DarkOliveGreen1", \
                               "DarkOliveGreen3","Chartreuse", \
                               "Chartreuse3","Green","ForestGreen", \
                               "Yellow","Orange","Red","Violet"/)
; Add map resources
  res@mpDataBaseVersion     = "MediumRes"          ; Default is LowRes
  res@mpOutlineDrawOrder    = "PostDraw"           ; Draw map outlines last
  res@mpGridAndLimbOn       = False                ; Turn off lat/lon lines
  res@pmTickMarkDisplayMode = "Always"             ; Turn on map tickmarks
  res = set_mp_wrf_map_resources(f[0],res)
  res@mpLimitMode        = "Corners"               ; Portion of map to zoom
  res@mpLeftCornerLatF   = min (lat) - 2.0
  res@mpLeftCornerLonF   = min (lon) - 5.0
  res@mpRightCornerLatF  = max (lat) + 2.0
  res@mpRightCornerLonF  = max (lon) + 5.0

; Add label bar resources
  res@lbLabelAutoStride = True 
  res@lbBoxMinorExtentF        = 0.13
  res@lbLabelFontHeightF       = 0.012


  res@gsnLeftStringFontHeightF = 0.01
  res@gsnRightStringFontHeightF = 0.008
  res@gsnMaximize     = True    ; Maximize plot in frame

; set up arrays to we can figure out how the nest is moving
  px  = new ( numFILES , integer )
  py  = new ( numFILES , integer )
  x1s = new ( numFILES , integer )
  x1e = new ( numFILES , integer )
  y1s = new ( numFILES , integer )
  y1e = new ( numFILES , integer )
  x2s = new ( numFILES , integer )
  x2e = new ( numFILES , integer )
  y2s = new ( numFILES , integer )
  y2e = new ( numFILES , integer )

  px(:) = f[:]@I_PARENT_START
  py(:) = f[:]@J_PARENT_START

; calculate the movement of the nests
  do ip = 1,numFILES-1
    dx = (px(ip)-px(ip-1))*3
    dy = (py(ip)-py(ip-1))*3
    if ( dx .ge. 0 ) then
      x1s(ip) = 0 + abs(dx)
      x1e(ip) = nx
      x2s(ip) = 0
      x2e(ip) = nx - x1s(ip) 
    else
      x2s(ip) = 0 + abs(dx)
      x2e(ip) = nx
      x1s(ip) = 0
      x1e(ip) = nx - x2s(ip) 
    end if
    if ( dy .ge. 0 ) then
      y1s(ip) = 0 + abs(dy)
      y1e(ip) = ny
      y2s(ip) = 0
      y2e(ip) = ny - y1s(ip) 
    else
      y2s(ip) = 1 + abs(dy)
      y2e(ip) = ny
      y1s(ip) = 1
      y1e(ip) = ny - y2s(ip) 
    end if
  end do

; For each time period, get the tendency of the rainfall in the overlapping regions
; As the nests move very slowly we can assume the overlay is big, so we can
; discard the information in no-overlapping regions
; BE CAREFUL - if domains do not overlay over large regions, this assumption is wrong
  do ip = 1,numFILES-1      ; only interested in time 1 omwards, as time 0 should have 0 rainfall

    print("Working on time : " + times(ip) )
    cnx = (x2e(ip)-x2s(ip))+1    ; the area we are plotting will be this big
    cny = (y2e(ip)-y2s(ip))+1
    var = new ( (/cny,cnx/) , typeof(var_in) )
    var = 0.0
    var@description = "1hr rainfall tendency ("+ times(ip) +")"
    var@units = "mm"
    ix1 = x1s(ip)-1
    ix2 = x2s(ip)-1
    do ix = 0,cnx-1          
      ix1 = ix1 + 1
      ix2 = ix2 + 1
      iy1 = y1s(ip)-1
      iy2 = y2s(ip)-1
      do iy = 0,cny-1          
        iy1 = iy1 + 1
        iy2 = iy2 + 1
        var(iy,ix) = var_in(ip,iy2,ix2) - var_in(ip-1,iy1,ix1)
      end do
    end do
    var@lat2d = lat(ip,y2s(ip):y2e(ip),x2s(ip):x2e(ip))
    var@lon2d = lon(ip,y2s(ip):y2e(ip),x2s(ip):x2e(ip))
    map = gsn_csm_contour_map(wks, var, res)
    delete(var)
  
  end do

end
