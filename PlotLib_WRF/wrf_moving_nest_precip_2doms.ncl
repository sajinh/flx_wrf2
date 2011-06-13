; Script to overlay RAINNC tendencies from dom 1 and 2
; Domain 2 is a moving nest
; November 2009

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "./WRFUserARW_add.ncl"

begin

  wks = gsn_open_wks("x11", "moving_nest_precip_2dom")  ; Open graphics file

  root="/mnt/disk1/work/WRF/wrf-typhoon/"
  dir="out/2010/8/9/"
  FILES1=systemfunc("ls -1 "+root+dir+"wrfout_d01_2010-08*")+".nc"
  FILES2_0=systemfunc("ls -1 "+root+dir+"wrfout_d02_2010-08*")+".nc"
  FILES2=FILES2_0(::3)

  numFILES = dimsizes(FILES1)

  print("numFILES = " + numFILES)
  print(FILES1)
  print(" ")
  print(FILES2)

  f1 = addfiles(FILES1,"r")
  f2 = addfiles(FILES2,"r")

  var_in1 = wrf_user_getvar_from_files(FILES1,"RAINNC",-1,True)
  lat1 = wrf_user_getvar_from_files(FILES1,"XLAT",-1,True)
  lon1 = wrf_user_getvar_from_files(FILES1,"XLONG",-1,True)
  var_in2 = wrf_user_getvar_from_files(FILES2,"RAINNC",-1,True)
  lat2 = wrf_user_getvar_from_files(FILES2,"XLAT",-1,True)
  lon2 = wrf_user_getvar_from_files(FILES2,"XLONG",-1,True)
  times = wrf_user_getvar_from_files(FILES1,"times",-1,True)

  dom_dims = dimsizes(var_in2)
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

  res@gsnDraw         = False
  res@gsnFrame        = False

  res2 = res

; Add map resources
  res@mpDataBaseVersion     = "MediumRes"          ; Default is LowRes
  res@mpOutlineDrawOrder    = "PostDraw"           ; Draw map outlines last
  res@mpGridAndLimbOn       = False                ; Turn off lat/lon lines
  res@pmTickMarkDisplayMode = "Always"             ; Turn on map tickmarks
  res = set_mp_wrf_map_resources(f1[0],res)
  res@mpLimitMode        = "Corners"               ; Portion of map to zoom
  res@mpLeftCornerLatF   = min (lat1) 
  res@mpLeftCornerLonF   = min (lon1) 
  res@mpRightCornerLatF  = max (lat1) 
  res@mpRightCornerLonF  = max (lon1) 

; Add label bar resources
  res@lbLabelAutoStride = True 
  res@lbBoxMinorExtentF        = 0.13
  res@lbLabelFontHeightF       = 0.012


  res@gsnLeftStringFontHeightF = 0.01
  res@gsnRightStringFontHeightF = 0.008
  res@gsnMaximize     = True    ; Maximize plot in frame

  res2@gsnLeftString = ""
  res2@gsnRightString = ""
  res2@lbLabelBarOn = False  ; Labelbar already created in 1st plot
  res2@gsnMaximize  = False  ; Use maximization from original plot

; we need these to later draw boxes for the location of the nest domain
  xbox_out = new(5,float)
  ybox_out = new(5,float)
  lnres = True

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

  px(:) = f2[:]@I_PARENT_START
  py(:) = f2[:]@J_PARENT_START

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

    var1 = var_in1(ip,:,:) - var_in1(ip-1,:,:)
    var1@description = "3hr rainfall tendency ("+ times(ip) +")"
    var1@units = "mm"
    var1@lat2d = lat1(ip,:,:)
    var1@lon2d = lon1(ip,:,:)

    cnx = (x2e(ip)-x2s(ip))+1    ; the area we are plotting will be this big
    cny = (y2e(ip)-y2s(ip))+1
    var2 = new ( (/cny,cnx/) , typeof(var_in2) )
    var2@_FillValue = -999.
    var2 = -999.
    var2@description = "3hr rainfall tendency ("+ times(ip) +")"
    var2@units = "mm"
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
        var2(iy,ix) = var_in2(ip,iy2,ix2) - var_in2(ip-1,iy1,ix1)
      end do
    end do
    var2@lat2d = lat2(ip,y2s(ip):y2e(ip),x2s(ip):x2e(ip))
    var2@lon2d = lon2(ip,y2s(ip):y2e(ip),x2s(ip):x2e(ip))
    var2 = where(var2.eq.0. , -999. , var2 )
    map = gsn_csm_contour_map(wks, var1, res)
    plot = gsn_csm_contour(wks, var2, res2)
    overlay(map, plot)
    draw(map)

    xbox = (/lon2(ip,0,0),lon2(ip,0,nx),lon2(ip,ny,nx),lon2(ip,ny,0),lon2(ip,0,0)/)
    ybox = (/lat2(ip,0,0),lat2(ip,0,nx),lat2(ip,ny,nx),lat2(ip,ny,0),lat2(ip,0,0)/)
    datatondc(map, xbox, ybox, xbox_out, ybox_out)
    gsn_polyline_ndc(wks, xbox_out, ybox_out, lnres)
    frame(wks)
    delete(var1)
    delete(var2)
  
  end do

end
