procedure mymap(res:logical,opt:logical)
begin
	if (opt.and.isatt(opt,"latlon")) then
	res@gsnAddCyclic = False
	res@mpMinLatF      = opt@latlon(0)
	res@mpMaxLatF             = opt@latlon(1)
	res@mpMinLonF             = opt@latlon(2)
	res@mpMaxLonF             = opt@latlon(3)
        res@mpCenterLonF = (res@mpMinLonF+res@mpMaxLonF)/2. 
	end if
end

;:::::::::::::::::::::::::::::::::::::::::::::::::::::::
; Resources to make line contours on top of shaded contours
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::

procedure myfillcontour(res:logical, opt:logical)


begin
res                 = True                ; make plot mods
        if (opt.and.isatt(opt,"clevs")) then
        res@cnLevelSelectionMode = "ExplicitLevels"
        res@cnLevels             = opt@clevs
        end if

res@lbLabelBarOn    = False
res@cnSmoothingOn   = True
res@cnFillOn              = True
res@cnLinesOn             = True
res@cnLineColor = "black"
res@cnLineLabelsOn        = False    ; turn off line labels
res@cnInfoLabelOn         = False
res@gsnContourNegLineDashPattern = 1 ;
res@gsnSpreadColors = True

end
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::

;:::::::::::::::::::::::::::::::::::::::::::::::::::::::
; Resources to make only line contours
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::

procedure mycontour(res:logical, opt:logical)


begin
res                 = True                ; make plot mods
        if (opt.and.isatt(opt,"clevs")) then
        res@cnLevelSelectionMode = "ExplicitLevels"
        res@cnLevels             = opt@clevs
        end if


res@lbLabelBarOn    = False
res@cnSmoothingOn   = True
res@cnFillOn              = False
res@cnLinesOn             = True
res@cnLineLabelsOn        = False    ; turn off line labels
res@cnInfoLabelOn         = False
res@cnLineColor = "black"
res@gsnContourNegLineDashPattern = 1 ;

end
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::

;:::::::::::::::::::::::::::::::::::::::::::::::::::::::
; Resources to make only shaded contours
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::

 procedure myfill(res:logical, opt:logical)


begin
res                 = True                ; make plot mods
        if (opt.and.isatt(opt,"clevs")) then
        res@cnLevelSelectionMode = "ExplicitLevels"
        res@cnLevels             = opt@clevs
        end if


res@lbLabelBarOn    = False
res@cnSmoothingOn   = True
res@cnFillOn              = True
res@cnLinesOn             = False
res@cnLineLabelsOn        = False    ; turn off line labels
;res@cnLineColor = "black"
res@cnInfoLabelOn         = False
;res@gsnContourNegLineDashPattern = 1 ;
res@gsnSpreadColors = True

end
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::


;:::::::::::::::::::::::::::::::::::::::::::::::::::::::
; Resources to make vectors plus line contours on top of shaded contours
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::

procedure myvectorcontour(res:logical, opt:logical)

begin
res                 = True                ; make plot mods
        if (opt.and.isatt(opt,"clevs")) then
        res@cnLevelSelectionMode = "ExplicitLevels"
        res@cnLevels             = opt@clevs
        end if

res@lbLabelBarOn    = False
res@cnSmoothingOn   = True
res@cnFillOn              = True
res@cnLinesOn             = True
res@cnLineColor = "black"
res@cnLineLabelsOn        = False    ; turn off line labels
res@cnInfoLabelOn         = False
res@gsnScalarContour     = True
res@gsnContourNegLineDashPattern = 1 ;
res@gsnSpreadColors = True

res@mpMonoFillColor       = True
res@mpFillOn  = False
res@mpDefaultFillColor = "white"
res@mpDefaultFillPattern=-1
res@mpOutlineOn           = True
res@mpOutlineDrawOrder = "PostDraw"
end


procedure myvectorfill(res:logical, opt:logical)

begin
res                 = True                ; make plot mods
        if (opt.and.isatt(opt,"clevs")) then
        res@cnLevelSelectionMode = "ExplicitLevels"
        res@cnLevels             = opt@clevs
        end if

res@lbLabelBarOn    = False
res@cnSmoothingOn   = True
res@cnFillOn              = True
res@cnLinesOn             = False
res@cnLineColor = "black"
res@cnLineLabelsOn        = False    ; turn off line labels
res@cnInfoLabelOn         = False
res@gsnScalarContour     = True
res@gsnSpreadColors = True

res@mpMonoFillColor       = True
res@mpFillOn  = False
res@mpDefaultFillColor = "white"
res@mpDefaultFillPattern=-1
res@mpOutlineOn           = True
res@mpOutlineDrawOrder = "PostDraw"
end

procedure setvcskip(res,skip)
begin
res@vcMinDistanceF = skip ; to thin out the arrows
end


procedure setvcscl(res,mag,len)
begin
res@vcRefMagnitudeF           = mag
res@vcRefLengthF              = len
end

procedure setvc(res,color,style,thick)
;# Styles are "FillArrow", "LineArrow", "WindBarb", "CurlyVector"
begin
res@vcGlyphStyle            = style
res@vcLineArrowColor        = color           ; change vector color
res@vcLineArrowThicknessF   = thick               ; change vector thickness
res@vcVectorDrawOrder       = "PostDraw"  
end

procedure setvcanno(res,color,up)
begin
res@vcRefAnnoOrthogonalPosF   = up
res@vcRefAnnoArrowLineColor   = color
res@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref
end

procedure novcanno(res)
begin
res@vcRefAnnoOn=False
end

procedure setcn(res,clr,thick)
begin
res@cnLineColor=clr
res@cnLineThicknessF=thick
end

procedure setcstyle(res,style)
begin
res@xyDashPatterns    = style
end

procedure setstring(res,fhgt,left,right)
begin
res@gsnStringFontHeightF = fhgt
res@gsnLeftString=left
res@gsnRightString=right
end

function openwks(format:string,filename:string,colormap:string)
begin
wks  = gsn_open_wks(format,filename)              ; open a ps plot
gsn_define_colormap(wks,colormap)              ; choose colormap
return wks
end

function open_wks(format:string,filename:string,colormap:string)
begin
  return(openwks(format,filename,colormap))
end

procedure sety1y2(res,ymin,ymax)
begin
  res@trYMinF  = ymin
  res@trYMaxF  = ymax
end

procedure xt_minmax(res,xmin,xmax,units)
begin
hour=0
minute=0
second=0
yr1=xmin(0)
mo1=xmin(1)
dy1=xmin(2)
yr2=xmax(0)
mo2=xmax(1)
dy2=xmax(2)

x1=ut_inv_calendar(yr1,mo1,dy1,hour,minute,second,units,0)
x2=ut_inv_calendar(yr2,mo2,dy2,hour,minute,second,units,0)
res@trXMinF=x1
res@trXMaxF=x2
end

procedure yt_minmax(res,ymin,ymax,units)
begin
yr1=ymin(0)
mo1=ymin(1)
dy1=ymin(2)
yr2=ymax(0)
mo2=ymax(1)
dy2=ymax(2)

hour=0
minute=0
second=0
y1=ut_inv_calendar(yr1,mo1,dy1,hour,minute,second,units,0)
y2=ut_inv_calendar(yr2,mo2,dy2,hour,minute,second,units,0)
res@trYMinF=y1
res@trYMaxF=y2
end

procedure x_minmax(res,x1,x2)
begin
res@trXMinF=x1
res@trXMaxF=x2
end

procedure y_minmax(res,y1,y2)
begin
res@trYMinF=y1
res@trYMaxF=y2
end

procedure drawmarker(res,mark:string)
;mark can be "Markers" or "Lines" or "MarkLines"
begin
res@xyMarkLineModes=mark
end

procedure yrefline(res,refline,color1,color2)
begin
res@gsnYRefLine=refline
res@gsnAboveYRefLineColor=color1
res@gsnBelowYRefLineColor=color2
end

procedure noframe(res)
begin
  res@gsnFrame=False
end
procedure nodraw(res)
begin
  res@gsnDraw=False
end

procedure DebugPlot(var)
begin
  if isatt(var,"type")
    type=var@type
  else
    type="x11"
  end if
  wks = openwks(type, "test", "testcmap")
  res = True
  opt = True
  if isatt(var,"clevs")
    opt@clevs=var@clevs
  end if
    mycontour(res,opt)
  if isatt(var,"shaded")
    myfill(res,opt)
  end if
  if isatt(var,"raster")
    myfill(res,opt)
    res@cnFillMode = "RasterFill"
  end if
  if isatt(var,"title")
    res@tiMainString=var@title
  end if
  if isatt(var,"colors")
     gsn_define_colormap(wks,var@colors)
     res@gsnSpreadColors = False
  end if
  res@lbLabelBarOn=True
  res@lbLabelAutoStride=True
  if isatt(var,"map")
    lat = var&lat
    lon = var&lon
    if isatt(var,"latlon")
      opt@latlon=var@latlon
    else
      opt@latlon=(/min(lat),max(lat),min(lon),max(lon)/)
    end if
    mymap(res,opt)
    plot = gsn_csm_contour_map(wks,var,res)
  else
    plot = gsn_csm_contour(wks,var,res)
  end if
end

procedure DebugGraph(var)
begin
  wks  = gsn_open_wks("x11","a graph")
  plot = gsn_csm_y(wks,var,False)
end
