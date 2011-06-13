procedure draw_pref(wks,plot,latlon)
begin
  lat=str_split(latlon(0),",")
  lon=str_split(latlon(1),",")
  flat=stringtofloat(lat)
  flon=stringtofloat(lon)
  plres = True
 ; plres@gsEdgesOn   = True       ; draw border around polygons
 ; plres@gsEdgeColor = "black"
  gsn_polyline(wks, plot, flon, flat, plres)
end

  ; a function that returns lat and lon coords of each feature
  ; in the shapefile
  function get_feature_latlon(fshp,feat_idx)
  begin
    shp_lon = fshp->x
    shp_lat = fshp->y
    segments = fshp->segments
    geometry = fshp->geometry
    segsDims = dimsizes(segments)
    geomDims = dimsizes(geometry)
    ;
    ; Read global attributes
    ;
    geom_segIndex = fshp@geom_segIndex
    geom_numSegs  = fshp@geom_numSegs
    segs_xyzIndex = fshp@segs_xyzIndex
    segs_numPnts  = fshp@segs_numPnts
    numFeatures = geomDims(0)
    startSegment = geometry(feat_idx, geom_segIndex)
    numSegments  = geometry(feat_idx, geom_numSegs)
    lon_str=""
    lat_str=""
    do seg=startSegment, startSegment+numSegments-1
      startPT = segments(seg, segs_xyzIndex)
      endPT = startPT + segments(seg, segs_numPnts) - 1
      lon_str=lon_str+str_join(shp_lon(startPT:endPT),",")
      lat_str=lat_str+str_join(shp_lat(startPT:endPT),",")
    end do
    latlon=(/lat_str,lon_str/)
    return(latlon)
  end

procedure draw_pref_on_plot(wks,plot)
begin
  
  shp1 = "/home/saji/AizuWeather/ShapeFiles/CHGIS/jp_grid_ken_pgn.shp"
  shp2 = "/home/saji/AizuWeather/ShapeFiles/JGS/fks_adm_area.shp"
  if str_upper(plot@resolution) .eq. "PREFECTURES"
    fshp = fopen_read(shp1)
  else
    fshp = fopen_read(shp2)
  end if
  geometry = fshp->geometry
  geomDims = dimsizes(geometry)
  nfeat=geomDims(0)
  do ifeat=0,nfeat-1
    latlon=get_feature_latlon(fshp,ifeat)
    draw_pref(wks,plot,latlon)
  end do
end

procedure draw_pref_on_plot_v2(wks,plot)
begin
  shp1 = "/home/saji/AizuWeather/ShapeFiles/CHGIS/jp_grid_ken_pgn.shp"
  shp2 = "/home/saji/AizuWeather/ShapeFiles/JGS/fks_adm_area.shp"
  if isatt(plot,"resolution") .and. \
      str_upper(plot@resolution) .eq. "PREFECTURES"
    fshp = fopen_read(shp1)
  else
    fshp = fopen_read(shp2)
  end if
  shp_lon = fshp->x
  shp_lat = fshp->y
  segments = fshp->segments
  geometry = fshp->geometry
  geomDims = dimsizes(geometry)
  segsDims = dimsizes(segments)
  geom_segIndex = fshp@geom_segIndex
  geom_numSegs  = fshp@geom_numSegs
  segs_xyzIndex = fshp@segs_xyzIndex
  segs_numPnts  = fshp@segs_numPnts
  nfeat=geomDims(0)
  plres = True
 ; plres@gsEdgesOn   = True       ; draw border around polygons
 ; plres@gsEdgeColor = "black"
  do ifeat=0,nfeat-1
    startSegment = geometry(ifeat, geom_segIndex)
    numSegments  = geometry(ifeat, geom_numSegs)
    do seg=startSegment, startSegment+numSegments-1
       startPT = segments(seg, segs_xyzIndex)
       endPT = startPT + segments(seg, segs_numPnts) - 1
       gsn_polyline(wks, plot, shp_lon(startPT:endPT), \
                               shp_lat(startPT:endPT), plres)
    end do
  end do
end
