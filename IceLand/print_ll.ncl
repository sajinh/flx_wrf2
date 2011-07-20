load "$AFS2/lib/ncl/helper_libs.ncl"
load "$WRF/WRFUserARW_add.ncl"

IN_TMP=True
begin
  yr=2011
  imo=5
  idy=22
  mo=5
  dy=22
  mon=sprinti("%02i",mo)
  hr="00:00:00"
  root="/home/saji/IceLand/"
  dir="out/"+yr+"/"+imo+"/"+idy
  if IN_TMP
    dir="./"
  end if
  dom="d01"

  fin = addfile(root+dir+"/wrfout_"+dom+"_"+yr+"-"+mon+"-"+dy+"_"+hr+".nc", "r")
  lat1 = wrf_user_getvar(fin,"XLAT",0)
  lon1 = wrf_user_getvar(fin,"XLONG",0)

  dom_dims = dimsizes(lat1)
  dom_rank = dimsizes(dom_dims)
  nx1 = dom_dims(dom_rank - 1) - 1
  ny1 = dom_dims(dom_rank - 2) - 1

  locs = wrf_user_ij_to_ll(fin, 0, 0, True)
  loce = wrf_user_ij_to_ll(fin, nx1, ny1, True)

  lats  = locs(1)
  lons  = locs(0)
  late  = loce(1)
  lone  = loce(0)

  pp("Original Bounds")
  pp("start "+lats+" : "+lons)
  pp("end "+late+" : "+lone)
  pp("      .        ")

  locs = wrf_user_ij_to_ll(fin, 1, 1, True)
  loce = wrf_user_ij_to_ll(fin, nx1-1, ny1-1, True)

  lats1  = locs(1)
  lons1  = locs(0)
  late1  = loce(1)
  lone1  = loce(0)
  pp("Bounds moved 1 grid inwards")
  pp("start "+lats1+" : "+lons1)
  pp("end "+late1+" : "+lone1)
  pp("      .        ")

  pp("grid resolution")
  dx_ll=lons1-lons
  dx_ur=lone-lone1

  dy_ll=lats1-lats
  dy_ur=late-late1

  pp("x resolution : ["+dx_ll+" ,"+dx_ur)
  pp("y resolution : ["+dy_ll+" ,"+dy_ur)

  pp("New Corners are:")
  pp((lats+dy_ll)+" , "+(lons+dx_ll))
  pp((late-dy_ur)+" , "+(lone-dx_ur))
  ;pp((late-dy_ll)+" , "+(lone-dx_ll))

  fout=fopen_write("cord2d."+dom+".nc")
  xlat2 = lat1(1:(ny1-1),1:(nx1-1))
  xlon2 = lon1(1:(ny1-1),1:(nx1-1))
  fout->lat2d=xlat2
  fout->lon2d=xlon2

end
