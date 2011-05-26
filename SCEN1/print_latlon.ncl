load "$AFS2/lib/ncl/helper_libs.ncl"
load "$WRF/WRFUserARW_add.ncl"

IN_TMP=False
begin
  yr=2011
  imo=3
  idy=12
  mo=3
  dy=14
  mon=sprinti("%02i",mo)
  hr="18:00:20"
  root="/home/saji/WRF-FKS/"
  dir="out/"+yr+"/"+imo+"/"+idy
  if IN_TMP
    dir="run"
  end if
  dom="d03"

  fin = addfile(root+dir+"/wrfout_"+dom+"_"+yr+"-"+mon+"-"+dy+"_"+hr+".nc", "r")
  lat1 = wrf_user_getvar(fin,"XLAT",0)
  lon1 = wrf_user_getvar(fin,"XLONG",0)

  dom_dims = dimsizes(lat1)
  dom_rank = dimsizes(dom_dims)
  nx1 = dom_dims(dom_rank - 1) - 1
  ny1 = dom_dims(dom_rank - 2) - 1

  lats  = lat1(0,0)
  lons  = lon1(0,0)
  late  = lat1(ny1,nx1)
  lone  = lon1(ny1,nx1)

  info(lat1)
  pp(lats+" : "+lons)
  pp(late+" : "+lone)

  lats1  = lat1(1,0)
  lons1  = lon1(0,1)
  late1  = lat1(ny1-1,nx1)
  lone1  = lon1(ny1,nx1-1)
  pp(lats1+" : "+lons1)
  pp(late1+" : "+lone1)

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

  fout=fopen_write("cord2d."+dom+".nc")
  xlat2 = lat1(1:(ny1-1),1:(nx1-1))
  xlon2 = lon1(1:(ny1-1),1:(nx1-1))
  fout->lat2d=xlat2
  fout->lon2d=xlon2

end
