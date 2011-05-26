load "$AFS2/lib/ncl/helper_libs.ncl"
load "./conc_plot_helpers.ncl"
load "$WRF/draw_pref.ncl"

begin
  yr=2011
  mo=3
  dy=11
  mon=sprinti("%02i",mo)
  hr="18:00:00"
  dom="d03"
  hora=str_split(hr,":")
  plot_root = getenv("WRF")+"/WRF-FKS-PLOT2/"
  plot_root="./"
  plot_type = "x11"
  vname = "drydep"
  ;vname = "conc"
  ilev = 3
  nt=1
  istart=0
  iend=nt-1
  path="../"+vname+".dat"

  ngx1 = 0;  40
  ngy1 = 0;  70
  ngx2 = 0;  80
  ngy2 = 0;  60

  true_lat = 37.561
  true_lon = 148.43
  lat1 = 17.496
  lon1 = 129.961
  lat2 = 52.861
  lon2 = 178.697


  fin=fopen_read("cord2d."+dom+".nc")
  lat2d=fin->lat2d
  lon2d=fin->lon2d
  dom_dims = dimsizes(lat2d)
  dom_rank = dimsizes(dom_dims)
  nx1 = dom_dims(dom_rank - 1) - 1
  ny1 = dom_dims(dom_rank - 2) - 1

  lats  = lat2d(ngy1,ngx1)
  lons  = lon2d(ngy1,ngx1)
  late  = lat2d(ny1-ngy2,nx1-ngx2)
  lone  = lon2d(ny1-ngy2,nx1-ngx2)


  npp_lat = 37.4230
  npp_lon = 141.0329

  lat1 = lats
  lon1 = lons
  lat2 = late
  lon2 = lone

  conc = read_conc(vname,path,nx1,ny1,ngx1,ngx2,ngy1,ngy2,ilev)
  nt = conc@nt
  Times = asciiread("../time_list",(/nt,4/),"integer")
  print(Times(nt-1,:))

  res                        = True               ; plot mods desired
  set_plot_resources(res,true_lat, true_lon, lat1,lon1,lat2,lon2)

  if vname.eq."drydep"
    res@tiMainString           = "Dry Deposition (normalized by source conc) %"
  else
    res@tiMainString           = "Concentration at "+flx_levs(ilev)+"m (normalized by source conc) %"
  end if
  opt = True
  opt@clevs = (/     0.1, 1, 10, 30, 60, 100, 150, 300, 600/)
  myfill(res,opt)
  res@cnFillColors = (/2, 3,  5, 11, 13, 14, 15,   16, 17, 18/) 
  res@gsnSpreadColors        = False               ; use total colormap
  res@lbLabelBarOn = True

  plot = new(1, graphic)
  ;plot@resolution="prefectures"
  do it=istart,nt-1
    time=Times(it,:)
    plot_dir = plot_root+"init."+yr+"-"+mo+"-"+dy+"-"+hora(0)+"/"+time(0)+"/"+time(1)+"/"
    plot_dir = plot_dir+time(2)+"/"+time(3)+"/"
    ;system("mkdir -p "+plot_dir)
    plot_file = plot_dir+vname
    wks  = gsn_open_wks (plot_type, plot_file)          ; open workstation
    gsn_define_colormap (wks,"precip3_16lev")         ; choose color map

    conc_sfc=conc(it,:,:)
    res@gsnLeftString = str_join(Times(it,:2), "-")+"--"+Times(it,3)+":00hr"
    plot = gsn_csm_contour_map(wks,conc_sfc,res)     ; Draw contours over a map.
    mark_point(wks,plot,npp_lat, npp_lon)
    plot@resolution="prefectures"
    draw_pref_on_plot_v2(wks,plot)
    frame(wks)
    delete(wks)
  end do
end
