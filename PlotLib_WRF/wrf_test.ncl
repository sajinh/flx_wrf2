load "$AFS2/lib/ncl/helper_libs.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

function wrf_aave(f,vnam,itim)
begin
  lt1=36.5978892330702
  lt2=40.9798979696202
  ln1=135.0000001
  ln2=143.4374999
  var=wrf_user_getvar(f,vnam,itim)
  lat2d  = f->XLAT(0,:,:)     
  lon2d  = f->XLONG(0,:,:) 
  lat = fspan(lt1,lt2,20)
  lon = fspan(ln1,ln2,30)
  var_r  = rcm2rgrid(lat2d,lon2d,var,lat,lon,0)
  return(wgt_areaave(var_r,1,1,0))
end

begin
  ;wrf_dir="../wrf-oper/test_out/2010/7/20/"
  wrf_dir="../wrf-demo/run/"
  do idom = 1,1
    f = addfile(wrf_dir+"wrfout_d0"+idom+"_2010-07-20_12:00:00.nc","r")
    ofil="wrf_dom_"+idom+".temp"
    type = "x11"
    wks = gsn_open_wks(type,ofil)
    gsn_define_colormap(wks,"temp_19lev")
    times  = wrf_user_list_times(f)  ; get times in the file
    ntimes = dimsizes(times)         ; number of times in the file
    var = new(ntimes,float)
    do i = 0,ntimes-1
      var(i) = wrf_aave(f,"T2",i)
    end do
    var = var - 273.16

    DebugGraph(var)
    delete(var)
    delete(times)
exit
  end do
end

