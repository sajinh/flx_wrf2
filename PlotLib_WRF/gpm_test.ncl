load "$AFS2/lib/ncl/helper_libs.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

function gpm_aave(f)
begin
  lt1=36.5978892330702
  lt2=40.9798979696202
  ln1=135.0000001
  ln2=143.4374999
 lt1=35.
  lt2=39.0
  ln1=-60
  ln2=-100
  var=f->TMP_212_SFC({lt1:lt2},{ln1:ln2})
  return(wgt_areaave(var,1,1,0))
end

begin
  ;wrf_dir="../test_get_data/"
  wrf_dir="../wrf-demo/JAN00/"
  fctime=(/"2000012412","2000012418","2000012500"/)
  fnam=".AWIPSF.grb"
  var = new(3,float)
  type = "x11"
  wks = gsn_open_wks(type,"ofil")
  gsn_define_colormap(wks,"temp_19lev")
  do itim = 0,2
    f = addfile(wrf_dir+fctime(itim)+fnam,"r")
    var(itim) = gpm_aave(f)
  end do
  var = var - 273.16
  DebugGraph(var)
end

