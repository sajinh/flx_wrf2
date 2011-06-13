load "$AFS2/lib/ncl/helper_libs.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

function gpm_aave(f)
begin
  lt1=36.5978892330702
  lt2=40.9798979696202
  ln1=135.0000001
  ln2=143.4374999
  var=f->TMP_P0_L103_GLL0({lt1:lt2},{ln1:ln2})
  return(wgt_areaave(var,1,1,0))
end

begin
  wrf_dir="../wrf-demo/"
  fctime=(/"18","24","30"/)
  fnam="10071918.gfs.t18z.pgrb2f"
  var = new(3,float)
  type = "x11"
  wks = gsn_open_wks(type,"ofil")
  gsn_define_colormap(wks,"temp_19lev")
  do itim = 0,2
    f = addfile(wrf_dir+fnam+fctime(itim)+".grb","r")
    var(itim) = gpm_aave(f)
  end do
  var = var - 273.16
  DebugGraph(var)
end

