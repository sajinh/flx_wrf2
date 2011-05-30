load "$AFS2/lib/ncl/helper_libs.ncl"
load "$AFS2/lib/ncl/AFS.Definitions"
load "$AFS2/lib/ncl/readers/obs_reader.ncl"


opt=True
opt@season="ASO"
data_dir="/home/saji/AFS_OBS"
var_name="prec"
years=ispan(1982,2005,1)
start_mon=8
nmon=1
;opt@months=(/"aug","sep","oct","Nov"/)
var2=get_obs_data(data_dir,var_name,start_mon,nmon,years,opt)
info(var2)
exit
var2@map=True
var2@shaded=True
DebugPlot(var2(0,0,:,:))
var3=get_obs_anom(data_dir,var_name,start_mon,nmon,years,opt)
info(var3)
var3@map=True
var3@shaded=True
DebugPlot(var3(0,0,:,:))
