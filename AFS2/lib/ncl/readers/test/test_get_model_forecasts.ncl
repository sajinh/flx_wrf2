load "$AFS2/lib/ncl/helper_libs.ncl"
load "$AFS2/lib/ncl/AFS.Definitions"
load "$AFS2/lib/ncl/readers/mme_in_helper.ncl"

models=(/"SINT","NCEP","POAMA"/)
var_name="prec"
opt=True
db_root="/apcc01/OPER/AFS/DATA"
data_type="Forecast"
opt@ensemble_average=True
opt@season="JFM"
opt@season="NONE"
;opt@months=(/"Feb","Mar"/)
lead_time=6
start_mon=12
start_year=2008
var1=single_year_predictions(db_root,data_type,models,var_name, \
                              start_mon,lead_time,start_year,opt)
info(var1)
exit
;var1@map=True
;var1@shaded=True
;DebugPlot(var1(0,0,:,:,0))

var2=single_year_predictions_as_anom(db_root,data_type,models,var_name, \
                              start_mon,lead_time,start_year,opt)
info(var2)
;var2@map=True
;var2@shaded=True
;DebugPlot(var2(0,0,:,:,0))

