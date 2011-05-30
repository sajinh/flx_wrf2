load "$AFS2/lib/ncl/helper_libs.ncl"
load "$AFS2/lib/ncl/AFS.Definitions"

models=(/"SINT","NCEP","POAMA"/)
var_name="prec"

opt=True
db_root="/apcc01/OPER/AFS/DATA"
data_type="Forecast"
opt@ensemble_average=True
opt@season="JFM"

lead_time=6
start_mon=12
start_year=2008
years=find_hindcast_years(db_root,models,var_name,start_mon,lead_time,opt)
print(years)

