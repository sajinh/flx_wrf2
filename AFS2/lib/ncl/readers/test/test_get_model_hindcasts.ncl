load "$AFS2/lib/ncl/helper_libs.ncl"
load "$AFS2/lib/ncl/AFS.Definitions"
load "$AFS2/lib/ncl/readers/mme_in_helper.ncl"

models=(/"COLA","NCEP","POAMA"/)
var_name="prec"

opt=True
db_root="/apcc11/AFS_DB"
opt@season="JAS"

opt@season="NONE"
lead_time=3
start_mon=7
;opt@months=8
hindcast_years=find_hindcast_years(db_root,models,var_name,start_mon,lead_time,opt)
start_years = ispan(hindcast_years(0),hindcast_years(1),1)
var1=multi_year_hindcasts_as_anom(db_root,models,\
                     var_name,start_mon,lead_time,hindcast_years,opt)

info(var1)
var1@shaded=True
var1@map=True
DebugPlot(var1(0,0,:,:,0))

