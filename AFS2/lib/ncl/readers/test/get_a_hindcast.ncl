load "$AFS2/lib/ncl/helper_libs.ncl"
load "$AFS2/lib/ncl/AFS.Definitions"
load "$AFS2/lib/ncl/readers/mme_in_helper.ncl"

models="COLA"
var_name="prec"
opt=True
db_root="/apcc11/AFS_DB"
lead_time=3
start_mon=8
hindcast_year=2002
var=multi_year_hindcasts_as_anom(db_root,models,\
                     var_name,start_mon,lead_time,hindcast_year,opt)

info(var)
delete(var)

opt@season="ASO"
var=multi_year_hindcasts_as_anom(db_root,models,\
                     var_name,start_mon,lead_time,hindcast_year,opt)
info(var)
delete(var)
opt@season="NONE"
opt@months="AUG"
var=multi_year_hindcasts_as_anom(db_root,models,\
                     var_name,start_mon,lead_time,hindcast_year,opt)
info(var)
delete(var)

hindcast_years=ispan(1982,2002,1)
var1=multi_year_hindcasts_as_anom(db_root,models,\
                     var_name,start_mon,lead_time,hindcast_years,opt)

var1@shaded=True
var1@map=True
var1@title="Simple Plot w/ friendly AFS libraries"
var1@type="eps"
DebugPlot(var1(1,0,:,:,0))

