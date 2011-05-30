load "$AFS2/lib/ncl/helper_libs.ncl"
load "$AFS2/lib/ncl/AFS.Definitions"
load "$AFS2/lib/ncl/readers/mme_in_helper.ncl"

model="COLA"
var_name="prec"

opt=True
db_root="/apcc01/OPER/AFS/DATA"
data_type="Hindcast"
data_period=(/6,8,1990,1994/)
opt@ensemble_average=False
opt@season="JJA"

var1=mme_in_reader(db_root,data_type,model,var_name,data_period,opt)
print(var1&year)
info(var1)

exit

opt@nlev=3
;var2=mme_in_reader(db_root,data_type,model,var_name,data_period,opt)
;info(var2)

opt@ensemble_average=True
;var3=mme_in_reader(db_root,data_type,model,var_name,data_period,opt)
;info(var3)

data_type="Forecast"
data_period=(/6,8,2006,2007/)
opt@season="NONE"
;var4=mme_in_reader(db_root,data_type,model,var_name,data_period,opt)
;info(var4)


opt@months=(/6,7/)
print(opt)
var5=mme_in_reader(db_root,data_type,model,var_name,data_period,opt)
info(var5)

exit

data_period=(/6,8,1901,1901/)
opt@terminate_if_no_data=True
var105=mme_in_reader(db_root,data_type,model,var_name,data_period,opt)
info(var105)

exit

opt@terminate_if_no_data=False
var106=mme_in_reader(db_root,data_type,model,var_name,data_period,opt)
info(var106)
