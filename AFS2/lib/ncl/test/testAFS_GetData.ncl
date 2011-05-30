load "$AFS_NCL/AFS.Definitions"
load "$AFS_NCL/AFS_GetData.ncl"


start_mon=5 ; December
lead_time=6  ; 3-month lead
end_mon=start_mon+lead_time-1

;if end_mon.gt.12
;  end_mon=end_mon-12
;end if

var_name="t850" 

opt=True

opt@DataPeriod=(/start_mon,end_mon,2006,2007/)
opt@Lat=(/-5,5/)
opt@Lon=(/210,280/)

; Get the Data from Data base using provided API
var=Get_OBS_Data(var_name, opt)

printVarSummary(var)
  
