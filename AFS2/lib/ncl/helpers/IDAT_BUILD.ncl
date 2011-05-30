; =============================================================================
;
;  File Name: IDAT_BUILD.ncl
;
; -----------------------------------------------------------------------------
;
;   Author: Soo-Jin Sohn
;     APCC(APEC Climate Center), Busan, Republic of Korea
; 
;   Date: Fri Mar 2 2007
; 
;   Description:
;     This script defines all of the basic functions and procedures used in the
;     APCC seasonal forecast and hindcast based on MRG and SSE methods.
; 
;     To use the functions and procedures in this scripts,
;     you must have the line:
;   
;     load "IDAT_BUILD.ncl"
;
;     at the top of your NCL scripts, before the begin statement.
;    
;   Code Description:
;     NCAR Command Language (NCL) Version  4.2.0.a035
;
; -----------------------------------------------------------------------------


; -----------------------------------------------------------------------------
;  For every function and procedure definded in this scripts, undefine it
;  with a call to "undef" so it doesn't clash with other functions and 
;  procedures with this same name.
; -----------------------------------------------------------------------------

undef ("seasonal_mean")
undef ("fill_up_missing_val")
undef ("shuffle")
undef ("bias_correction")


; -----------------------------------------------------------------------------
;  procedure : seasonal_mean
;                ARR[*][*][*][*][*]:float
;                OBS[*][*][*][*]:float
;                xARR[*][*][*][*][*]:float
;                xOBS[*][*][*][*]:float
;
;  Calculates seasonal means of model and observation datasets and add them to 
;  monthly mean data.
; -----------------------------------------------------------------------------

procedure seasonal_mean(ARR[*][*][*][*][*]:float, OBS[*][*][*][*]:float , \
                        xARR[*][*][*][*][*]:float, xOBS[*][*][*][*]:float)
local sARR, sOBS, dims, nmon
begin

  sARR = dim_avg_Wrap(ARR(model|:,time|:,latitude|:,longitude|:,month|:))
  sOBS = dim_avg_Wrap(OBS(time|:,latitude|:,longitude|:,month|:))

  dims = dimsizes(xARR)
  nmon = dims(2)

  xARR(:,:,0:nmon-2,:,:) = ARR(:,:,:,:,:)
  xARR(:,:,nmon-1,:,:) = sARR(:,:,:,:)
  xOBS(:,0:nmon-2,:,:) = OBS(:,:,:,:)
  xOBS(:,nmon-1,:,:) = sOBS(:,:,:)

end


; -----------------------------------------------------------------------------
;  function : fill_up_missing_val
;                OBS[*][*][*][*]:float
;                ARR[*][*][*][*][*]:float
;         
;               
;  Calculates observed climatology except target year and fill up missing 
;  values in model datasets using it. 
; -----------------------------------------------------------------------------

function fill_up_missing_val(OBS[*][*][*][*]:float, ARR[*][*][*][*][*]:float)
local ARR, OBS, dims, ntim, AVGo, cClim
begin

  dims = dimsizes(ARR)
  ntim = dims(1)

  AVGo = dim_avg_Wrap(OBS(month|:,lat|:,lon|:,year|:ntim-2))
  cClim = conform( ARR, AVGo, (/2,3,4/) )
  ARR=where(ismissing(ARR),cClim,ARR)

  return(ARR)

end


; -----------------------------------------------------------------------------
;  procedure : shuffle
;                ARR[*][*][*][*][*]:float
;                OBS[*][*][*][*]:float
;                iyr:integer
;
;  Suffles model and observation datasets for cross-validation.
; -----------------------------------------------------------------------------

procedure shuffle(ARR[*][*][*][*][*]:float, OBS[*][*][*][*]:float, iyr, \
	  cvARR[*][*][*][*][*]:float, cvOBS[*][*][*][*]:float)
local dims, arr_temp, obs_temp
begin

  dims = dimsizes(ARR)
  
  arr_temp = new((/dims(0), dims(1), dims(2), dims(3), dims(4)/), \
             typeof(ARR), getFillValue(ARR))
  obs_temp = new((/dims(1), dims(2), dims(3), dims(4)/), typeof(ARR), \
             getFillValue(ARR))

  if (iyr .lt. dims(1)-1) then
     arr_temp(:,dims(1)-1-iyr:dims(1)-1,:,:,:) = ARR(:,0:iyr,:,:,:)
     arr_temp(:,0:dims(1)-1-iyr-1,:,:,:)    = ARR(:,iyr+1:dims(1)-1,:,:,:)
     obs_temp(dims(1)-1-iyr:dims(1)-1,:,:,:)   = OBS(0:iyr,:,:,:)
     obs_temp(0:dims(1)-1-iyr-1,:,:,:)      = OBS(iyr+1:dims(1)-1,:,:,:)
  else
     arr_temp  = ARR
     obs_temp  = OBS
  end if

  cvARR  = arr_temp
  cvOBS  = obs_temp

end


; -----------------------------------------------------------------------------
;  function : bias_correction
;                ARR[*][*][*][*][*]:float
;                OBS[*][*][*][*]:float
;
;
;  Calculates model and observed climatology except target year and corrects   
;  the bias of model climatolgy using observed climatology. 
;  2. Calculates calibrated output SCM
; -----------------------------------------------------------------------------

function bias_correction(ARR[*][*][*][*][*]:float, OBS[*][*][*][*]:float)
local dims, AVGf, AVGo, COMf, COMo, new_ARR
begin

  dims = dimsizes(ARR)

  AVGf = \
     dim_avg_Wrap(ARR(model|:,month|:,lat|:,lon|:,year|:dims(1)-2))
  AVGo = dim_avg_Wrap(OBS(month|:,lat|:,lon|:,year|:dims(1)-2))

  COMf = conform(ARR, AVGf, (/0,2,3,4/))
  COMo = conform(ARR, AVGo, (/2,3,4/))

  new_ARR = ARR - COMf + COMo

  copy_VarMeta(ARR, new_ARR)

  return(new_ARR)

end
