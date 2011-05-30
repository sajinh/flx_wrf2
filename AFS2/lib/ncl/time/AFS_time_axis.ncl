function make_month_points(opt)
local opt
begin
  inc=opt@inc
  Date1=opt@_t1
  Date2=opt@_t2
  yr1=floattointeger(Date1(0,0))
  mo1=floattointeger(Date1(0,1))
  dy1=floattointeger(Date1(0,2))
  tunits=opt@_tunits
  nyrs=opt@_nyrs
  nmos=opt@_nmos

  newpoints=new(nmos,"double")
  newmonths=new(nmos,"integer")
  yr=yr1
  mo=mo1
  dy=15 ; center a  month on the 15th day
  do i = 0,nmos-1
    newpoints(i)=ut_inv_calendar(yr,mo,dy,0,0,0,tunits,0)
    newmonths(i)=mo
    mo=mo+1
    if mo.gt.12
      yr=yr+1
      mo=mo-12
    end if
  end do
newpoints@months=newmonths(::inc)
return(newpoints(::inc))
end

function make_month_labels(points,opt)
begin
 nmos=dimsizes(points)
 labels=new(nmos,"string")
 months=points@months

 if opt.eq.0
   labels=points@months
 else
   do i = 0,nmos-1
   labels(i)=integer_to_month(months(i))
   end do
 end if
 return(labels)
end


function make_year_points(opt)
begin
  inc=opt@inc
  yr_month=1
  Date1=opt@_t1
  Date2=opt@_t2
  yr1=floattointeger(Date1(0,0))
  yr2=floattointeger(Date2(0,0))
  mo1=Date1(0,1)
  mo2=Date2(0,1)
  tunits=opt@_tunits
  nyrs=opt@_nyrs
  if mo1.gt.yr_month
    yr1=yr1+1
    nyrs=nyrs-1
  end if
  if mo2.lt.yr_month
     yr2=yr2-1
     nyrs=nyrs-1
  end if
  newpoints=new(nyrs,"double")
  yr=yr1
  mo=yr_month
  dy=15 ; center a  month on the 15th day
  do i = 0,nyrs-1
    newpoints(i)=ut_inv_calendar(yr,mo,dy,0,0,0,tunits,0)
    yr=yr+1
  end do
  newpoints@years=ispan(yr1,yr2,inc)
  return(newpoints(::inc))
end

function make_year_labels(points)
begin
  labels=points@label
  nlabs=dimsizes(labels)
  new_labels=new(nlabs,"string")
  do i = 0,nlabs-1
    new_labels(i)=labels(i)
  end do
  return(new_labels)
end

procedure set_time_on_xaxis(res,opt)
begin
  if opt@_nyrs.gt.1
    res@tmXBMode         = "Explicit"
    res@tmXBValues       = opt@yr_values
    res@tmXBLabels       = opt@yr_labels
    res@tmXBLabelDeltaF  = 0.6     ; Move label away from tickmarks.
  else
    res@tmXBMode         = "Explicit"
    res@tmXBValues       = opt@mo_values
    res@tmXBLabels       = opt@mo_labels
    res@tmXBLabelDeltaF  = 0.6     ; Move label away from tickmarks.
  end if
  nodraw(res)
end

procedure set_time_on_yaxis(res,opt)
begin
  if opt@_nyrs.gt.1
  res@tmYLMode         = "Explicit"
  res@tmYLValues       = opt@yr_values
  res@tmYLLabels       = opt@yr_labels
  res@tmYLLabelAngleF  = 90
  res@tmYLLabelJust    = "CenterCenter"
  res@tmYLLabelDeltaF  = 0.9     ; Move label away from tickmarks.
  end if
  nodraw(res)
end

procedure draw_a_blank_plot(wks,plot,opt)
begin
    ; plot months as minor ticks
  if isatt(opt,"MonthInc")
    opt@inc=opt@MonthInc
  else
    opt@inc=3
  end if
  if opt@_nyrs.lt.1
    opt@inc=1
  end if
  mo_values=make_month_points(opt)
  mo_values@label="a"
  if isatt(opt,"MonthFmt")
    opt@fmt=opt@MonthFmt
  else
    opt@fmt=0
  end if
  mo_labels=make_month_labels(mo_values,opt@fmt)

  if opt@NoMinor.eq.False 
    getvalues plot 
      "vpXF"      : vpxf
      "vpYF"      : vpyf
      "vpHeightF" : vpheightf
      "vpWidthF"  : vpwidthf
      "tmXBLabelFontHeightF" : tickheightf
    end getvalues

    res2           = True
    res2@vpXF      = vpxf
    res2@vpYF      = vpyf
    res2@vpHeightF = vpheightf
    res2@vpWidthF  = vpwidthf

    if opt@_timeOrdinate
      res2@tmYLMode                = "Explicit"
      res2@tmYLValues              = mo_values 
      res2@tmYLLabels              = mo_labels 
      res2@tmYLLabelFontHeightF    = tickheightf/1.5 
      res2@tmYLMajorOutwardLengthF = 0.015     ; Don't draw tickmarks b/c they
      res2@tmYLMajorLengthF        = 0.008     ; were drawn on previous plot.
      res2@tmYLLabelDeltaF         = -0.2     ; Move label away from tickmarks.
      res2@tmYROn = False            ; Turn off right tickmarks.
      res2@tmXBOn = False            ; Turn off top tickmarks.
      res2@tmXTOn = False            ; Turn off left tickmarks.
    else
      res2@tmXBMode                = "Explicit"
      res2@tmXBValues              = mo_values 
      res2@tmXBLabels              = mo_labels 
      res2@tmXBLabelFontHeightF    = tickheightf/1.5 
      res2@tmXBMajorOutwardLengthF = 0.015     ; Don't draw tickmarks b/c they
      res2@tmXBMajorLengthF        = 0.008     ; were drawn on previous plot.
      res2@tmXBLabelDeltaF         = -0.2     ; Move label away from tickmarks.
      res2@tmYROn = False            ; Turn off right tickmarks.
      res2@tmXTOn = False            ; Turn off top tickmarks.
      res2@tmYLOn = False            ; Turn off left tickmarks.
    end if

    blank = gsn_blank_plot(wks,res2)    ; Create a blank plot.
    overlay(plot,blank)                 ; Overlay on existing plot.
  end if
;  end if
end

function monthly_timeseries(wks,time,data,res,opt)
begin
  ;   if nyr .gt. 1, plot the graph with years as major tickmarks
  ;   add months as minor tickmarks
  ;if opt@_nyrs.gt.1
    ; plot years as major ticks
    plot=gsn_csm_xy(wks,time,data,res)
  ;end if
  return(plot)
end

function monthly_contour(wks,data,res,opt)
begin
  plot=gsn_csm_contour(wks,data,res)
return(plot)
end
procedure time_plot_controller(time,_res,_opt)
local plot
begin
  print("Time series plot")
  ; retrieve the resources into a separate resource variable, __res perhaps
  ; check if x-axis data is valid time object
  ; find data boundaries from time object or from resource settings
  ;                  trXMaxF and trXMinF
  ; depending on data resolution or space constraints, plot
  ; yearly, monthly or daily tickmarks
  ;    e.g: if yearly, plot years as major tickmarks at app. intervals
  ;    e.g: if monthly, plot years as major ticks, months as minor
  ;    e.g: if daily, plot months as major ticks, days as minor : optionally
  ;           the year can be indicated in one of the month (major) ticks
  
  ; return plot data and pass it onto caller sub-program

  Date1=get_yr_mo_dy(time(0))
  yr1=Date1@year
  mo1=Date1@month
  dy1=Date1@day

  ntim=dimsizes(time)
  Date2=get_yr_mo_dy(time(ntim-1))
  yr2=Date2@year
  mo2=Date2@month
  dy2=Date2@day

  xm1=0
  xm2=0
  if mo1.gt.1
    xm1=-1*(mo1-1) ; reduce this from total no of months
  end if
  if mo2.lt.12
    xm2=-1*(12-mo2) ; reduce this from total no of months
  end if
  nyrs=yr2-yr1+1
  nmos=nyrs*12+xm1+xm2

 ; _res=res
 ; _opt=True
  _opt@_t1=Date1
  _opt@_t2=Date2
  _opt@_nyrs=nyrs
  _opt@_nmos=nmos
  _opt@_tunits=time@units
   if isatt(_res,"NoMinor")
     _opt@NoMinor=True
     delete(_res@NoMinor)
   else
     _opt@NoMinor=False
   end if
   if isatt(_res,"YearInc")
     _opt@YearInc=_res@YearInc
     delete(_res@YearInc)
   end if
   if isatt(_res,"MonthInc")
     _opt@MonthInc=_res@MonthInc
     delete(_res@MonthInc)
   end if
  nyrs=_opt@_nyrs
  nmos=_opt@_nmos

  if isatt(_opt,"YearInc")
    _opt@inc=_opt@YearInc
  else
    _opt@inc=1
  end if

  opt = _opt
  if nyrs.gt.1
    yr_values=make_year_points(_opt)
    _opt@yr_values=yr_values
    _opt@yr_labels=yr_values@years
  else
    mo_values=make_month_points(opt)
    mo_values@label="a"
    if isatt(opt,"MonthFmt")
      opt@fmt=opt@MonthFmt
    else
      opt@fmt=0
    end if
    mo_labels=make_month_labels(mo_values,opt@fmt)
    _opt@mo_values=mo_values
    _opt@mo_labels=mo_labels
  end if
end
procedure draw_the_plot(res,plot,wks)
begin
  if isatt(res,"gsnDraw") 
    if res@gsnDraw
      draw(plot)
    end if
  else
    draw(plot)
  end if
  if isatt(res,"gsnFrame") 
    if res@gsnFrame
      frame(wks)
    end if
  else
    frame(wks)
  end if
end

function time_series_plot(wks,time,data,res)
; options
; _opt@NoMinor=False

begin
  _opt=True
  _res=res
  _opt@_timeOrdinate=False
  time_plot_controller(time,_res,_opt)
  set_time_on_xaxis(_res,_opt)
  plot=monthly_timeseries(wks,time,data,_res,_opt)
  draw_a_blank_plot(wks,plot,_opt)
  draw_the_plot(res,plot,wks)
  return(plot)
end

function contour_with_time_ordinate(wks,time,data,res)
begin
  _opt=True
  _res=res
  time_plot_controller(time,_res,_opt)
  set_time_on_yaxis(_res,_opt)
  plot=monthly_contour(wks,data,_res,_opt)
  _opt@_timeOrdinate=True
  draw_a_blank_plot(wks,plot,_opt)
  draw_the_plot(res,plot,wks)
;  draw(plot)
;  frame(wks)
  return(plot)
end
