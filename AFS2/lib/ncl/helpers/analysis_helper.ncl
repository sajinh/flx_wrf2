procedure info(var)
begin
  printVarSummary(var)
end

function cr()
begin
  return(inttochar(10))
end

function dectobin(dec_num:integer)
; Mary's response to question on
; NCL list on how to convert a
; decimal digit to binary
local dq, cmd
begin
;---Double quote character
  dq = str_get_dq()

;---UNIX command that converts decimal to binary
  cmd = "echo " + dq + "ibase=10;obase=2;" + dq + dec_num + " | bc "
  print(cmd)

;---Execute UNIX command and return an integer
  return(toint(systemfunc(cmd)))
end


function rank_of(var)
local dims
begin
  dims=dimsizes(var)
  return(dimsizes(dims))
end

function chop(astring)
begin
  chr_str=stringtochar(astring)
  nchr=dimsizes(chr_str)
  if nchr.le.2
    return(astring)
  end if
  return(chartostring(chr_str(:nchr-3)))
end

function boolean_return_value(rvalue)
begin
  if (typeof(rvalue) .eq. "logical")
    return(True)
  end if
  return(False)
end
function not_boolean_return_value(rvalue)
begin
  return(.not. boolean_return_value(rvalue))
end
function returns_False(rvalue)
begin
  if boolean_return_value(rvalue)
     if rvalue .eq. False
       return(True)
     end if
   end if
   return(False)
end


function _seasons()
begin
return( (/"JFM","FMA","MAM","AMJ","MJJ","JJA",  \
          "JAS","ASO","SON","OND","NDJ","DJF"/) )
end
function _months()
begin
return( (/"JAN","FEB","MAR","APR","MAY","JUN", \
         "JUL", "AUG","SEP","OCT","NOV","DEC"/) )
end
  
function start_and_end_months_of_season(org_seas)                               
begin                                                                           
  seas=str_upper(org_seas)                                                      
  c_seas=stringtochar(seas)                                                     
  seas=c_seas(0:2)                                                              
  nmon=dimsizes(c_seas) - 1                                                     
  iseas=ind(_seasons.eq.seas)                                                   
  seas@iseas=iseas                                                              
  seas@start_month=iseas                                                        
  seas@end_month=iseas+nmon-1                                                   
return(seas)                                                                    
end               
function months_for_season(season)
begin
  seas=start_and_end_months_of_season(season)
  return(ispan(seas@start_month,seas@end_month,1))
end
procedure make_the_dir(dirName)
local dirName
begin
  system("mkdir -p "+dirName)
end
procedure pp(aString)
begin
  print(""+aString)
end
procedure TerminateProgram()
begin
  print("Program Terminated")
  exit
end
procedure terminate_program(aString)
begin
  print("-----------%< -------")
  pp(aString)
  print("Program Terminated")
  print("-----------%< -------")
  exit
end
procedure abort(aString)
begin
  terminate_program(aString)
end
procedure add_dimensions(var,att_array)
begin
 ndims=dimsizes(att_array)
  do i = 0,ndims-1
   var!i=  att_array(i)
  end do
end

function to_integer(value)
begin
  dims=dimsizes(value)
  dvalue=new(dims,"double")
  dvalue=value
  int_value=doubletointeger(dvalue)
  return(int_value)
end

function to_float(value)
begin
  if typeof(value) .eq. "string"
    return(stringtofloat(value))
  end if
  if typeof(value) .eq. "character"
    return(charactertofloat(value))
  end if
  if typeof(value) .eq. "double"
    return(doubletofloat(value))
  end if
  return(value)
end
procedure start_end_time(var,stim,etim)
begin

ndims=dimsizes(var(:,0)) ; how many time arrays
mnt=new(ndims,typeof(var))  ; array to hold min time in each array
mxt=mnt                  ; array to hold max time in each array

	do itims = 0,ndims-1
	mnt(itims)=min(var(itims,:))
	mxt(itims)=max(var(itims,:))
	end do

stim=max(mnt)
etim=min(mxt)

; We find the minimum time index in each data set and store as mnt
; The maximum of mnt is the starting time
; Then we find the maximum time index in each data set and store as mxt
; The minimum of mxt is the ending time
; Thus all data considered will include the starting and ending time
end

function m2s(var)

local var, rm_var

begin
ndims=dimsizes(var)
nmon=ndims(0)
var!0="time"
var!1="lat"
var!2="lon"
        rm_var=runave_Wrap(var(lat|:,lon|:,time|:),3,0)
        rm_var(:,:,0)=dim_avg_Wrap(var(lat|:,lon|:,time|0:1))
        rm_var(:,:,nmon-1)=dim_avg_Wrap(var(lat|:,lon|:,time|nmon-2:nmon-1))
return (rm_var)
end

function m2s1d(var)

local var, rm_var

begin
ndims=dimsizes(var)
nmon=ndims(0)
var!0="time"
        rm_var=runave_Wrap(var,3,0)
        rm_var(0)=dim_avg_Wrap(var(0:1))
        rm_var(nmon-1)=dim_avg_Wrap(var(nmon-2:nmon-1))
return (rm_var)
end


procedure part_regress(y,x1,x2,a1,a2)

begin
r1y=escorc(x1,y)
r2y=escorc(x2,y)
r12=escorc(x1,x2)

c=(1-r12^2)
a1=r1y-r2y*r12
a1=a1/c

a2=r2y-r1y*r12
a2=a2/c

end
procedure part_corr(y,x1,x2,a1,a2,s1,s2)

begin
r1y=escorc(x1,y)
r2y=escorc(x2,y)
r12=escorc(x1,x2)

s1=r1y
s2=r2y

c=(1-r12^2)
a1=r1y-r2y*r12
a1=a1/sqrt(c)

a2=r2y-r1y*r12
a2=a2/sqrt(c)

end

;;;;;;;;;;;;;;;;;;

procedure partl_corr(y,x1,x2,a1,a2)

begin
r1y=escorc(x1,y)
r2y=escorc(x2,y)
r12=escorc(x1,x2)

c=(1-r12^2)
a1=r1y-r2y*r12
a1=a1/sqrt(c*(r2y^2))

a2=r2y-r1y*r12
a2=a2/sqrt(c*(r1y^2))

end

function isodd(var)
local odd
begin
  y=round(var/2.0,0)
  x=fabs(var-(y+y))
  if (var.eq.1 .or. x.gt.0) then
    odd=True
  else
    odd=False
  end if
  return(odd)
end


procedure start_and_end(yy,mm,tunits,start_time,end_time)
local start_time,end_time,e_month,yy,mm,nday,tunits

begin
	  start_time=ut_inv_calendar(yy,mm,1,0,0,0,tunits,0)
e_month=mm+2
	if (e_month.gt.12) then
	  nday=days_in_month(yy+1,e_month-12)
	  end_time=ut_inv_calendar(yy+1,e_month-12,nday,0,0,0,tunits,0)
	else
	  nday=days_in_month(yy,e_month)
	  end_time=ut_inv_calendar(yy,e_month,nday,0,0,0,tunits,0)
	end if
end

function dim_other_seasons(var,org_seas)
; this will designed as a function that when
; given a variable and a season will diminish the
; amplitude of the variable at times other than the
; specified season (with undefined values)
; a cosine bell will be fitted

begin

; We should parse the string called "seas" and
; find out the starting and ending months of the
; season under consideration

; We have an array with all possible 3-month seasons
; These are arranged in the order that the first
; month of each season is also their index in the
; array. Then we can deduce the first and last months
; of the desired season.

seas=org_seas
seasons=(/"JFM","FMA","MAM","AMJ","MJJ","JJA",  \
          "JAS","ASO","SON","OND","NDJ","DJF"/)

; The code below will match the required season
; against the array and find the corresponding index
seas=changeCase(seas,"up")
iseas=ind(seasons.eq.seas)

; Thus we find the starting and ending months of the season
month1=iseas+1  ; starting month of the season
month2=month1+2 ; ending month of the season

;Find the temporal endpoints of the data
time=var&time
if (time(1)-time(0) .gt. 1 ) then
print("Time is not in daily resolution: You cannot use this function")
print(" E  X  I  T ")
exit
end if
ntim=dimsizes(time)
tunits=time@units

; Calculate value of first&last time index in Years, Months and Days
year1=ut_calendar(time(0),0)
year2=ut_calendar(time(ntim-1),0)
yr1=floattointeger(year1(:,0))
yr2=floattointeger(year2(:,0))
mo1=floattointeger(year1(:,1))
mo2=floattointeger(year2(:,1))
dy1=floattointeger(year1(:,2))
dy2=floattointeger(year2(:,2))

; If month1 = 11 or 12, we should shift yr1, 1 year back

if (month1.eq.11 .or. month1.eq.12) then
yr1 = yr1-1
end if

; Introduce a masking variable that has same time co-ordinates
; as the variable

tmask=new(ntim,"float")
add_dimensions(tmask,"time")
tmask&time=time
tmask=var@_FillValue
tmask@_FillValue=var@_FillValue

; For the required season, replace tmask with 1.0, so that
; the data is left undisturbed

; We treat the first and last years specially
; Check if starting time includes the last day of the final month
; of the required season
; If not do not operate

start_time=new(1,"double")
end_time=start_time
start_and_end(yr1,month1,tunits,start_time,end_time)

; Check if end_time is within time(0)
if (time(0).lt.end_time ) then

     t2=ind(time.eq.end_time)

; Check if start_time is within time(0)
  if (time(0).lt.start_time) then
     t1=ind(time.eq.start_time)
     tmask(t1:t2)=1.0
start_and_end(yr1-1,month1,tunits,start_time,end_time)
    t1old=doubletoint(end_time-time(0))
; Do the cosine bell fitting
     nlen=(t1-t1old+1)
     it1 = ispan(0,nlen-1,1)
     window=new(nlen,"float")
     window=cos(3.14*it1/nlen)
     window=fabs(window)
     tmask(0:t1)=(/window(0-t1old:nlen-1)/)
  else
; Otherwise start_time is time(0)
     t1=0
     tmask(t1:t2)=1.0
  end if
else
    t2=doubletoint(end_time-time(0))
end if

;     if (isdefined("it1")) then
;     delete(it1)
;     end if
;     if (isdefined("window")) then
;     delete(window)
;     end if

; For the years away from the end-points, the calculation is
; fairly simple

     t1old=t2
     
do iyear = yr1+1,yr2-1
	start_and_end(iyear,month1,tunits,start_time,end_time)
     t1=ind(time.eq.start_time)
     t2=ind(time.eq.end_time)
     tmask(t1:t2)=1.0
; Do the cosine bell fitting
     nlen=(t1-t1old+1)
     it1 = ispan(0,nlen-1,1)
     window=new(nlen,"float")
     window=cos(3.14*it1/nlen)
     window=fabs(window)
     if (t1old.lt.0) then
	tw=0
      else
	tw=t1old
      end if
     tmask(tw:t1)=(/window(tw-t1old:nlen-1)/)
     t1old=t2
     delete(it1)
     delete(window)
     
end do

; Now the last year
start_and_end(yr2,month1,tunits,start_time,end_time)

; Check if start_time is within time(ntim-1)
if (time(ntim-1).gt.start_time ) then

  t1=ind(time.eq.start_time)
; Do the cosine bell fitting
     nlen=(t1-t1old+1)
     it1 = ispan(0,nlen-1,1)
     window=new(nlen,"float")
     window=cos(3.14*it1/nlen)
     window=fabs(window)
     tw = t1old
     tmask(tw:t1)=(/window(tw-t1old:nlen-1)/)

; Check if end_time is within time(ntim-1)
  if (time(ntim-1).gt.end_time) then
     t2=ind(time.eq.end_time)
     tmask(t1:t2)=1.0
; Since we have more data available after the season, we have
; to taper it
start_and_end(yr2+1,month1,tunits,start_time,end_time)
t1old=t2
delete(it1)
delete(window)
     nlen=doubletoint(start_time-time(t1old))+1
     it1=new(nlen,"integer")
     it1 = ispan(0,nlen-1,1)
     window=new(nlen,"float")
     window=cos(3.14*it1/nlen)
     window=fabs(window)
     tmask(t1old:ntim-1)=(/window(0:ntim-1-t1old)/)


  else
; Otherwise end_time is time(ntim-1)
     t2=ntim-1
     tmask(t1:t2)=1.0
  end if
else
     nlen=doubletoint(start_time-time(t1old))+1
     it1=new(nlen,"integer")
     it1 = ispan(0,nlen-1,1)
     window=new(nlen,"float")
     window=cos(3.14*it1/nlen)
     window=fabs(window)
     tmask(t1old:ntim-1)=(/window(0:ntim-1-t1old)/)
end if

  return (tmask)
end

function mask_other_seasons(var,org_seas)
; this will designed as a function that when
; given a variable and a season will mask the
; values of the variable at times other than the
; specified season (with undefined values)

begin
seas=org_seas
; We should parse the string called "seas" and
; find out the starting and ending months of the
; season under consideration

; We have an array with all possible 3-month seasons
; These are arranged in the order that the first
; month of each season is also their index in the
; array. Then we can deduce the first and last months
; of the desired season.

seasons=(/"JFM","FMA","MAM","AMJ","MJJ","JJA",  \
          "JAS","ASO","SON","OND","NDJ","DJF"/)

; The code below will match the required season
; against the array and find the corresponding index
seas=changeCase(seas,"up")
iseas=ind(seasons.eq.seas)

; Thus we find the starting and ending months of the season
month1=iseas+1  ; starting month of the season
month2=month1+2 ; ending month of the season
;if (month2.gt.12) then
;month2=month2-12
;end if

;Find the temporal endpoints of the data
time=var&time
ntim=dimsizes(time)
tunits=time@units

; Calculate value of first&last time index in Years, Months and Days
year1=ut_calendar(time(0),0)
year2=ut_calendar(time(ntim-1),0)
yr1=floattointeger(year1(:,0))
yr2=floattointeger(year2(:,0))
mo1=floattointeger(year1(:,1))
mo2=floattointeger(year2(:,1))
dy1=floattointeger(year1(:,2))
dy2=floattointeger(year2(:,2))


; Introduce a masking variable that has same time co-ordinates
; as the variable

tmask=new(ntim,"float")
add_dimensions(tmask,"time")
tmask&time=time
tmask=var@_FillValue
tmask@_FillValue=var@_FillValue

; For the required season, replace tmask with 1.0, so that
; the data is left undisturbed

; We treat the first and last years specially
; Check if starting time includes the last day of the final month
; of the required season
; If not do not operate

start_time=new(1,"double")
end_time=start_time
start_and_end(yr1,month1,tunits,start_time,end_time)

; Check if end_time is within time(0)
if (time(0).lt.end_time ) then

     t2=ind(time.eq.end_time)

; Check if start_time is within time(0)
  if (time(0).lt.start_time) then
     t1=ind(time.eq.start_time)
     tmask(t1:t2)=1.0
  else
; Otherwise start_time is time(0)
     t1=0
     tmask(t1:t2)=1.0
  end if
end if

; For the years away from the end-points, the calculation is
; fairly simple

do iyear = yr1+1,yr2-1
	start_and_end(iyear,month1,tunits,start_time,end_time)
     t1=ind(time.eq.start_time)
     t2=ind(time.eq.end_time)
     tmask(t1:t2)=1.0
end do

; Now the last year
start_and_end(yr2,month1,tunits,start_time,end_time)

; Check if start_time is within time(ntim-1)
if (time(ntim-1).gt.start_time ) then

  t1=ind(time.eq.start_time)

; Check if start_time is within time(0)
  if (time(ntim-1).gt.end_time) then
     t2=ind(time.eq.end_time)
     tmask(t1:t2)=1.0
  else
; Otherwise start_time is time(0)
     t2=ntim-1
     tmask(t1:t2)=1.0
  end if
end if
  return (tmask)
end

procedure parse_time(time,year,month,day)
local tunits,tstring,year,mon_day,mon,day
begin
        tunits=time@units
        tstring=systemfunc("echo "+tunits+" | awk -Fsince '{print $2}'")
	print(tstring+":")
        year1=systemfunc("echo "+tstring+" | awk -F- '{print $1}'")
        mon_day=systemfunc("echo "+tstring+" | awk -F"+year1+"- '{print $2}'")
        month1=systemfunc("echo "+mon_day+" | awk -F- '{print $1}'")
        day1=systemfunc("echo "+mon_day+" | awk -F- '{print $2}'")
	year=stringtoint(year1)
	month=stringtoint(month1)
	day=stringtoint(day1)
end

procedure convert_time(oldtime,year2,month2,day2)

local oldtime,newtime
begin
mdays=(/31,28,31,30,31,30,31,31,30,31,30,31/)
; We assume that people write months as time since year-month-day1(1st day)


base_year=new(1,"integer")
base_month=new(1,"integer")
base_day=new(1,"integer")
parse_time(oldtime,base_year,base_month,base_day)
print(base_month+""+base_day)

if ( oldtime.and.isatt(oldtime,"calendar") ) then
calendar=oldtime@calendar
else
calendar="standard"
end if
print(calendar+"")

; These are for the newtime
tunits="days since 1900-01-01"
hour=0
minute=0
second=0

if (calendar.eq."gregorian" .or. calendar.eq."standard") then
oldtime@units="days since "+base_year+"-"+base_month+"-"+base_day
utc_date=ut_calendar(oldtime,0)
year2=floattointeger(utc_date(:,0))
month2=floattointeger(utc_date(:,1))
day2=floattointeger(utc_date(:,2))
end if

if (calendar.eq."noleap" \
.or. calendar.eq."365_day" \
.or. calendar.eq."no_leap") then
smonth=base_month-1
emonth=smonth+11
k=-1
nyrs=0
tim0=0.0
do while (k.le.-1)
nmos=0
  do imon=smonth,emonth
   tim0=tim0+mdays(imon)
      if (tim0.gt.oldtime) then
      k=0
	break
      end if
    nmos=nmos+1
   end do
nyrs=nyrs+1
end do

year2=base_year+nyrs-1
month2=base_month+nmos
day2=doubletoint(tim0-oldtime)

end if


if (calendar.eq."360_day" )  then
oldtime2=doubletoint(oldtime)
day_offset=oldtime2%30   ; check on what day the month is centered
nmos=floattoint((oldtime2-day_offset)/30.0) ; find out number of months
mo_offset=nmos%12  ; find if a full year is present; if not how many months
nyrs=(nmos-mo_offset)/12 ; how many years
year2=base_year+nyrs
month2=base_month+mo_offset
day2=base_day+day_offset
end if

;newtime=ut_inv_calendar(year2,month2,day2,hour,minute,second,tunits,0)
;return(newtime)
end

procedure make_month_time_in_std_format(time,yr1,yr2,mo1,mo2)

begin

;        (CALCULATE TIME AS DAYS SINCE 1900-01-01)
         
         time@units="days since 1900-01-01"
         time@calendar="standard"
         if isatt(yr1,"day")
          day=yr1@day
         else
           day=15
         end if
         hr=0
         minute=0
         sec=0
         k =0
	;
        nyr=yr2-yr1+1
        pp(nyr)
        if nyr.eq.1
          _mo2=mo2
        else
          _mo2=12
        end if
	do month=mo1,_mo2
         time(k)=ut_inv_calendar(yr1,month,day,hr,minute,sec,time@units,0)
         k=k+1
 	end do
	;
	if (yr1.ne.yr2) then
	;
       	 do year = yr1+1,yr2-1
          ;pp("here"+year)
       	 do month = 1,12
       	  time(k)=ut_inv_calendar(year,month,day,hr,minute,sec,time@units,0)
       	  k=k+1
       	 end do
       	 end do
	;
	 do month=1,mo2
       	  time(k)=ut_inv_calendar(yr2,month,day,hr,minute,sec,time@units,0)
	  k=k+1
	 end do
	;
	end if

end

function months_in_a_year()
begin
  months=(/"JAN","FEB","MAR","APR","MAY","JUN", \
         "JUL", "AUG","SEP","OCT","NOV","DEC"/)
  return(months)
end

function integer_value_of_month(org_month[1]:string)
local month
begin
  month=org_month
  months=(/"JAN","FEB","MAR","APR","MAY","JUN", \
         "JUL", "AUG","SEP","OCT","NOV","DEC"/)
  month=str_upper(month)
  return(ind(months.eq.month)+1)
end

function month_to_int(months)
begin
 return(get1Dindex(_months,str_upper(months))+1) 
end
function month_to_string(i_org_months)
begin
  i_months=i_org_months
  all_months=_months
  i_months=where(i_months .gt. 12, i_months - 12, i_months)
  return(all_months(i_months-1))
end

function month_to_name(month,opt)
begin 

        months1=(/"January","February","March","April","May","June", "July","August","September","October","November","December"/)

        month_chars=stringtochar(months1)
        months2=chartostring ( month_chars(:,:2) )
        months3=chartostring ( month_chars(:,:0) )

        
        if (opt.eq.0)
          return(months1(month-1)) 
          print(months1(month-1))
        end if
        if (opt.eq.1)
          return(months2(month-1))
        end if
        if (opt.eq.2)
          return(months3(month-1))
        end if
end


;       -------------------------------------

function integer_to_month(org_imon[1]:integer)
begin
  imon=org_imon
  months=(/"JAN","FEB","MAR","APR","MAY","JUN", \
         "JUL", "AUG","SEP","OCT","NOV","DEC"/)
  long_months=(/"JANUARY","FEBRUARY","MARCH","APRIL","MAY","JUNE", \
         "JULY", "AUGUST","SEPEMBER","OCTOBER","NOVEMBER","DECEMBER"/)
  if imon.gt.12
    imon=imon-12
  end if

  if isatt(imon,"longform")
    return(long_months(imon-1))
  else
    return(months(imon-1))
  end if
end
function integer_to_season(org_imon[1]:integer)
begin
  imon=org_imon
  months=_seasons
  if imon.gt.12
    imon=imon-12
  end if
return(months(imon-1))
end

function first_char_of_month(month[1]:string)
begin
  char_month=stringtochar(month)
  return(chartostring(char_month(0)))
end

function construct_season(start_season:integer,end_season:integer)
begin
  nmon=end_season-start_season+1
  seas=start_season
  myseason=first_char_of_month(integer_to_month(seas))
  do mon=1,nmon-1
    seas=seas+1
    myseason=myseason+first_char_of_month(integer_to_month(seas))
  end do
  return(myseason)
end

function array_to_string(var_arr)
begin
 dims=dimsizes(var_arr)
 rank=dimsizes(dims)
   if rank.gt.1
     print("This function is for 1-D arrays")
     exit
   end if
 var_str=""
if isatt(var_arr,"join")
spc=var_arr@join
else
spc=""
end if

   do i = 0,dims-2
     var_str=var_str+var_arr(i)+spc
   end do
   var_str=var_str+var_arr(dims-1)
return(var_str)
end

function array_join(arr,delimiter)
begin
  arr@join=delimiter
  return(array_to_string(arr))
end

function file_join(arr)
begin
  arr@join="/"
  return(array_to_string(arr))
end


function ChangeMissingValue(org_var,newMissVal)
begin
  var=org_var
  var=where(ismissing(var),newMissVal,var)
  var@_FillValue=newMissVal
  var@missing_value=newMissVal
  return(var)
end

function to_string(args)                                                        
begin                                                                           
  ndims=dimsizes(args)                                                          
  str_args=new(ndims,string)                                                    
  str_args=args
  return(str_args)
  ; the stuff below is unnecessary ;)
  if rank_of(args) .gt. 1                                                       
    terminate_program("Only accepts 1 D arrays")                                
  end if                                                                        
  ndims=dimsizes(args)                                                          
  str_args=new(ndims,string)                                                    
  do i = 0,ndims-1                                                              
    str_args(i) = ""+args(i)                                                    
  end do                                                                        
  return(str_args)                                                              
end
function numeral_to_string(var)
begin
return(systemfunc("echo "+var))
end

;function string_to_array(str)
;begin
;print(str)
;exit
;  if .not.isatt(str,"separator")
;    str@separator=""
;  end if
;  nstr=systemfunc("echo "+str+" | tr "+str@separator+" \\n | wc -l")
;  istr=stringtoint(nstr)
;  str_arr=new((/istr/),"string")
;
;  do i = 0, istr-1
;    CMD="echo "+str+" | awk -F"
;    CMD=CMD+str@separator+" '{print $"+(i+1)+"}'"
;    str_arr(i)=systemfunc(CMD)
;  end do
;  return(str_arr)
;end

function string_to_array(org_str)
begin
  str=org_str
  if .not.isatt(str,"separator")
    str@separator=","
    cmd="echo "+str+" | tr [:blank:] ',' "
    str=systemfunc(cmd)
  end if
  if str@separator .eq. " "
    sep="\ "
  else
    sep=str@separator
  end if
  istr = stringtoint(str_fields_count(str,sep))
  str_arr=new((/istr/),"string")

  do i = 0, istr-1
    CMD="echo "+str+" | awk -F"
    CMD=CMD+sep+" '{print $"+(i+1)+"}'"
    str_arr(i)=systemfunc(CMD)
  end do
  return(str_arr)
end

;function str_split(str,splitter)
; till we have string split in NCL
;begin
;  str@separator=splitter
;  return(string_to_array(str))
;end

function capitalize(org_aString)
begin
  aString=org_aString
  aString=changeCase(aString,"low")
  stringChars=stringtocharacter(aString)
  firstChar=systemfunc("echo "+stringChars(0)+" | tr a-z A-Z")
  return(firstChar+charactertostring(stringChars(1:)))
end

function Cel2Kel(org_var)
begin
  var=org_var
  var=var+273.15
  var@units="degK"
  return(var)
end

function make_time_in_std_format(yr1,yr2,mo1,mo2,dy1,dy2,ntim)
begin
 ; does not work properly
        terminate_program("does not work as expected: retired")
        tunits="days since 1900-01-01 00:00:00"
        start_time=ut_inv_calendar(yr1,mo1,dy1,0,0,0,tunits,0)
        end_time  =ut_inv_calendar(yr2,mo2,dy2,0,0,0,tunits,0)
        ;ntim=doubletoint(end_time-start_time+1)
        time=fspan(start_time,end_time,ntim)*1.d0
        time!0="time"
        time&time=time
        time@calendar="gregorian"
        time@units=tunits
        return(time)
end

function get_yr_mo_dy(time)
local time, tDat, year,month,day
begin
  tDat=ut_calendar(time,0)
  tDat@year=floattointeger(tDat(0,0))
  tDat@month=floattointeger(tDat(0,1))
  tDat@day=floattointeger(tDat(0,2))
return(tDat)
end

function to_std_time(oldtime)
begin
  ntim=dimsizes(oldtime)
  newtime=new(ntim,"double")
  newtime!0="time"
  ntunits = "days since 1900-01-01"
  newtime@units=ntunits
  do i = 0,ntim-1
    tmp=ut_calendar(oldtime(i),0)
    yr=to_integer(tmp(0,0))
    mo=to_integer(tmp(0,1))
    dy=to_integer(tmp(0,2))
    hr=to_integer(tmp(0,3))
    mn=to_integer(tmp(0,4))
    sc=to_integer(tmp(0,5))
    newtime(i)=ut_inv_calendar(yr,mo,dy,hr,mn,sc,ntunits,0)
  end do
  newtime&time=newtime
  newtime@calendar="standard"
  return(newtime)
end

function make_yearly_time(yr1,mo1,dy1,yr2,mo2,dy2,tunits,ntim)
begin
  newtime=new(ntim,"double")
  newtime!0="time"
  newtime@units=tunits
  do i = yr1,yr2
    newtime(i)=ut_inv_calendar(iyr,mo1,dy1,0,0,0,tunits,0)
  end do
  newtime&time=newtime
  newtime@calendar="standard"
return(newtime)
end

function make_monthly_time(yr1,mo1,dy1,yr2,mo2,dy2,tunits,ntim)
begin
  newtime=new(ntim,"double")
  newtime!0="time"
  newtime@units=tunits
  yr1@day=dy1
  make_month_time_in_std_format(newtime,yr1,yr2,mo1,mo2)
  newtime&time=newtime
  newtime@calendar="standard"
return(newtime)
end

function make_daily_time(yr1,mo1,dy1,yr2,mo2,dy2,tunits,ntim)
begin
  newtime=new(ntim,"double")
  newtime!0="time"
  newtime@units=tunits
  hr=0
  mn=0
  sc=0
  yr=yr1
  mo=mo1
  dy=dy1
  do i = 0,ntim-1
    newtime(i)=ut_inv_calendar(yr,mo,dy,hr,mn,sc,tunits,0)
    dy=dy+1
    if (dy.gt.days_in_month(yr,mo))
      dy=1
      mo=mo+1
    end if
    if (mo.gt.12)
      mo=1
      yr=yr+1
    end if
  end do
  newtime&time=newtime
  newtime@calendar="standard"
return(newtime)
end 

procedure fillMissingValues(var,dims,yindex)
local var,dims,yindex
begin
  rank=rank_of(var)
  if isatt(var,"opt")
    opt=var@opt
  else
    opt=-1
  end if

  ; create an array of dimensions
  ; move required dimension to the end
  newdims=dims
  newdims(yindex)=dims(rank-1)
  newdims(rank-1)=dims(yindex)
  if rank.eq.5
    var2=var($newdims(0)$|:,$newdims(1)$|:,$newdims(2)$|:,\
                            $newdims(3)$|:,$newdims(4)$|:)
    var2=linmsg(var2,opt)
    var=var2($dims(0)$|:,$dims(1)$|:,$dims(2)$|:,\
                                     $dims(3)$|:,$dims(4)$|:)
  end if
  if rank.eq.4
    var2=var($newdims(0)$|:,$newdims(1)$|:,$newdims(2)$|:,$newdims(3)$|:)
    var2=linmsg(var2,opt)
    var=var2($dims(0)$|:,$dims(1)$|:,$dims(2)$|:,$dims(3)$|:)
  end if
  if rank.eq.3
    var2=var($newdims(0)$|:,$newdims(1)$|:,$newdims(2)$|:)
    var2=linmsg(var2,opt)
    var=var2($dims(0)$|:,$dims(1)$|:,$dims(2)$|:)
  end if
  if rank.eq.2
    var2=var($newdims(0)$|:,$newdims(1)$|:)
    var2=linmsg(var2,opt)
    var=var2($dims(0)$|:,$dims(1)$|:)
  end if
end

procedure fillup_along_lon(var)
local dims,mydim,yindex
begin
  dims=getvardims(var)

  mydim="lon"
  yindex=ind(dims.eq.mydim)
  if ismissing(yindex)
    terminate_program("Missing the dimension of "+mydim)
  end if
  fillMissingValues(var,dims,yindex)
end

procedure fillup_along_lat(var)
local dims,mydim,yindex
begin
  dims=getvardims(var)

  mydim="lat"
  yindex=ind(dims.eq.mydim)
  if ismissing(yindex)
    terminate_program("Missing the dimension of "+mydim)
  end if
  fillMissingValues(var,dims,yindex)
end

procedure fillup_along_time(var)
local dims,mydim,yindex
begin
  dims=getvardims(var)

  mydim="time"
  yindex=ind(dims.eq.mydim)
  if ismissing(yindex)
    terminate_program("Missing the dimension of "+mydim)
  end if
  fillMissingValues(var,dims,yindex)
end

function next_time_of(time,incUnit)
local time, tdat, mon, year, day, incUnit
begin
  if incUnit.eq."day"
    return(time+1)
  end if
  if incUnit.eq."week"
    return(time+1)
  end if

  tdat=get_yr_mo_dy(time)
  mon=tdat@month
  year=tdat@year
  day=tdat@day
  if incUnit.eq."month"
   mon= mon + 1
   if mon.gt.12
     year=year+1
     mon=mon-12
    end if
  return(ut_inv_calendar(year,mon,day,0,0,0,time@units,0))
  end if
end
function next_month_of(time)
begin
 return(next_time_of(time,"month"))
end
function next_day_of(time)
begin
 return(next_time_of(time,"day"))
end
function next_week_of(time)
begin
 return(next_time_of(time,"week"))
end
function date2time(date_array,tunits)
; a trivial function, not well thought out
; restrains date to be specified as
; ddMMYYYY
begin
  ndates=dimsizes(date_array)
  tarray=new( ndates, double)
  do idat = 0,ndates-1
    datec=stringtochar(date_array(idat))
    dy=stringtoint(chartostring(datec(0:1)))
    mo=integer_value_of_month(chartostring(datec(2:4)))
    yr=stringtoint(chartostring(datec(5:8)))
    hr=0
    mn=0
    sc=0
    tarray(idat)=ut_inv_calendar(yr,mo,dy,hr,mn,sc,tunits,0)
  end do
return(tarray)
end
procedure any_undefined(var)
begin
  if any(ismissing(var))
    print("We have some missing values")
  end if
end
procedure exit_if_any_undefined(var)
begin
  if any(ismissing(var))
    print("We have some missing values")
    TerminateProgram
  end if
end
function get_obs_mon_year(start_mon,month,years)
begin
  if month .gt. 12
    needed_month=month-12
    needed_month@years=years+1
  else
    needed_month=month
    needed_month@years=years
  end if
  if month .lt. start_mon
    needed_month=month
    needed_month@years=years+1
  end if
  return(needed_month)
end
