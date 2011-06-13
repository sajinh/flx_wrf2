function reformat_wrf_time(time,opt)
local time
begin
  split_time=str_split(time,"_")
  tim_left=split_time(0)
  tim_right=split_time(1)
  ymd=stringtoint(str_split(tim_left,"-")) ; year-mon-dy
  hms=stringtoint(str_split(tim_right,":")) ; hr:mi:se
  hoffset=9 ; 9 hrs
  tunits="days since 1900-01-01"
  new_time=ut_calendar(ut_inv_calendar(ymd(0),ymd(1),ymd(2),hms(0)+hoffset,hms(1),hms(2),tunits,0),0)
  ntime=str_join(new_time(:,:2),"-")+"---"+new_time(:,3)+":00hr"
  ntime@yr=new_time(:,0)
  ntime@mo=new_time(:,1)
  ntime@dy=new_time(:,2)
  ntime@hr=new_time(:,3)
  return(ntime)
end

