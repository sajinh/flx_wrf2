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
a="I,am,fine"
print(string_to_array(a))
b="I am fine"
print(string_to_array(b))
c=" I am fine "
print(string_to_array(c))
