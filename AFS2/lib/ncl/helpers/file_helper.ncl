function fopen_read(fname)
begin
  return(addfile(fname,"r"))
end
function fopen_write(fname)
begin
  return(addfile(fname,"c"))
end
