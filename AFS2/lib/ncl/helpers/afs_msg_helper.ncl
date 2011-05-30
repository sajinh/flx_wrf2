function afs_top_signature()
begin
  return("_____afs=starting=up_____")
end
function afs_bottom_signature()
begin
  return("_____afs=winding=down_____")
end
procedure start_script(msg)
begin
  pp(afs_top_signature)
  pp(msg)
end
procedure stop_script(msg)
begin
  pp(msg)
  pp(afs_bottom_signature)
  exit
end
procedure afs_abort(msg)
begin
  pp("abort::"+msg)
  exit
end

