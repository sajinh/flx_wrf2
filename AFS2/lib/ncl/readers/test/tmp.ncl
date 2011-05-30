function reorder_dims(var,newdims)
begin
  rank=rank_of(var)
  if rank.eq.5
    return(var($newdims(0)$|:,$newdims(1)$|:,$newdims(2)$|:,\
                            $newdims(3)$|:,$newdims(4)$|:))
  end if
  if rank.eq.4
    return(var($newdims(0)$|:,$newdims(1)$|:,$newdims(2)$|:,$newdims(3)$|:))
  end if
  if rank.eq.3
    return(var($newdims(0)$|:,$newdims(1)$|:,$newdims(2)$|:))
  end if
  if rank.eq.2
    return(var($newdims(0)$|:,$newdims(1)$|:))
  end if
end

function ensure_last_dim(var,last_dim)
begin
  rank=rank_of(var)
  dims=getvardims(var)
  dindex=ind(dims.eq.last_dim)
  if ismissing(dindex)
    terminate_program("Missing the dimension of "+last_dim)
  end if

  newdims=dims
  newdims(dindex)=dims(rank-1)
  newdims(rank-1)=dims(dindex)
  return(reorder_dims(var,newdims))
end

function reorder_dims_as_in_template(from,to)
begin
  org_dims=getvardims(from)
  dims=getvardims(to)

  re_dims=new( dimsizes(dims), typeof(dims))
  idim=0
  do i = 0, dimsizes(org_dims) - 1
    if isdim(to,org_dims(i))
       re_dims(idim) = org_dims(i)
       idim=idim+1
    end if
  end do
  return(reorder_dims(to,re_dims))
end

