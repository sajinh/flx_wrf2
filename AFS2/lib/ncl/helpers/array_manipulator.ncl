function flatten_array(array)
begin
  igood=ind(.not. ismissing(array))
  return(array(igood))
end


procedure add_coordinates(var,axes)
begin
  ndims=dimsizes(axes)
  do i = 0,ndims-1
    axis=axes(i)
    units=axis+"_units"
    long_name=axis+"_long_name"
    var!i=axis
    tmp=axes@$axis$
    tmp!0=axis
    if isatt(axes,units)
      tmp@units=axes@$units$
    end if
    if isatt(axes,long_name)
      tmp@long_name=axes@$long_name$
    end if
    tmp&$axis$=tmp
    var&$axis$=tmp
    delete(tmp)
  end do
end
; + - +
procedure remove_unnecessary_attribute_of(var,att)
begin
 if isatt(var,att)
   delete(var@$att$)
 end if
end

function removeRedundantDimsOf(var)
begin
  var_dims=dimsizes(var)
  dim_names=getvardims(var)
  new_indcs=ind(var_dims.gt.1)
  newdims=var_dims(new_indcs)
  new_dimnames=(dim_names(new_indcs))
  ndims=dimsizes(newdims)
  newvar=new(newdims,typeof(var))
  newvar=var
  do i=0,ndims-1
    newvar!i=new_dimnames(i)
  end do
return(newvar)
end

function remove_redundant_dims_of(var)
begin
  return(removeRedundantDimsOf(var))
end

function remove_empty_dims_of(var)
begin
  return(removeRedundantDimsOf(var))
end

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

function slice_one(var,newdims,slice)
begin
  rank=rank_of(var)
  if rank.eq.5
    return(var($newdims(0)$|:,$newdims(1)$|:,$newdims(2)$|:,\
                            $newdims(3)$|:,$newdims(4)$|slice:slice))
  end if
  if rank.eq.4
    return(var($newdims(0)$|:,$newdims(1)$|:,$newdims(2)$|:,$newdims(3)$|slice:slice))
  end if
  if rank.eq.3
    return(var($newdims(0)$|:,$newdims(1)$|:,$newdims(2)$|slice:slice))
  end if
  if rank.eq.2
    return(var($newdims(0)$|:,$newdims(1)$|slice:slice))
  end if
end

function slice_many(var,newdims,slice)
begin
  rank=rank_of(var)
  if rank.eq.5
    return(var($newdims(0)$|:,$newdims(1)$|:,$newdims(2)$|:,\
                            $newdims(3)$|:,$newdims(4)$|slice))
  end if
  if rank.eq.4
    return(var($newdims(0)$|:,$newdims(1)$|:,$newdims(2)$|:,$newdims(3)$|slice))
  end if
  if rank.eq.3
    return(var($newdims(0)$|:,$newdims(1)$|:,$newdims(2)$|slice))
  end if
  if rank.eq.2
    return(var($newdims(0)$|:,$newdims(1)$|slice))
  end if
end

function exclude_from_last_dim(var,slice)
begin
  newdims=getvardims(var)
  rank=rank_of(var)
  if rank.eq.5
    slice2=get1Dindex_Exclude( var&$newdims(4)$, slice)
    return(var($newdims(0)$|:,$newdims(1)$|:,$newdims(2)$|:,\
                            $newdims(3)$|:,$newdims(4)$|slice2))
  end if
  if rank.eq.4
    slice2=get1Dindex_Exclude( var&$newdims(3)$, slice)
    return(var($newdims(0)$|:,$newdims(1)$|:,$newdims(2)$|:,$newdims(3)$|slice2))
  end if
  if rank.eq.3
    slice2=get1Dindex_Exclude( var&$newdims(2)$, slice)
    return(var($newdims(0)$|:,$newdims(1)$|:,$newdims(2)$|slice2))
  end if
  if rank.eq.2
    slice2=get1Dindex_Exclude( var&$newdims(1)$, slice)
    return(var($newdims(0)$|:,$newdims(1)$|slice2))
  end if
end

function cut_from_last_dim(var,slice,opt)
begin
  newdims=getvardims(var)
  if isatt(opt,"preserve_coordinates") .and. opt@preserve_coordinates
    if dimsizes(slice) .eq. 1
      return (slice_one(var,newdims,slice))
    else
      return (slice_many(var,newdims,slice))
    end if
  else
    return (slice_many(var,newdims,slice))
  end if
end
function ensure_last_dim(var,last_dim)
begin
  rank=rank_of(var)
  dims=getvardims(var)
  dindex=ind(dims.eq.last_dim)
  if dindex .eq. dimsizes(dims)-1
    return(var)
  end if
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
function reorder_dims_like_this(from,to)
begin
  org_dims=from
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
function average_along_given_dim(var,dim)
begin
  vclim = dim_avg_Wrap(ensure_last_dim(var,dim))
  return(reorder_dims_as_in_template(var,vclim))
end

function average_along_dim(var,dim,opt)
begin
  new_var=ensure_last_dim(var,dim)
  if isatt(opt,"preserve_coordinates") .and. opt@preserve_coordinates
    vclim = cut_from_last_dim(new_var,0,opt)
  end if 
  vclim = dim_avg_Wrap(new_var)
  return(reorder_dims_as_in_template(var,vclim))
end

function anomaly_along_dim(var,dim)
begin
  var1=ensure_last_dim(var,dim)
  vclm=dim_avg_Wrap(var1)
  rank=rank_of(var)
  if rank.eq.5
    vcli = conform_dims(dimsizes(var1), vclm, (/0,1,2,3/))
  end if
  if rank.eq.4
    vcli = conform_dims(dimsizes(var1), vclm, (/0,1,2/))
  end if
  if rank.eq.3
    vcli = conform_dims(dimsizes(var1), vclm, (/0,1/))
  end if
  if rank.eq.2
    vcli = conform_dims(dimsizes(var1), vclm, (/0/))
  end if
  vano = var1
  vano = vano - vcli
  return(reorder_dims_as_in_template(var,vano))
end

function remove_specific_dim(var,dim)
begin
  new_var=cut_from_last_dim(ensure_last_dim(var,dim),0,False)
  rank=rank_of(new_var)
  if rank .eq. 1
    return(new_var)
  end if
  return(reorder_dims_as_in_template(var,new_var))
end

function remove_specific_empty_dim(var, dim)
begin
  return(remove_specific_dim(var,dim))
end
function length_of(var)
begin
  dims=dimsizes(var)
  ndims=dimsizes(dims)
  length=1
  do i=0,ndims-1
    length=length*dims(i)
  end do
  return(length)
end
function split_at_last_dim(var)
begin
  dims=dimsizes(var)
  ndim=dimsizes(dims)
  last_dim=dims(ndim-1)
  blk_siz  = length_of(var)/last_dim
  return(onedtond(ndtooned(var), (/blk_siz,last_dim/)))
end
function reduce_to_2D(var,opt)
begin
  if opt .and. isatt(opt,"last_dim")
    new_var=ensure_last_dim(var,opt@last_dim)
  else
    new_var=var
  end if
  rank=rank_of(new_var)
  if rank .eq. 2
    pure_data=(/ new_var /)
    return(pure_data)
  end if
  return(split_at_last_dim(new_var))
end
function map_to_multidim(map,dims)
begin
  ndims=dimsizes(dims)
  nvar_dims=new(ndims,integer)
  do i=0,ndims-1
    nvar_dims(i)=dimsizes(dims@$dims(i)$)
  end do
  multi_map=onedtond(map,nvar_dims)
  add_coordinates(multi_map,dims)
  if isatt(dims,"reorder_dims")
    return(reorder_dims_like_this(dims@reorder_dims,multi_map)) 
  end if
  return(multi_map)
end

function does_vector_belong_to_set(vec,set)
  ; vector = 2d vector(x,y)
  ; set    = a collection of 2d vectors (/ (x1,y1), (x2,y2), .., (xn,yn)/)
begin
  dim_set=dimsizes(set)
  nset=dim_set(0)
  do i = 0,nset-1
    if all(set(i,:).eq.vec)
      return(True)
    end if
  end do
  return(False)
end

function common_elements_of_two_vectors(vec1,vec2)
begin
  nvec1=dimsizes(vec1)
  nvec2=dimsizes(vec2)
  str=""
  do i = 0,nvec1-1
    i2=ind(vec2.eq.vec1(i))
    if .not. ismissing(i2)
      str=str+i2+","
    end if
  end do
  str=chop(str)
  str@separator=","
  if str_is_blank(str)
    return(False)
  end if
  arr=stringtoint(string_to_array(str))
  return(vec2(arr))
end
