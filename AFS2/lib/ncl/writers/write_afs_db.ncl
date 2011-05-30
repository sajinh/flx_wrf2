load "$AFS2/lib/ncl/writers/writem.ncl"

function rank_check(var)
begin
  rank=rank_of(var)
  if rank.lt.3 .or. rank.gt.4
    terminate_program("rank < 3 OR  > 4: cannot process")
  end if
  return(rank)
end

function make_outdir(db_root,lead,start_mon,store,cabinet,model,start_year,opt)
begin
  root_dir=db_root
  outdir=lead+"-MON"+"/"+store
  if isatt(opt,"season") .and. str_upper(opt@season) .ne. "NONE"
    mon2=start_mon+"/"+opt@season
  else
    mon2=start_mon
  end if
  if isatt(opt,"mme_ver")
    if .not. ismissing(opt@mme_ver)
      mme_ver =opt@mme_ver
      outdir=outdir+"/"+"VER-"+mme_ver+"/"
    end if
  end if
  outdir=outdir+"/"+cabinet+"/"+model+"/"+mon2+"/"+start_year
  outdir=root_dir+"/"+str_upper(outdir)
  return(outdir)
end

function mme_out_dims_match(var,rank)
begin
  vardims=getvardims(var)
  dims3d=(/"time", "lat", "lon"/)  
  dims4d= (/"time", "level", "lat", "lon"/)  
  if rank .eq. 3
    match=all(vardims .eq. dims3d)
    match@dims=dims3d
  end if
  if rank .eq. 4
    match=all(vardims .eq. dims4d)
    match@dims=dims4d
  end if
  return(match)
end

function afs_db_writer(db_root,var,vname,lead,store,cabinet,model,start_mon,start_year,opt)
begin
  rank=rank_check(var)
  matching_dims = mme_out_dims_match(var,rank)
  if .not. matching_dims
    terminate_program("Output Dimensions must be "+matching_dims@dims)
  end if
  outdir=make_outdir(db_root,lead,start_mon,store,cabinet,model,start_year,opt)
  outfil=outdir+"/"+vname+".nc"


  if isfilepresent(outfil)
    if isatt(opt, "delete_output") .and. (opt@delete_output .eq. True)
      system("rm -rf "+outfil)
    else
      afs_abort("output file "+cr+outfil+cr+" exists : please delete to continue")
    end if
  end if

  if .not.isfilepresent(outdir)
    system("mkdir -p "+outdir)
  end if

  fout=addfile(outfil,"c")
  if str_upper(store) .eq. "MME_OUT"
    fout@mme_method=model
  else
    fout@model=model
  end if
  file_atts = getvaratts(model)
  if(.not.any(ismissing(file_atts))) then
    do i=0,dimsizes(file_atts)-1
      fout@$file_atts(i)$=model@$file_atts(i)$
    end do
  end if

  lat=var&lat
  lon=var&lon
  nlat=dimsizes(lat)
  nlon=dimsizes(lon)
  month=integer_value_of_month(start_mon) 
  tunits=AFS_TimeUnits()
  if isatt(opt,"season") .and. (str_upper(opt@season) .ne. "NONE")
    time=ut_inv_calendar(start_year,month,15,0,0,0,tunits,0)
  else
    time=new(lead,double)
    do i = 0,lead-1
      time(i)=ut_inv_calendar(start_year,month+i,15,0,0,0,tunits,0)
    end do
  end if

  ntim=dimsizes(time)

  if rank.eq.4
    nlev=dimsizes(var&level)
    vdat=new( (/ntim,nlev,nlat,nlon/), typeof(var), var@_FillValue)
    add_dimensions(vdat,(/"time","level","lat","lon"/))
    vdat&level=var&level
    vdat&time=time
    vdat&lat=var&lat
    vdat&lon=var&lon
    vdat= (/ var /)
    copy_VarAtts(var,vdat)
    write_TLLL(fout,vname,vdat)
  else
    vdat=new( (/ntim,nlat,nlon/), typeof(var), var@_FillValue)
    add_dimensions(vdat,(/"time","lat","lon"/))
    vdat&time=time
    vdat&lat=var&lat
    vdat&lon=var&lon
    vdat= (/ var /)
    copy_VarAtts(var,vdat)
    write_TLL(fout,vname,vdat)
  end if
  return(outfil)
end
