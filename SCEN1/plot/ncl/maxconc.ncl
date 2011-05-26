load "$AFS2/lib/ncl/helper_libs.ncl"
load "./conc_plot_helpers.ncl"

begin
  nx=36
  ny=68
  nz=43
  nt=17
  read_flxout_dims(nx,ny,nz,nt)
  path="../conc.dat"

  conc = cbinread( path,  (/nt,nz,ny,nx/), "float")
  mconc=max(conc(0,:,:,:))

  path="../drydep.dat"
  ddep = cbinread( path,  (/nt,ny,nx/), "float")
  mddep=max(ddep(0,:,:))
  
  path="../wetdep.dat"
  wdep = cbinread( path,  (/nt,ny,nx/), "float")
  mwdep=max(wdep(0,:,:))

  print(mconc+" "+mddep+" "+mwdep)
  max_conc=(mconc+mddep+mwdep)
  print(max_conc)
  asciiwrite("maxconc.dat",max_conc)
end


