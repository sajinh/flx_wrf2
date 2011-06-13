
;   Example script to plot some 3D fields from a single metgrid file
;   Note the change in map and wind vector colors and the addition of PlotLevel info
;   November 2008

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
;
  a = addfile("../wrf-oper/run/met_em.d01.2010-07-20_12:00:00.nc","r")  ; Open a file


; We generate plots, but what kind do we prefer?
  type = "x11"
; type = "pdf"
; type = "ps"
; type = "ncgm"
  wks = gsn_open_wks(type,"plt_metgrid_6")


  res = True                           ; Set up some basic plot resources
  res@MainTitle = "METGRID FILES"
  res@Footer = False

  pltres = True
  mpres = True
  mpres@mpGeophysicalLineColor      = "Black"   ; Overwrite map settings
  mpres@mpGridLineColor             = "Black"
  mpres@mpLimbLineColor             = "Black"
  mpres@mpNationalLineColor         = "Black"
  mpres@mpPerimLineColor            = "Black"
  mpres@mpUSStateLineColor          = "Black"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; First get the variables we will need        
; Note: we do not have "tc", "ua", "va" in the input field - but we
;       know how to calculate them


    tc = wrf_user_getvar(a,"tc",0)  ; Calculate tc from TT 
    ua = wrf_user_getvar(a,"ua",0)  ; Get U on mass points
    va = wrf_user_getvar(a,"va",0)  ; Get V on mass points

    p  = wrf_user_getvar(a,"PRES",0) 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    dims = dimsizes(tc)          ; Get dims for tc
    nd = dimsizes(dims)          ; Find out how many dimensions there are
    nl = dims(nd-3)              ; We know 3rd dim from right is nz

    do level = 1,nl-1,5      ; LOOP OVER LEVELS, plot every 5th (surface is 0, lets start one level up)
 
       res@PlotLevelID = p(level,0,0)*0.01 + "hPa"    ; Add level info to plot
 
     ; Temperature
       opts = res
       opts@cnFillOn = True
       opts@ContourParameters = (/ .5 /)
       opts@lbBoxLinesOn = False           ; Move lines around label bar
       contour = wrf_contour(a,wks,tc(level,:,:),opts)
       delete(opts)
    
     ; Wind 
       opts = res
       opts@FieldTitle = "Wind"            ; Overwrite Field Title
       opts@NumVectors = 35                ; wind barb density - higher is more vectors
       opts@vcWindBarbColor = "White"      ; Draw wins barbs in white
       vector = wrf_vector(a,wks,ua(level,:,:),va(level,:,:),opts)
       delete(opts)

       plot = wrf_map_overlays(a,wks,(/contour, vector/),pltres,mpres)

    end do      ; END OF LEVEL LOOP

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

end
