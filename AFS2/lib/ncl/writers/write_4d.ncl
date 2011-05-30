procedure write_TDLL(fout,vnam,var) 
begin
lat=var&lat
lon=var&lon
depth=var&depth
nlat=dimsizes(lat)
nlon=dimsizes(lon)
nlev=dimsizes(depth)
time=var&time
    ;===================================================================
    ; create global attributes of the file
    ;===================================================================
        fAtt               = True            ; assign file attributes
	fAtt@title         = "NCL Efficient Approach to netCDF Creation"  
	fAtt@Conventions   = "None"   
	fAtt@creation_date = systemfunc ("date")        
	fileattdef( fout, fAtt )            ; define file attributes    
    ;===================================================================
    ; predefine the coordinate variables and their dimensionality
    ; note: to get an UNLIMITED record dimension, we set the dimensionality
    ; to -1 and set the unlimited array to True.
    ;===================================================================
        dimNames = (/"time","depth","lat","lon"/)  
	dimSizes = (/-1,nlev,nlat,nlon/) 
	dimUnlim = (/True,False,False,False/)   
	filedimdef(fout,dimNames,dimSizes,dimUnlim)
    ;===================================================================
    ; predefine the the dimensionality of the variables to be written out
    ; filevardef(output_file, name of variable, type of variable, CV's of 
    ; variable) note: at this point, when you do an ncdump on the file, it 
    ; will look like the variable is there. If you do not output the values,
    ; however, then an ncdump -v xxx will reveal that there is no data.
    ;===================================================================
       filevardef(fout, "time" ,typeof(time),"time")
       filevardef(fout, "depth"  ,typeof(depth),"depth")
       filevardef(fout, "lat"  ,typeof(lat),"lat")
       filevardef(fout, "lon"  ,typeof(lon),"lon")
       filevardef(fout, vnam   ,typeof(var), (/"time","depth","lat","lon"/))
    ;===================================================================
    ; add some attributes
    ; note, since we are copying variables from one file to another,
    ; the original variables already have attributes, and we can copy
    ; them directly by just passing the variable to filevarattdef. You can 
    ; create your own attributes by assigning them to the variable as we
    ; do with u.
    ;====================================================================
       filevarattdef(fout,vnam,var)
       filevarattdef(fout,"lat",lat)
       filevarattdef(fout,"lon",lon)
       filevarattdef(fout,"depth",depth)
       filevarattdef(fout,"time",time)
    ;===================================================================
    ; output only the data values since the dimensionality and such have
    ; been predefined. the "(/", "/)" syntax tells NCL to only output the
    ; data values.
    ;====================================================================
       fout->time   = (/time/)     
       fout->depth   = (/depth/)     
       fout->lat    = (/lat/)
       fout->lon    = (/lon/) 
       fout->$vnam$      = (/var/)

end
