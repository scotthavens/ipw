load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


begin

;---------------------------------------------------------------------
;---Input file, weight file, and grid information (WRF data)
;---------------------------------------------------------------------

srcDirName  = "/Volumes/Drobo1/SCOTT/WRF/2009_1km/"
srcFileName = "wrfout_d03_WY2009.nc"
srcFilePath =  srcDirName + srcFileName
sfile = addfile(srcFilePath,"r")

;---------------------------------------------------------------------
;---Set up NetCDF output file
;---------------------------------------------------------------------

ofile = "test.nc"

system("/bin/rm -f " + ofile)                   ; remove any old version
setfileoption("nc","Format","NetCDF4Classic")
ncdf  = addfile(ofile,"c")                      ; create a new NetCDF
filedimdef(ncdf,"Time",-1,True)                 ; define the file dimensions

;--- Determine number of time steps and other dimension sizes
time  = wrf_user_getvar(sfile,"Times",-1)       ; get times in the file
times = time(1:10,:)                            ; Adjust for boundary condition
ntimes = dimsizes(times)                        ; number of times in the file

print("Time steps = " + ntimes(0))

;--- file definition mode
setfileoption(ncdf,"DefineMode",True)

;--- create global atributes of the file
fAtt               = True            ; assign file attributes
fAtt@title         = "Test NetCDF File"
fAtt@source_file   =  srcFileName
fAtt@Conventions   = "None"
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes

;--- predefine the coordinate variables and their dimensionality
dimNames = (/"Time", "DateStrLen", "south_north", "west_east"/)
dimSizes = (/ -1   , 19, 240,  345/)
dimUnlim = (/ True , False, False, False/)
filedimdef(ncdf,dimNames,dimSizes,dimUnlim)

;--- Time variable
filevardef(ncdf, "Times" ,typeof(time),(/"Time","DateStrLen"/))
filevarattdef(ncdf,"Times" ,time)                    ; copy time attributes
setfileoption(ncdf,"DefineMode",False)
ncdf->Times = times

vdimNames = (/"Time", "south_north", "west_east"/)

;---------------------------------------------------------------------
;--- VARIABLES DESIRED
;---------------------------------------------------------------------

;--- Terrain height, actually only need one of them
print("Getting variable HGT")

var = wrf_user_getvar(sfile,"HGT",1)
filevardef(ncdf, "HGT" ,typeof(var),(/"south_north", "west_east"/))       ; define dimensions
filevarattdef(ncdf,"HGT" ,var)                      ; copy attributes

ncdf->HGT = (/ var /)        ; output only the values
delete(var)


print("----------------------------")

;--- Get all the other variables in v

v = (/ "GLW", "T2", "rh2" /)
;v = (/ "GLW" /)
ds = dimsizes(v)

do i = 0, ds(0)-1

    print("Getting variable " + v(i))

    do it = 0, ntimes(0)-1
    var = wrf_user_getvar(sfile,v(i),it+1)

    if (it .eq. 0) then
        filevardef(ncdf, v(i) ,typeof(var),vdimNames)       ; define dimensions
        filevarattdef(ncdf, v(i) ,var)                      ; copy attributes
    end if

    ;var_regrid = ESMF_regrid_with_weights(var,wfile,Opt)

    ncdf->$v(i)$(it,:,:) = (/ var /)

    delete(var)

end do        ; END OF TIME LOOP

print("----------------------------")

end do          ; END OF VARIABLE LOOP

end

