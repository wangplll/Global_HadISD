; Calculate satuated vapor pressure
filename = "/home/wangplll/paper/Global_HPT/data3_preprocessing/tmean.nc"
ds = addfile(filename, "r")
t = ds->tmean

;Calculate satuated vapor pressure
esat_1  = satvpr_water_bolton(t, (/0,1/))   ; Pa

outname0 = "/home/wangplll/paper/Global_HPT/vapor pressure/esat_reshaped_global.nc"
es = addfile(outname0, "c")
es->esat = esat_1

