; ACE2-ERA5 vertical interpolation to pressure levels using vinth2p
; Author: Assistant (test version for single member)

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"

begin

  ;--- Set output format
  setfileoption("nc", "Format", "NetCDF4Classic")
  dates = (/"01-01", "01-11", "01-21", "01-31", "02-10", "02-20", "03-01", "03-11", "03-21", "03-31", "04-10", "04-20", "04-30", "05-10", "05-20", "05-30", "06-09", "06-19", "06-29", "07-09", "07-19", "07-29", "08-08", "08-18", "08-28", "09-07", "09-17", "09-27", "10-07", "10-17", "10-27", "11-06", "11-16", "11-26", "12-06", "12-16", "12-26"/)
  N = 37
  do m = 0, N-1  ; Loop over members (adjust upper bound if more members are available)

    i = m
    ;--- Load input files
      finT = addfile("air_temperature_1980-"+dates(i)+".nc", "r")   ; contains T (member_id,time,level,lat,lon)


    ;--- Read variables
    Tin  = finT->air_temperature       ; (member_id,time,level,lat,lon)
    PS   = finT->PS                    ; (member_id,time,lat,lon)
    ak   = finT->ak
    bk   = finT->bk
    time = finT->time
    lat  = finT->lat
    lon  = finT->lon
    tbot = finT->tbot
    phis = finT->zsfc

    ;hyam = 0.5 * (ak(0:6) + ak(1:7))
    ;hybm = 0.5 * (bk(0:6) + bk(1:7))

    hyam = 0.5 * (ak(0:7) + ak(1:8))
    hybm = 0.5 * (bk(0:7) + bk(1:8))

    ;hyam = ak(1:8)
    ;hybm = bk(1:8)

    printVarSummary(hyam)
    printVarSummary(hybm)
    ;--- Set constants
    P0mb   = 1000.0     ; reference pressure in Pa
    hyam   = hyam / (P0mb * 100.)
    interp = 1            ; linear interpolation
    extrap = True        ; allow extrapolation
    varflg = 1            ; variable to interpolate

    pnew = (/   1.0,    2.0,    3.0,    5.0,    7.0,   10.0,   20.0,   30.0,   50.0,   70.0, 100.0,  125.0,  150.0,  175.0,  200.0,  225.0,  250.0,  300.0,  350.0,  400.0, 450.0,  500.0,  550.0,  600.0,  650.0,  700.0,  750.0,  775.0,  800.0,  825.0, 850.0,  875.0,  900.,  925.,  950.0,  975.0, 1000.0 /)  ; Target pressure levels (Pa)

    print("Interpolating member: " + m)

    ; Extract data
    temp_m = Tin(:,1:,:,:)    ; (time,level,lat,lon)
    ps_m   = PS(:,:,:)       ; (time,lat,lon)
    tbot   = tbot(:,:,:)
    zsfc   = phis(:,:,:)
    
    ; Diagnostics
    ;printVarSummary(temp_m)
    ;printVarSummary(ps_m)

    ; Interpolation from hybrid-sigma to pressure levels
    Tm_interp = vinth2p_ecmwf(temp_m, hyam, hybm, pnew, ps_m, interp, P0mb, 1, extrap, varflg, tbot, phis)  ; (time,plev,lat,lon)
    printVarSummary(Tm_interp)

    ; Output file name
    out_file = "/glade/derecho/scratch/katyr/ai-models/ACE2/ACE2_37_Interpolation_t/ace2-era5.member" + m + ".air_temperature_month.1980-2022.plevs.nc"

    ; Delete existing file if present
    if (isfilepresent(out_file)) then
      system("rm -f " + out_file)
      print("Existing file " + out_file + " deleted.")
    end if

    ; Open NetCDF file for writing
    ncid = addfile(out_file, "c")

    ; Get dimension sizes
    ntime = dimsizes(Tm_interp(:,0,0,0))
    nlev  = dimsizes(pnew)
    nlat  = dimsizes(Tm_interp(0,0,:,0))
    nlon  = dimsizes(Tm_interp(0,0,0,:))

    ; Define dimensions
    filedimdef(ncid, "time", -1, True)      ; unlimited
    filedimdef(ncid, "lev_p", nlev, False)
    filedimdef(ncid, "lat", nlat, False)
    filedimdef(ncid, "lon", nlon, False)

    ; Define coordinate variables
    filevardef(ncid, "time", "float", (/"time"/))
    filevardef(ncid, "lev_p", typeof(pnew), (/"lev_p"/))
    filevardef(ncid, "lat", typeof(lat), (/"lat"/))
    filevardef(ncid, "lon", typeof(lon), (/"lon"/))

    ; Define main variable
    filevardef(ncid, "T", typeof(Tm_interp), (/"time", "lev_p", "lat", "lon"/))

    ; Assign coordinate values
    ;ncid->time  = tofloat(ispan(0, ntime-1, 1))  ; fake time axis
    ;ncid->lev_p = pnew
    ;ncid->lat   = lat
    ;ncid->lon   = lon

    ; Assign temperature variable
    ncid->T = Tm_interp
    ncid->T@long_name = "Air Temperature"
    ncid->T@units = "K"
    ncid->T@_FillValue = Tm_interp@_FillValue

    ; Close file
    delete(ncid)

  end do
end
