function xx = state_from_si(xx_si, mdata)
    xx = xx_si;
    xx(1:3) = xx(1:3) / mdata.ulength_km;
    xx(4:6) = xx(4:6) / mdata.ulength_km * mdata.utime_s;
end