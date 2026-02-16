function xx_si = state_to_si(xx, mdata)
    xx_si = xx;
    xx_si(1:3) = xx_si(1:3) * mdata.ulength_km;
    xx_si(4:6) = xx_si(4:6) * mdata.ulength_km / mdata.utime_s;
end