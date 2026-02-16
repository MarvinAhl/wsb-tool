function sim_data = init_sim_data(mdata)
    sim_data.dt = 3600 / mdata.utime_s;  % 1h in dimensionless units
    sim_data.rel_tol = 1e-9;
    sim_data.abs_tol = 1e-10;
end