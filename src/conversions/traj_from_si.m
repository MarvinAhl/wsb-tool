function traj = traj_from_si(traj_si, mdata)
    traj = traj_si;
    traj(:, 1:3) = traj(:, 1:3) / mdata.ulength_km;
    traj(:, 4:6) = traj(:, 4:6) / mdata.ulength_km * mdata.utime_s;
end