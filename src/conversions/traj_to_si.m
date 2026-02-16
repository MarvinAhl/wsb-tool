function traj_si = traj_to_si(traj, mdata)
    traj_si = traj;
    traj_si(:, 1:3) = traj_si(:, 1:3) * mdata.ulength_km;
    traj_si(:, 4:6) = traj_si(:, 4:6) * mdata.ulength_km / mdata.utime_s;
end