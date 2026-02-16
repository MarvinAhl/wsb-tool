function traj_o2 = origin_change_traj(tspan, traj_o1, frame, origin1, origin2, mdata)
    traj_o2 = zeros(size(traj_o1, 1), size(traj_o1, 2));    
    for ii = 1 : length(tspan)
        traj_o2(ii, :) = traj_o1(ii, :) + origin_change(tspan(ii), frame, origin1, origin2, mdata)';
    end
end