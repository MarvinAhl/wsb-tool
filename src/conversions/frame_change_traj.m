function traj_f2 = frame_change_traj(tspan, traj_f1, frame1, frame2, mdata)
    traj_f2 = zeros(size(traj_f1, 1), size(traj_f1, 2));
    for ii = 1 : length(tspan)
        traj_f2(ii, :) = (frame_change(tspan(ii), frame1, frame2, mdata) * traj_f1(ii, :)')';
    end
end