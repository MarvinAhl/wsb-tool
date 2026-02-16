clear variables
close all
clc

% Add src and all subfolders to path
addpath(genpath("../src"));

mdata = init_model_data(false);
sdata = init_sim_data(mdata);

% Generate grid of raans and eccentricities for capture check
N_raans = 8;
N_eccs = 8;
raans = linspace(0, 2*pi, N_raans+1);
raans = raans(1:end-1);
eccs = linspace(0.92, 1, N_eccs+1);
eccs = eccs(1:end-1);

% Propagation time
tf_back_si = 3600 * 24 * 60;  % s, 60 days to see how far they go
tf_back = tf_back_si / mdata.utime_s;
tf_forw_si = 3600 * 24 * 10;  % s 10 day for one orbit seems enough
tf_forw = tf_forw_si / mdata.utime_s;
dt_si = 1 * 60;  % s
dt = dt_si / mdata.utime_s;
tspan_back = 0 : -dt : -tf_back;
tspan_forw = 0 : dt : tf_forw;
tspan_full = [tspan_back(end:-1:1), tspan_forw];

% For saving the trajectories
trajs = zeros(N_raans, N_eccs, length(tspan_full), 6);
traj_success = zeros(N_raans, N_eccs);

% Capture height
h_c = 600;  % Optimally this would be 100 but that would be too dangerous (don't wanna crash into the moon)
% min height for second periapsis
h_c_p2_min = 100;
h_c_p2_max = 5000;

% Loops to compute trajectories
wb = waitbar(0, "Grid search (0%)");
wb_N_tot = N_raans * N_eccs;
for raan_idx = 1 : N_raans
    for ecc_idx = 1 : N_eccs
        wb_N = (raan_idx-1) * N_eccs + ecc_idx;
        waitbar(wb_N / wb_N_tot, wb, "Grid search (" + floor(wb_N / wb_N_tot * 100) + "%)");

        % Generate initial keplerian state
        rp = h_c + mdata.R_m;
        kep0.e = eccs(ecc_idx);
        kep0.a = rp / (1 - kep0.e);
        kep0.inc = deg2rad(90);
        kep0.raan = raans(raan_idx);
        kep0.aop = pi/2;
        kep0.ta = 0;

        % Transform to cartesian coordinates
        xx0 = kep2cars(kep0, mdata.GM_m);

        % Transform to em frame
        xx0 = state_from_si(xx0, mdata);
        xx0 = xx0 + origin_change(0, "in", "m", "em", mdata);
        xx0 = frame_change(0, "in", "em", mdata) * xx0;

        % Basic idea: Propagate forward to periapsis (to make sure that
        % it's actually captured) and backward until SOI exit to make sure
        % that it escapes. Then test if it's inside or outside the moon
        % radius

        % Propagate forward and check if another periapsis is reached
        ode_opts = odeset("RelTol", sdata.rel_tol, "AbsTol", sdata.abs_tol, "Events", @(t, xx) event_periapsis(t, xx, "m", mdata, false));
        [~, traj_forw, t_ev_p, xx_ev_p] = ode113(@(t, xx) BCRFBP_dyn(t, xx, mdata), tspan_forw, xx0, ode_opts);
        % Save trajectory
        trajs(raan_idx, ecc_idx, end-length(tspan_forw)+1:end, :) = traj_forw;

        % Indicates if the orbit is a capture orbit
        is_capture = ~isempty(t_ev_p) && sum(t_ev_p > 2*dt);
        % Makes sure that periapsis pass is in right altitude range
        % (between 100km and 800km)
        if is_capture
            peri_idxs = find(t_ev_p > 2*dt);
            peri_1 = peri_idxs(1);
            d_m = norm(xx_ev_p(peri_1, 1:3)' - [1-mdata.mu; 0; 0]) * mdata.ulength_km;
            
            if d_m > (h_c_p2_max + mdata.R_m) || d_m < (h_c_p2_min + mdata.R_m)
                is_capture = false;
            end
        end

        % Propagate backward and check if the orbit escapes
        ode_opts = odeset("RelTol", sdata.rel_tol, "AbsTol", sdata.abs_tol, "Events", @(t, xx) event_escape(t, xx, "m", mdata.GM_m, mdata.r_soi_m, mdata, false));
        [~, traj_back, t_ev_e, xx_ev_e] = ode113(@(t, xx) BCRFBP_dyn(t, xx, mdata), tspan_back, xx0, ode_opts);
        % Save trajectory
        trajs(raan_idx, ecc_idx, 1:length(tspan_back), :) = traj_back(end:-1:1, :);

        is_escape = ~isempty(t_ev_e);

        % Check if escapes to inside or outside of moon orbit
        is_outside_moon = false;
        if is_escape
            % Transform to earth centred
            xx_ev_e_in = frame_change(t_ev_e(1), "em", "in", mdata) * xx_ev_e(1, :)';
            xx_ev_e_in = xx_ev_e_in + origin_change(t_ev_e(1), "in", "em", "e", mdata);
            xx_ev_e_in = state_to_si(xx_ev_e_in, mdata);

            is_outside_moon = norm(xx_ev_e_in(1:3)) > mdata.ulength_km;
        end

        % Counts as successful ballistic capture if these conditions all
        % hold
        traj_success(raan_idx, ecc_idx) = is_capture && is_escape && is_outside_moon;
    end
end

close(wb);

% Print successful captures
disp(" ");
disp("Successful captures for index pairs:")
disp("RAAN idx    ECC idx")
[succ_raan_idxs, succ_ecc_idxs] = find(traj_success == 1);
for ii = 1:length(succ_raan_idxs)
    disp(succ_raan_idxs(ii) + "           " + succ_ecc_idxs(ii));
end


%% Plot a single one of the trajectories
close all

% Choose which one to plot here:
raan_idx = 7;
ecc_idx = 4;

disp(" ");
disp("Plotting trajectory with indexes (" + raan_idx + ", " + ecc_idx + ")");
disp("RAAN = " + raans(raan_idx) + "    ECC = " + eccs(ecc_idx));

[Xsph, Ysph, Zsph] = sphere;

R_e_3b = mdata.R_e / mdata.ulength_km;
R_m_3b = mdata.R_m / mdata.ulength_km;

moon_orbit = zeros(100, 3);
for ii = 1 : 100
    moon_orbit(ii, 1) = mdata.ulength_km * cos(2*pi / 100 * ii);
    moon_orbit(ii, 2) = mdata.ulength_km * sin(2*pi / 100 * ii);
end

traj = reshape(trajs(raan_idx, ecc_idx, :, :), length(tspan_full), 6);
% Transform to moon inertial and (later) se rotating frame
traj_in_m = frame_change_traj(tspan_full, traj, "em", "in", mdata);
traj_in_m = origin_change_traj(tspan_full, traj_in_m, "in", "em", "m", mdata);
traj_in_m = traj_to_si(traj_in_m, mdata);

traj_se = frame_change_traj(tspan_full, traj, "em", "se", mdata);
traj_se = origin_change_traj(tspan_full, traj_se, "se", "em", "e", mdata);
traj_se = traj_to_si(traj_se, mdata);

E2b_m = zeros(length(tspan_full), 1);
for ii = 1 : length(tspan_full)
    xx_in_m = traj_in_m(ii, :)';
    E2b_m(ii) = E_2B(xx_in_m, mdata.GM_m);
end

J = zeros(length(tspan_full), 1);
for ii = 1 : length(tspan_full)
    xx = traj(ii, :)';
    J(ii) = J_CRTBP(mdata.mu, xx(1), xx(4), xx(2), xx(5), xx(3), xx(6));
end

% Plot in Moon centred inertial frame
f_in_m = figure();
plot3(traj_in_m(:, 1), traj_in_m(:, 2), traj_in_m(:, 3), "k-", "LineWidth", 1);
hold on;
% Plot moon
surf(Xsph * mdata.R_m, Ysph * mdata.R_m, Zsph * mdata.R_m, 'FaceColor', '#ffffff');
% Adjust plot
axis equal;
grid on;
xlim([-5e4, 5e4]);
ylim([-5e4, 5e4]);
zlim([-1e5, 5e4]);
%title("Moon-Centred Inertial Frame");
xlabel("x / km");
ylabel("y / km");
zlabel("z / km");

% Plot in Earth centred sun-earth rotating frame
f_se = figure();
plot3(traj_se(:, 1), traj_se(:, 2), traj_se(:, 3), "k-", "LineWidth", 1);
hold on;
% Plot L1 and L2
plot3(-mdata.L1_se, 0, 0, "k+", "LineWidth", 1);  % L1
plot3(mdata.L2_se, 0, 0, "k+", "LineWidth", 1);  % L2
% Plot earth
surf(Xsph * mdata.R_e, Ysph * mdata.R_e, Zsph * mdata.R_e, 'FaceColor', '#ffffff');
% Plot moon orbit
plot3(moon_orbit(:, 1), moon_orbit(:, 2), moon_orbit(:, 3), "k--", "LineWidth", 1);
% Adjust plot
axis equal;
grid on;
xlim([-mdata.L1_se, mdata.L2_se] * 1.1);
%ylim([min(traj_se(:, 2)), max(traj_se(:, 2))] * 1.2);
title("Sun-Earth Rotating Frame (Earth-Centred)")
xlabel("x / km");
ylabel("y / km");

% Plot Moon 2 body energy
f_E2b = figure();
plot(tspan_full, E2b_m, "k-", "LineWidth", 1);
%hold on;
grid on;
title("Moon 2-Body Energy");
xlabel("t / s");
ylabel("E_m / km^2/s^2");

% Plot Jacobi constant
f_J = figure();
plot(tspan_full, J, "k-", "LineWidth", 1);
%hold on;
grid on;
title("Jacobi Energy")
xlabel("t / s");
ylabel("J");