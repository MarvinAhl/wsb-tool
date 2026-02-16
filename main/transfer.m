clear variables
close all
clc

% Add src and all subfolders to path
addpath(genpath("../src"));

mdata = init_model_data(true);
sdata = init_sim_data(mdata);

% Ballistic capture at time 0 moon centred inertial
% INSERT VALUES FROM capture.m HERE
kep0.raan = 4.7124;
kep0.e = 0.95;

% Capture altitude
h_c = 600;

rp = h_c + mdata.R_m;  % 600km periapsis;
kep0.a = rp / (1 - kep0.e);
kep0.inc = deg2rad(90);
kep0.aop = pi/2;
kep0.ta = 0;

% Transform to em frame
xx0 = kep2cars(kep0, mdata.GM_m);
xx0 = state_from_si(xx0, mdata);
xx0 = xx0 + origin_change(0, "in", "m", "em", mdata);
xx0 = frame_change(0, "in", "em", mdata) * xx0;

tf_back_si = 3600 * 24 * 300;  % s, 300 days to see how far they go
tf_back = tf_back_si / mdata.utime_s;
dt_si = 10 * 60;  % s
dt = dt_si / mdata.utime_s;
tspan_back = 0 : -dt : -tf_back;

% Min distance from earth for the apoapsis to count (in case it's circling
% around moon a few times first)
min_r = mdata.ulength_km * 1.2;


%% Find the periapsis altitude at earth departure closest to 180km via grid search

% Earth departe target orbit
target_alt = 180;

N_th0 = 200;  % Increase this potentially
th0_gs_range = linspace(0, 2*pi, N_th0);
h_rps = zeros(N_th0, 1);
wb = waitbar(0, "Grid search (0%)");
for ii = 1 : N_th0
    waitbar(ii / N_th0, wb, "Grid search (" + floor(ii / N_th0 * 100) + "%)");
    h_rps(ii) = full_ballistic_zerof(th0_gs_range(ii), tspan_back, xx0, mdata, sdata);
end
close(wb);
[h_rps_global, th0_global_idx] = min(abs(h_rps - target_alt));
th0_global = th0_gs_range(th0_global_idx);

disp("Found transfer from earth altitude h = " + h_rps_global);
disp("th0 = " + th0_global);

%% Refine grid search with fzero

% Tweak this value a bit maybe
fzero_range = 1/1000;

th0_guess_range = [th0_global - fzero_range, th0_global + fzero_range];
opts = optimset("Display", "iter");
th0_hp_opt = fzero(@(th0) full_ballistic_zerof(th0, tspan_back, xx0, mdata, sdata) - target_alt, th0_guess_range, opts);

%% Plot solution

% If fzero worked
mdata.th0 = th0_hp_opt;
% Otherwise
%mdata.th0 = th0_global;

% Propagate back to earth apoapsis
max_r = mdata.ulength_km * 0.8;
ode_opts = odeset("RelTol", sdata.rel_tol, "AbsTol", sdata.abs_tol, "Events", @(t, xx) event_periapsis_e_maxh(t, xx, max_r, mdata, true));
[t_back, traj_back, t_peri_e] = ode113(@(t, xx) BCRFBP_dyn(t, xx, mdata), tspan_back, xx0, ode_opts);

% Convert trajectory to the one in sun-earth rotating (centred at earth)
traj_back_se = frame_change_traj(t_back, traj_back, "em", "se", mdata);
traj_back_se = origin_change_traj(t_back, traj_back_se, "se", "em", "e", mdata);  % Change origin to earth
traj_back_se = traj_to_si(traj_back_se, mdata);  % Change to km and km/s

% Convert trajectory to earth-centred inertial
traj_back_in_e = frame_change_traj(t_back, traj_back, "em", "in", mdata);
traj_back_in_e = origin_change_traj(t_back, traj_back_in_e, "in", "em", "e", mdata);  % Change origin to earth
traj_back_in_e = traj_to_si(traj_back_in_e, mdata);  % Change to km and km/s

% Plot setup
[Xsph, Ysph, Zsph] = sphere;

R_e_3b = mdata.R_e / mdata.ulength_km;
R_m_3b = mdata.R_m / mdata.ulength_km;

moon_orbit = zeros(100, 3);
for ii = 1 : 100
    moon_orbit(ii, 1) = mdata.ulength_km * cos(2*pi / 100 * ii);
    moon_orbit(ii, 2) = mdata.ulength_km * sin(2*pi / 100 * ii);
end

% Plot in Earth centred sun-earth rotating frame
f_se = figure();
% Plot Trajectory
plot3(traj_back_se(:, 1), traj_back_se(:, 2), traj_back_se(:, 3), "k-", "LineWidth", 1);  % Second part
hold on;
% Plot L1 and L2
plot3(-mdata.L1_se, 0, 0, "k+", "LineWidth", 1, "MarkerSize", 10);  % L1
plot3(mdata.L2_se, 0, 0, "k+", "LineWidth", 1, "MarkerSize", 10);  % L2
% Plot earth
surf(Xsph * mdata.R_e, Ysph * mdata.R_e, Zsph * mdata.R_e, 'FaceColor', '#ffffff');
% Plot moon orbit
plot3(moon_orbit(:, 1), moon_orbit(:, 2), moon_orbit(:, 3), "k--", "LineWidth", 1);
% Adjust plot
axis equal;
grid on;
%title("Sun-Earth Rotating Frame (Earth-Centred)")
xlabel("x / km");
ylabel("y / km");
zlabel("z / km");

% Plot in Earth centred inertial frame
f_in_e = figure();
% Plot Trajectory
plot3(traj_back_in_e(:, 1), traj_back_in_e(:, 2), traj_back_in_e(:, 3), "k-", "LineWidth", 1);
hold on;
% Plot earth
surf(Xsph * mdata.R_e, Ysph * mdata.R_e, Zsph * mdata.R_e, 'FaceColor', '#ffffff');
% Plot moon orbit
plot3(moon_orbit(:, 1), moon_orbit(:, 2), moon_orbit(:, 3), "k--", "LineWidth", 1);
% Adjust plot
axis equal;
grid on;
%title("Earth-Centred Inertial Frame");
xlabel("x / km");
ylabel("y / km");
zlabel("z / km");

% Save some data about the transfer
rp_e = norm(traj_back_in_e(end, 1:3));
ballistic_hp_e = rp_e - mdata.R_e;  % Earth periapsis
dv_rp_e = norm(traj_back_in_e(end, 4:6)) - sqrt(mdata.GM_e / rp_e);
kep_rp_e = car2keps(traj_back_in_e(end, :)', mdata.GM_e);
inc_rp_e_deg = rad2deg(kep_rp_e.inc);
tof_ballistic = -t_peri_e * mdata.utime_s;  % s

disp(" ");
disp("Periapsis alt: " + ballistic_hp_e + "km");
disp("Time of Flight:: " + tof_ballistic + " s (" + tof_ballistic/3600/24/29.530589 + " months)");
disp("Departure dv: " + dv_rp_e + " km/s");
disp("Inclination: " + inc_rp_e_deg + " deg");

%% Functions

function dhp = full_ballistic_zerof(th0, tspan_back, xx0, mdata, sdata)
    mdata.th0 = th0;

    % Propagate back to earth apoapsis and stop at earth periapsis
    max_r = mdata.ulength_km * 0.8;  % Only detect peripasis below this
    ode_opts = odeset("RelTol", sdata.rel_tol, "AbsTol", sdata.abs_tol, "Events", @(t, xx) event_periapsis_e_maxh(t, xx, max_r, mdata, true));
    [~, ~, t_peri_e, xx_peri_e] = ode113(@(t, xx) BCRFBP_dyn(t, xx, mdata), tspan_back, xx0, ode_opts);

    if isempty(t_peri_e)
        dhp = inf;
    else
        dr_e = (xx_peri_e(1, 1:3)' - [-mdata.mu; 0; 0]) * mdata.ulength_km;
        dhp = norm(dr_e) - mdata.R_e;
    end
end