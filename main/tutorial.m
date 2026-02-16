clear variables
close all
clc

% Add src and all subfolders to path
addpath(genpath("../src"));

mdata = init_model_data(true);
sdata = init_sim_data(mdata);

% Initial state LEO
% xx0_inert = [6500; 0; 0; 0; sqrt(3.986004418e5/6500); 0];
% xx0 = state_from_si(xx0_inert, mdata);
% xx0 = xx0 + origin_change(0, "in", "e", "em", mdata);
% xx0 = frame_change(0, "in", "em", mdata) * xx0;

% Initial state LLO
% xx0_inert = [1837; 0; 0; 0; sqrt(4.902800118e3/1837); 0];
% xx0 = state_from_si(xx0_inert, mdata);
% xx0 = xx0 + origin_change(0, "in", "m", "em", mdata);
% xx0 = frame_change(0, "in", "em", mdata) * xx0;

% Sun-Earth L2
% xx0_inert = [2.165e6; 0; 0; 0; 0; 0];
% xx0 = state_from_si(xx0_inert, mdata);
% xx0 = frame_change(0, "se", "em", mdata) * xx0;

% Initial state to something looking like a low energy transfer
rp = 6500;
ra = 0.7785e6;
let_kep.a = (ra + rp)/2;
let_kep.e = 1 - rp/let_kep.a;
let_kep.inc = 0;
let_kep.raan = pi - pi/16;
let_kep.aop = 0;
let_kep.ta = 0;
xx0_inert = kep2cars(let_kep, mdata.GM_e);
xx0 = state_from_si(xx0_inert, mdata);
xx0 = xx0 + origin_change(0, "in", "e", "em", mdata);
xx0 = frame_change(0, "in", "em", mdata) * xx0;

% Integrate
tf_si = 3600 * 24 * 85;  % s
tspan = 0 : sdata.dt : tf_si / mdata.utime_s;
ode_opts = odeset("RelTol", sdata.rel_tol, "AbsTol", sdata.abs_tol);
[~, traj] = ode113(@(t, xx) BCRFBP_dyn(t, xx, mdata), tspan, xx0, ode_opts);

% Convert trajectory to inertial one around earth
traj_in_e = frame_change_traj(tspan, traj, "em", "in", mdata);
traj_in_e = origin_change_traj(tspan, traj_in_e, "in", "em", "e", mdata);  % Change origin to earth
traj_in_e = traj_to_si(traj_in_e, mdata);  % Change to km and km/s

% Convert trajectory to inertial one around moon
traj_in_m = frame_change_traj(tspan, traj, "em", "in", mdata);
traj_in_m = origin_change_traj(tspan, traj_in_m, "in", "em", "m", mdata);  % Change origin to moon
traj_in_m = traj_to_si(traj_in_m, mdata);  % Change to km and km/s

% Convert trajectory to the one in sun-earth rotating (centred at earth)
traj_se = frame_change_traj(tspan, traj, "em", "se", mdata);
traj_se = origin_change_traj(tspan, traj_se, "se", "em", "e", mdata);  % Change origin to earth
traj_se = traj_to_si(traj_se, mdata);  % Change to km and km/s

% Compute Jacobi constant and Earth and Moon 2 body energy
J = zeros(length(tspan), 1);
E2b_e = zeros(length(tspan), 1);
E2b_m = zeros(length(tspan), 1);
for ii = 1 : length(tspan)
    xx = traj(ii, :)';
    J(ii) = J_CRTBP(mdata.mu, xx(1), xx(4), xx(2), xx(5), xx(3), xx(6));

    % Change frame to inertial
    xx_in = frame_change(tspan(ii), "em", "in", mdata) * xx;
    % Move origin to earth and moon
    xx_in_e = xx_in + origin_change(tspan(ii), "in", "em", "e", mdata);
    xx_in_m = xx_in + origin_change(tspan(ii), "in", "em", "m", mdata);
    % Transform them to SI units
    xx_in_e = state_to_si(xx_in_e, mdata);
    xx_in_m = state_to_si(xx_in_m, mdata);

    E2b_e(ii) = E_2B(xx_in_e, mdata.GM_e);
    E2b_m(ii) = E_2B(xx_in_m, mdata.GM_m);
end

tspan_si = tspan * mdata.utime_s;

%% Plots
[Xsph, Ysph, Zsph] = sphere;

R_e_3b = mdata.R_e / mdata.ulength_km;
R_m_3b = mdata.R_m / mdata.ulength_km;

moon_orbit = zeros(100, 3);
for ii = 1 : 100
    moon_orbit(ii, 1) = mdata.ulength_km * cos(2*pi / 100 * ii);
    moon_orbit(ii, 2) = mdata.ulength_km * sin(2*pi / 100 * ii);
end

% Plot in 3B frame
figure();
hold on;
plot3(traj(:, 1), traj(:, 2), traj(:, 3), "k-", "LineWidth", 1);
% Plot earth and moon
surf(Xsph * R_e_3b - mdata.mu, Ysph * R_e_3b, Zsph * R_e_3b, 'FaceColor', '#ffffff');
surf(Xsph * R_m_3b + 1 - mdata.mu, Ysph * R_m_3b, Zsph * R_m_3b, 'FaceColor', '#ffffff');
% Adjust plot
xlim([min([-1; traj(:, 1)]), max([1; traj(:, 1)])] * 1.2)
ylim([min([-1; traj(:, 2)]), max([1; traj(:, 2)])] * 1.2)
axis equal;
grid on;
title("Earth-Moon Rotating Frame");
xlabel("x");
ylabel("y");

% Plot in Earth centred inertial frame
figure();
hold on;
plot3(traj_in_e(:, 1), traj_in_e(:, 2), traj_in_e(:, 3), "k-", "LineWidth", 1);
% Plot earth
surf(Xsph * mdata.R_e, Ysph * mdata.R_e, Zsph * mdata.R_e, 'FaceColor', '#ffffff');
% Plot moon orbit
plot3(moon_orbit(:, 1), moon_orbit(:, 2), moon_orbit(:, 3), "k--", "LineWidth", 1);
% Adjust plot
xlim([min([-mdata.ulength_km; traj_in_e(:, 1)]), max([mdata.ulength_km; traj_in_e(:, 1)])] * 1.2)
ylim([min([-mdata.ulength_km; traj_in_e(:, 2)]), max([mdata.ulength_km; traj_in_e(:, 2)])] * 1.2)
axis equal;
grid on;
title("Earth-Centred Inertial Frame");
xlabel("x / km");
ylabel("y / km");

% Plot in Moon centred inertial frame
figure();
hold on;
plot3(traj_in_m(:, 1), traj_in_m(:, 2), traj_in_m(:, 3), "k-", "LineWidth", 1);
% Plot moon
surf(Xsph * mdata.R_m, Ysph * mdata.R_m, Zsph * mdata.R_m, 'FaceColor', '#ffffff');
% Adjust plot
xlim([-mdata.R_m, mdata.R_m] * 5);
axis equal;
grid on;
title("Moon-Centred Inertial Frame");
xlabel("x / km");
ylabel("y / km");

% Plot in Earth centred sun-earth rotating frame
figure();
hold on;
plot3(traj_se(:, 1), traj_se(:, 2), traj_se(:, 3), "k-", "LineWidth", 1);
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
xlim([min([-mdata.L1_se; traj_se(:, 1)]), max([mdata.L2_se; traj_se(:, 1)])] * 1.1);
%ylim([min(traj_se(:, 2)), max(traj_se(:, 2))] * 1.2);
title("Sun-Earth Rotating Frame (Earth-Centred)")
xlabel("x / km");
ylabel("y / km");

% Plot Jacobi constant
figure();
hold on;
plot(tspan_si, J, "k-", "LineWidth", 1);
grid on;
title("Jacobi Energy")
xlabel("t / s");
ylabel("J");

% Plot Earth 2 body energy
figure();
hold on;
plot(tspan_si, E2b_e, "k-", "LineWidth", 1);
grid on;
title("Earth 2-Body Energy")
xlabel("t / s");
ylabel("E_e / km^2/s^2");

% Plot Moon 2 body energy
figure();
hold on;
plot(tspan_si, E2b_m, "k-", "LineWidth", 1);
grid on;
title("Moon 2-Body Energy")
xlabel("t / s");
ylabel("E_m / km^2/s^2");

%% Appendix: Spice part
init_spice;

%% Transform to earth-centred J2000 or earth_fixed

% this et0 doesn't correspond the the moon and sun angle orientation in
% this example!! It's only for demonstrating the conversion. To make moon
% and sun angles accurate use init_model_data_epoch in the beginning
et0 = string_to_epoch("2025-11-02 18:19:00.000 TDB");

traj_J2000 = traj_ecl_to_equ(et0 + tspan_si, traj_in_e, false);
traj_IAUEARTH = traj_ecl_to_equ(et0 + tspan_si, traj_in_e, true);

% Plot them both

% Plot in Earth centred J2000 frame
figure();
hold on;
plot3(traj_J2000(:, 1), traj_J2000(:, 2), traj_J2000(:, 3), "k-", "LineWidth", 1);
% Plot earth
surf(Xsph * mdata.R_e, Ysph * mdata.R_e, Zsph * mdata.R_e, 'FaceColor', '#ffffff');
% Adjust plot
xlim([min([-mdata.ulength_km; traj_J2000(:, 1)]), max([mdata.ulength_km; traj_J2000(:, 1)])] * 1.2)
ylim([min([-mdata.ulength_km; traj_J2000(:, 2)]), max([mdata.ulength_km; traj_J2000(:, 2)])] * 1.2)
axis equal;
grid on;
title("Earth-Centred J2000 Frame");
xlabel("x / km");
ylabel("y / km");

% Plot in Earth centred earth fixed frame
figure();
hold on;
plot3(traj_IAUEARTH(:, 1), traj_IAUEARTH(:, 2), traj_IAUEARTH(:, 3), "k-", "LineWidth", 1);
% Plot earth
surf(Xsph * mdata.R_e, Ysph * mdata.R_e, Zsph * mdata.R_e, 'FaceColor', '#ffffff');
% Adjust plot
xlim([min([-mdata.ulength_km; traj_IAUEARTH(:, 1)]), max([mdata.ulength_km; traj_IAUEARTH(:, 1)])] * 1.2)
ylim([min([-mdata.ulength_km; traj_IAUEARTH(:, 2)]), max([mdata.ulength_km; traj_IAUEARTH(:, 2)])] * 1.2)
axis equal;
grid on;
title("Earth-Centred Earth-Fixed Frame");
xlabel("x / km");
ylabel("y / km");

%% Get Azimuth and Elevation data of a ground station for the trajectory

% Some ground station codes
gs_codes = ground_stations;

% Transform trajectory to local one
[traj_rho, traj_az, traj_el] = traj_to_topo_gs(et0 + tspan_si, traj_in_e, gs_codes.ESA_New_Norcia);

% Plot elevation
figure();
hold on;
plot(tspan_si, rad2deg(traj_el), "k-", "LineWidth", 1);
yline(6, "r--", "LineWidth", 1);  % Assume 6deg min elevation angle
grid on;
title("Elevation angle of ESTRACK New Norcia Station");
xlabel("t / s");
ylabel("el / deg");

%% Unload all Spice Kernels
close_spice;