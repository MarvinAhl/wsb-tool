function model_data = init_model_data(is_4bp)
    % Non-dim units for model
    model_data.mu = 1.2150e-2;
    model_data.ms = 0;
    if is_4bp
        model_data.ms = 3.2890e5;  % Set this to 0 to just get CRTBP dynamics
    end
    model_data.as = 3.8881e2;
    model_data.omega_s = 0.9251;  % Ang vel is clockwise!!!

    % Initial parameters for relative angles
    model_data.phi0 = 0;  % Initial earth-moon-phase w.r.t. inertial frame (couter-clockwise)
    model_data.th0 = 0;  % Initial sun-phase (clockwise from earth-moon -x axis)

    % SI units for conversions
    model_data.ulength_km = 384399;  % km
    model_data.umass_kg = 5.9722e24 + 7.346e22;  % kg
    model_data.uangfreq_sinv = 2*pi / (27.321661 * 24 * 60 * 60);  % 1/s
    % Not a real unit value but still in here because so commonly used
    model_data.utime_s = 1 / model_data.uangfreq_sinv;  % s

    % Some physical parameters
    model_data.GM_e = 3.986004418e5;  % km^3/s^2 Standard gravitational parameter earth
    model_data.GM_s = 1.32712440042e11;  % km^3/s^2 Standard gravitational parameter moon
    model_data.GM_m = 4.902800118e3;  % km^3/s^2 Standard gravitational parameter sun
    model_data.R_e = 6371;  % km
    model_data.R_m = 1737.5;  % km
    model_data.L1_se = 1.4975e6;  % km from earth
    model_data.L2_se = 1.5079e6;  % km from earth
    model_data.r_soi_m = model_data.ulength_km * (model_data.GM_m/model_data.GM_e)^(2/5);
end