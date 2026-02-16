% Computes frame change matrix from frame1 to frame2 for states that have
% earth-moon barycentre as origin. Available frames are "em" (earth moon),
% "se" (sun-earth), "in" (inertial, note: not truly inertial because still
% moving with earth-moon barycentre). t is unitless

function M12 = frame_change(t, frame1, frame2, mdata)
    % Always converts from frame1 to em and then from em to frame2

    % in se frame, sun is always on the left!

    % Angles of the frames
    phi = mdata.phi0 + t;
    th = mdata.th0 + mdata.omega_s * t;

    % Position rotation matrices
    R_em_in = [cos(phi), -sin(phi), 0;
               sin(phi), cos(phi), 0;
               0, 0, 1];
    R_em_se = [cos(th), -sin(th), 0;
               sin(th), cos(th), 0;
               0, 0, 1];

    % Some crossproduct helper matrices
    w_em_cr = [0, -1, 0;
               1, 0, 0;
               0, 0, 0];
    w_s_cr = [0, -mdata.omega_s, 0;
              mdata.omega_s, 0, 0;
              0, 0, 0];

    % Initial conv matrices
    % in to em
    M_in_em = [R_em_in', zeros(3);
               -w_em_cr*R_em_in', R_em_in'];
    % se to em
    M_se_em = [R_em_se', zeros(3);
               -w_s_cr*R_em_se', R_em_se'];

    % Final conv matrices
    % em to in
    M_em_in = [R_em_in, zeros(3);
               R_em_in*w_em_cr, R_em_in];
    % em to se
    M_em_se = [R_em_se, zeros(3);
               R_em_se*w_s_cr, R_em_se];

    % Assemble transformation matrix
    M12 = eye(6);

    if strcmp(frame1, "in")
        M12 = M_in_em * M12;
    elseif strcmp(frame1, "se")
        M12 = M_se_em * M12;
    end

    if strcmp(frame2, "in")
        M12 = M_em_in * M12;
    elseif strcmp(frame2, "se")
        M12 = M_em_se * M12;
    end
end