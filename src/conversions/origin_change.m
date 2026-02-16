% Possible frames: "em", "se", "in"
% Possible origins: "e" earth, "m" moon, "em" earth-moon barycentre

function xx_shift = origin_change(t, frame, origin1, origin2, mdata)
    % Move origin from origin1 to em barycentre and then from there to
    % origin 2

    % Angles of the frames
    % phi = mdata.phi0 + t;
    % th = mdata.th0 + mdata.omega_s * t;
    % 
    % % Position rotation matrices
    % R_em_in = [cos(phi), -sin(phi), 0;
    %            sin(phi), cos(phi), 0;
    %            0, 0, 1];
    % R_em_se = [cos(th), -sin(th), 0;
    %            sin(th), cos(th), 0;
    %            0, 0, 1];

    % Positions of bodies with em in origin
    xx_earth = [-mdata.mu; 0; 0; 0; 0; 0];
    xx_moon = [1-mdata.mu; 0; 0; 0; 0; 0];

    if strcmp(frame, "in")
        xx_earth = frame_change(t, "em", "in", mdata) * xx_earth;
        xx_moon = frame_change(t, "em", "in", mdata) * xx_moon;
    elseif strcmp(frame, "se")
        xx_earth = frame_change(t, "em", "se", mdata) * xx_earth;
        xx_moon = frame_change(t, "em", "se", mdata) * xx_moon;
    end

    % Add this shift vector to the position to make it work
    xx_shift = zeros(6, 1);
    if strcmp(origin1, "e")
        xx_shift = xx_shift + xx_earth;
    elseif strcmp(origin1, "m")
        xx_shift = xx_shift + xx_moon;
    end

    if strcmp(origin2, "e")
        xx_shift = xx_shift - xx_earth;
    elseif strcmp(origin2, "m")
        xx_shift = xx_shift - xx_moon;
    end
end