function [value, isterminal, direction] = event_escape(t, xx, body, GM, r_soi, mdata, terminal)
    % Transform to inertial frame around body where periapsis is supposed
    % ot happen
    xx_in = frame_change(t, "em", "in", mdata) * xx;
    xx_in = xx_in + origin_change(t, "in", "em", body, mdata);
    xx_in = state_to_si(xx_in, mdata);

    rr = xx_in(1:3);

    E2B = E_2B(xx_in, GM);  % Positive if on a non-captured 2b orbit
    d_soi = norm(rr) - r_soi;  % Positive if outside SOI

    value = (E2B > 0 && d_soi > 0) * 2 - 1;  % Is only >= 0 if E2B and d_soi are both positive
    isterminal = terminal;  % If to stop if periapsis is reached
    direction = 1;  % Escapes if switch from neg to pos (increasing value)
end