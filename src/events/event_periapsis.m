function [value, isterminal, direction] = event_periapsis(t, xx, body, mdata, terminal)
    % Transform to inertial frame around body where periapsis is supposed
    % ot happen
    xx_in = frame_change(t, "em", "in", mdata) * xx;
    xx_in = xx_in + origin_change(t, "in", "em", body, mdata);
    xx_in = state_to_si(xx_in, mdata);

    rr = xx_in(1:3);
    vv = xx_in(4:6);

    value = dot(rr, vv);  % Periapsis has position and velocity parallel
    isterminal = terminal;  % If to stop if periapsis is reached
    direction = 1;  % At peri dot(rr, vv) switches from neg to pos (increasing value)
end