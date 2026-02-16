function xx_dot = BCRFBP_dyn(t, xx, mdata)
    as = mdata.as;
    ms = mdata.ms;
    mu = mdata.mu;
    th = mdata.th0 + mdata.omega_s * t;

    x = xx(1);
    y = xx(2);
    z = xx(3);
    xdot = xx(4);
    ydot = xx(5);
    zdot = xx(6);

    xx_dot = f_BCRFBP(as,ms,mu,th,x,xdot,y,ydot,z,zdot);
end