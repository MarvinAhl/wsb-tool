% Two body energy, make sure the state is in a planet-centred inertial
% frame in si units
function E = E_2B(xx, GM)
    E = norm(xx(4:6))^2 / 2 - GM / norm(xx(1:3));
end