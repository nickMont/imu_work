function R = R_pqr(e)
% Generates R_pqr for a 3-2-1 rotation such that
% omegaB = R*[phidot;thetadot;psidot]

R=[1 0 -sin(e(2))
    0 cos(e(1)) sin(e(1))*cos(e(2))
    0 -sin(e(1)) cos(e(1))*cos(e(2))];

end

