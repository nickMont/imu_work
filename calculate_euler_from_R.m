function e = calculate_euler_from_R(R)
% Calculates euler from R, see calculate_R_from_euler

e=zeros(3,1);
e(2) = asin(-R(1,3));
e(3) = atan2(R(1,2),R(1,1));
e(1) = atan2(R(2,3),R(3,3));


end

