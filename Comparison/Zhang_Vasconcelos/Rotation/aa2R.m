% axis-angle to rotation matrix (rodrigues formula)

function R = aa2R(r)

theta = sqrt(r(1)^2 + r(2)^2 + r(3)^2);
r     = r/theta;
kr    = kron(r,r.');
R     = kr + cos(theta)*(eye(3) - kr) + sin(theta)*skew_symetric_v(r);