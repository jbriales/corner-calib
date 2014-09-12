function v = quat2aa(q)

v(4) = 2*acos(q(4));

v(1) = q(1)/sqrt(1-q(4)^2); 
v(2) = q(2)/sqrt(1-q(4)^2);
v(3) = q(3)/sqrt(1-q(4)^2);