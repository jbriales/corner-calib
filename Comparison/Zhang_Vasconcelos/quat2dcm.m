function [DCM] = quat2dcm(q_in)
% %Compute the direction cosine matrix from Euler Parameters
% % [DCM] = quat2dcm(q_in)
% %
% % Author: Karl Ludwig Fetzer
% 
% DCM(1,1) = 1 - 2 * q_in(2)^2 - 2 * q_in(3)^2;
% DCM(1,2) = 2 * (q_in(1)*q_in(2) - q_in(3)*q_in(4));
% DCM(1,3) = 2 * (q_in(3)*q_in(1) + q_in(2)*q_in(4));
% DCM(2,1) = 2 * (q_in(1)*q_in(2) + q_in(3)*q_in(4));
% DCM(2,2) = 1 - 2 * q_in(3)^2 - 2 * q_in(1)^2;
% DCM(2,3) = 2 * (q_in(2)*q_in(3) - q_in(1)*q_in(4));
% DCM(3,1) = 2 * (q_in(3)*q_in(1) - q_in(2)*q_in(4));
% DCM(3,2) = 2 * (q_in(2)*q_in(3) + q_in(1)*q_in(4));
% DCM(3,3) = 1 - 2 * q_in(1)^2 - 2 * q_in(2)^2;
% 
% end

% MathWorks: http://www.mathworks.es/es/help/aeroblks/quaternionstodirectioncosinematrix.html

q0 = q_in(1);
q1 = q_in(2);
q2 = q_in(3);
q3 = q_in(4);

DCM(1,1) = q0^2+q1^2-q2^2-q3^2;
DCM(1,2) = 2*(q1*q2+q0*q3);
DCM(1,3) = 2*(q1*q3-q0*q2);
DCM(2,1) = 2*(q1*q2-q0*q3);
DCM(2,2) = q0^2-q1^2+q2^2-q3^2;
DCM(2,3) = 2*(q2*q3+q0*q1);
DCM(3,1) = 2*(q1*q3+q0*q2);
DCM(3,2) = 2*(q2*q3-q0*q1);
DCM(3,3) = q0^2+q1^2-q2^2-q3^2;

end