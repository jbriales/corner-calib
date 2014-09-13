% 
% this function recives a rotation and translation and returns the rigid body transformation
% If R is 3x3xn t must be 3xn and will return T 4x4xxn
%  

function T=Rt2T(R,t)
n=size(R,3);
if n==1 && ~exist('t','var');
 T=[R [0;0;0];0 0 0 1];
elseif n==size(t,2)
    T=zeros(4,4,n);
    for i=1:n
        T(:,:,i)=[R(:,:,i) t(:,i) ; 0 0 0 1];
    end
else
    fprintf(1,'If R is 3x3xn t must be a 3xn matrix');
end
        