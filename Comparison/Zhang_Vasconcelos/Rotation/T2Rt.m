% 
% this function recives a transformation matix and returns the rotation and
% translation
% If T is 4x4xn will return R 3x3xn and t 3xn
%  

function [R,t]=T2Rt(T)
n=size(T,1);
if n==size(T,2)
    k=size(T,3);
    R=zeros(3,3,k);
    t=zeros(3,k);
    for i=1:k
        R(:,:,i)=T(1:3,1:3,i);
        t(:,i)=T(1:3,4,i);
    end
else
    fprintf(1,'T must be a 4x4xn matrix');
    R=[];
    t=[];
end
        