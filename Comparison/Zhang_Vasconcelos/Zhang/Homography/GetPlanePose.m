%The function computes the pose of the plane with respect to the camera given the correspondences 
%between plane and canonical image plane (the intrinsics are compensated in advance)
function [Tp2i,dif] =GetPlanePose(pontos_ref, pontos_img)

% Estimate the homography using svd
Hi=EstimatePlaneHomography(pontos_ref,pontos_img);
% Normalize the results
H=NormalizePlaneHomography(Hi,pontos_ref,pontos_img);
% Compute the 4 solutions for the pose
[r,N,T]=FactorizePlaneHomography(H);

% Pick the solutions where the plane is ahead of the camera for which N=[0;0;1]
if N(3,1)>0
    N=N(:,[1 4]);
    T=T(:,[1 4]);
    R(:,:,1)=r(:,:,1);
    R(:,:,2)=r(:,:,4);
else
    N=N(:,[2 3]);
    T=T(:,[2 3]);
    R(:,:,1)=r(:,:,2);
    R(:,:,2)=r(:,:,3);
end;

% The closest to N=[0;0;1]
dif=sqrt(ones(1,3)*(N-[0 0;0 0; 1 1]).^2);
[dummy,i]=min(dif);
% Move the origin toward the plane
t=R(:,:,i)*[0;0;1]+T(:,i);
Tp2i=[R(:,:,i) t;0 0 0 1];

