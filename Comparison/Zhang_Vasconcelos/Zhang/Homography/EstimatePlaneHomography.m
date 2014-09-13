%This function computes a 3x3 planar homography H, where pontos_ref are points in the plane and pontos_img are points in the image
%
% function H=EstimatePlaneHomography(pontos_ref,pontos_img)

function H=EstimatePlaneHomography(pontos_ref,pontos_img)

[dummy,n]=size(pontos_img);

%Normalizacao
Kpontos_ref=normalizar(pontos_ref(1:2,:));
pontos_ref=Kpontos_ref*pontos_ref;
Kpontos_img=normalizar(pontos_img(1:2,:));
pontos_img=Kpontos_img*pontos_img;

pontos_ref=pontos_ref*diag(pontos_ref(3,:).^-1);
pontos_img=pontos_img*diag(pontos_img(3,:).^-1);

A=[];
for i=1:n
    A=[A; transpose( kron(pontos_ref(:,i), skew_symetric_v(pontos_img(:,i)) ) )];
end
[U,S,V]=svd(A);
H=[V(1:3,9) V(4:6,9) V(7:9,9)];
H=inv(Kpontos_img)*H*Kpontos_ref;
