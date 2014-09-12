function [R,N,T]=FactorizePlaneHomography(H)

 [V,S]=svd(transpose(H)*H);
 v1=V(:,1); v2=V(:,2); v3=V(:,3);
% v1=v1*sqrt(transpose(v1)*v1)^-1;
% v2=v2*sqrt(transpose(v2)*v2)^-1;
% v3=v3*sqrt(transpose(v3)*v3)^-1;

 s=diag(S);
 u1=(sqrt(1-s(3))*v1+sqrt(s(1)-1)*v3)/sqrt(s(1)-s(3));
 u2=(sqrt(1-s(3))*v1-sqrt(s(1)-1)*v3)/sqrt(s(1)-s(3));

 U1=[v2 u1 skew_symetric_v(v2)*u1];
 U2=[v2 u2 skew_symetric_v(v2)*u2];

 W1=[H*v2 H*u1 skew_symetric_v(H*v2)*H*u1];
 W2=[H*v2 H*u2 skew_symetric_v(H*v2)*H*u2];

 R(:,:,1)=W1*transpose(U1); N(:,1)=skew_symetric_v(v2)*u1; T(:,1)=(H-R(:,:,1))*N(:,1);
 R(:,:,2)=W2*transpose(U2); N(:,2)=skew_symetric_v(v2)*u2; T(:,2)=(H-R(:,:,2))*N(:,2);
 R(:,:,3)=R(:,:,1); N(:,3)=-N(:,1); T(:,3)=-T(:,1);
 R(:,:,4)=R(:,:,2); N(:,4)=-N(:,2); T(:,4)=-T(:,2);

