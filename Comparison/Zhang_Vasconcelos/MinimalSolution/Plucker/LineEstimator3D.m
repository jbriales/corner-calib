% function L=LineEstimator3D(Q)
%
% This function estimates the line L (Plucker coordinates) going through
% the set of points Q=[Q1 Q2 ... Qn]. Q is 3xn

function L=LineEstimator3D(Q)

 ERR=10^-12;
 [m,n]=size(Q);
 K=normalizar(Q);
 %Note: K is a similarity transformation which means that it does not
 Qn=K*[Q;ones(1,n)];
 %Colinear points (do not bother)
 if rank(Qn)<3
  Ln=PluckerActionVector(Qn(:,1),Qn(:,2));
 else
  A=[];
  for i=1:1:n
   A=[A; skew_symetric_v(Qn(1:3,i)) eye(3,3)];
  end;
  [U,S,V]=svd(A);
  Ln=V(:,6);
  %I do not have a formal proof, but the result is already in the Klein
  %quadric. This is a safety code.
  if abs(transpose(Ln(1:3))*Ln(4:6))>ERR
   %Plucker correction following Adrien's paper
   a=Ln(1:3);
   b=Ln(4:6);
   [U,S,V]=svd([a b],0);
   Z=S*transpose(V);
   T=[Z(2,1) Z(2,2);Z(1,2) -Z(1,1)];
   [Ut,St,Vt]=svd(T);
   Vh=[Vt(1,2) -Vt(2,2); Vt(2,2) Vt(1,2)];
   Ln=reshape(U*Vh*diag(diag(transpose(Vh)*S*transpose(V))),6,1);
  end;
 end;
 %Take outA normalization
 L=inv(PluckerActionMatrix(K))*Ln;
 