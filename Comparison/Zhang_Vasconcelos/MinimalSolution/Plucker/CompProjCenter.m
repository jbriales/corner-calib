% function Centeroid = ComputeLineSetCentroid(PluckerLines)
% 
% This function computes the 3D point that minimizes the sum of euclidean distances to
% a set of N 3D lines.
% The lines are in Plucker coordinates (input matrix 6xN)

function [Center CenterCompPts] = CompProjCenter(Raxels, esp)

[dummy,N]=size(PluckerLines);
U=PluckerLines(1:3,:);
V=PluckerLines(4:6,:);

num=0; den=0;
for i=1:1:N
 num=num+(cross(V(:,i),U(:,i))/(transpose(U(:,i))*U(:,i)));
 den=den+(eye(3)-(U(:,i)*transpose(U(:,i))/(transpose(U(:,i))*U(:,i))));
end;
Centeroid=inv(den)*num;


