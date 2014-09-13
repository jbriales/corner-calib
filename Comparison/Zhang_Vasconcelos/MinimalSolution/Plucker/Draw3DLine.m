%function result=Draw3DLine(L,h,color)
%
% This function draws line L in the 3D axis h (use gca). If h is zero then
% the function creates a 3D plot (cube with 20x20x20 centered in the
% origin). It returns the 6 intersection points of the Line with the 6
% planes of the box. It is not optimized (see slides to improve)

function result=Draw3DLine(L,h,color)

if h==0
 figure;
 axis([-10 10 -10 10 -10 10]);
 grid on;
 h=gca;
end;
aux=axis(h);
BoxBound=[1 1 0 0 0 0;0 0 1 1 0 0; 0 0 0 0 1 1;-aux];
PMat=PluckerCoordinates2PluckerMatrix(L);
Points=PMat*BoxBound;
[U,V,C,d]=GetLineInfo(L);
lambda=[];
dir=[U;0];
c=[C;1];
for i=1:1:6
 a=null([-Points(:,i) c dir]);
 a=a*a(1)^-1;
 if a(2) ~= 0
  lambda=[lambda a(3)/a(2)];
 end;
end;
lambda=sort(lambda);
n=max(size(lambda));
%This can be improved ... just using the maximum and minimum lambda
P1=C+lambda(1)*U;
P2=C+lambda(n)*U;
aux=[transpose(P1);transpose(P2)];
line(aux(:,1),aux(:,2),aux(:,3),'Color',color);
result=Points;
 

