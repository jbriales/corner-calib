% function Points = IntersectionLinePlane (Planes, Line)
% 
% This function computes the points of intersection between planes Planes
% and lines Lines (plucker coordinates).


function Points = IntersectionLinePlane (Planes, Lines)

if size(Lines,2)>1
    [m,n]=size(Lines);
    aux=[Planes(4)*eye(3) skew_symetric_v(Planes(1:3)); -Planes(1:3)' zeros(1,3)];
    Points=aux*Lines;
else
    aux=PluckerCoordinates2PluckerMatrix(Lines);
    Points=aux*Planes;
end
% Points=Points(1:3,:)*diag(Points(4,:).^-1);
for i=1:size(Points,2)
    Points(:,i)=Points(:,i)/Points(4,i);
end
Points=Points(1:3,:);