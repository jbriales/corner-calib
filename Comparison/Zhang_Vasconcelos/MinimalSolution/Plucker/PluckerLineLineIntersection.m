% Estimate the more suitable point of intersection given the Num_lines lines
function Point=PluckerLineLineIntersection(L1,L2)
num=zeros(3,1);
den=zeros(3,3);
L=[L1,L2];
for i=1:2
 if (~isnan(L(1,i)))
  U=L(1:3,i);
  V=L(4:6,i);
  num=num+(cross(V,U)/(transpose(U)*U));
  den=den+(eye(3)-(U*transpose(U)/(transpose(U)*U)));
 end
end

%Point=inv(den)*num;
Point=den\num;