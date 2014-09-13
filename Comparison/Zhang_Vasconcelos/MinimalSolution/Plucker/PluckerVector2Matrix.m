
function out=PluckerVector2Matrix(in)
[m n k]=size(in);
if isequal([m n],[4 4])
 out=zeros(6,k);
 for i=1:k
  out(:,i)=[in(1:3,4,i);in(2,3,i);-in(1,3,i);in(1,2,i)];
 end
elseif m==6
 out=zeros(4,4,n);
 for i=1:n
  out(:,:,i)=[-skew_symetric_v(in(4:6,i)) in(1:3,i);-in(1:3,i)' 0];
 end
else
 error('Must be a Plucker vector (6x1) or a Plucker matrix (4x4) ');
end
 
