function out=PluckerDual(M)

[m n k]=size(M);
if isequal([m n],[4 4])
 out=zeros(4,4,k);
 for i=1:k
  out(1,2,i)=M(3,4,i);
  out(1,3,i)=M(4,2,i);
  out(1,4,i)=M(2,3,i);
  out(2,3,i)=M(1,4,i);
  out(4,2,i)=M(1,3,i);
  out(3,4,i)=M(1,2,i);
  out(:,:,i)=out(:,:,i)-out(:,:,i)';
 end
elseif m==6
%  out=zeros(m,n);
%  for i=1:n
%   out(:,i)=[zeros(3) eye(3);eye(3) zeros(3)]*M(:,i);
%  end
 out=M([4:6 1:3],:);
else
 error('Input must be Plucker matrix (4x4xn) or Plucker coordinates (6xn)');
end