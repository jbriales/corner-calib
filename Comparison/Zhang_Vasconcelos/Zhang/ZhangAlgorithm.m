function  [T,estT ]=ZhangAlgorithm(PI,POINTS)

Nplanes=size(PI,2);
Npts=numel(POINTS);

if Npts~=Nplanes
 fprintf(1,'Must have same number of lines and planes\,');
end


stackNP=[];
stackd=[];
for i=1:Npts
 points=POINTS{i};
 for j=1:size(points,2)
  % they consider that the laser plane is y=0 however we consider that
  % the laser plane is z=0 !!!
  p=[points(1,j);points(2,j);1];
  stackNP=[stackNP ; kron(transpose(p),transpose(PI(:,i))) ];
  stackd=[stackd;transpose(PI(:,i))*PI(:,i)];
 end
end


vecH=pinv(stackNP)*stackd;
H=reshape(vecH,3,3);

estR=[H(:,1), H(:,2) , cross(H(:,1),H(:,2)) ]';
estt=-estR*H(:,3);

[u d v]=svd(estR);
R=u*transpose(v);

t=-R*H(:,3);

estT=Rt2T(estR,estt);
T=Rt2T(R,t);



 