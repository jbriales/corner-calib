function result=PluckerActionVector(A,B);

 m=max(size(A));
 %Put in homogeneous coordinates
 if m==3
  A=[A;1];
  B=[B;1];
 end;
 aux=reshape(A*transpose(B)-B*transpose(A),16,1);
 [P,D]=ComputePluckerMatrices;
 result=inv(D)*P*aux;