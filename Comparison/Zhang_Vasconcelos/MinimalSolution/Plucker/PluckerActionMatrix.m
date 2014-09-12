function result=PluckerActionMatrix(A);
 
 [P,D]=ComputePluckerMatrices;
 result=inv(D)*P*kron(A,A)*transpose(P);