%function r=DistancePointLine(P,L)
%
%This function computes the orthogonal distance between points P and line L
%(plucker coordinates). P=P1 P2 ... Pn] is 3xn.


function r=DistancePointLine(P,L)
[m,n]=size(P); 
[U,V,C,d]=GetLineInfo(L);
r=[];
 for i=1:1:n
  v=P(:,i)-C;
  aux=C+(transpose(v)*U)*U;
  r=[r norm(P(:,i)-aux)];
 end;
 