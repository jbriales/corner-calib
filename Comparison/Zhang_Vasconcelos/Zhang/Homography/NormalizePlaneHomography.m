
% Following Yi Ma's book (lemma 5.17) a plane homography encoding a rotation and a translation must have the second singular value unitary. This function performs such normalization
%
% function Hnorm=NormalizePlaneHomography(H,pontos_ref,pontos_img)

function Hnorm=NormalizePlaneHomography(H,pontos_ref,pontos_img)

 [U,S,V]=svd(H);
 H=H*S(2,2)^-1;
 if H(3,3)<0;
  H=-H;
 end;
 [dummy,N]=size(pontos_ref);
 aux=find(diag(transpose(pontos_img)*H*pontos_ref)<0);
 [dummy,Naux]=size(aux);
 if isempty(aux) 
  Hnorm=H;
 elseif N==Naux
  Hnorm=-H;
  fprintf('WARNING NormalizePlaneHomography: The Z axis in translation is negative \n');
 else
  Hnorm=H;
  fprintf('WARNING NormalizePlaneHomography: There are points with positive and negative depth \n')
 end;
