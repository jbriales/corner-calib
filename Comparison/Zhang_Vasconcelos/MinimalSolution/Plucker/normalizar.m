function K=normalizar(data);

 [m,n]=size(data);
 Center=data*ones(n,1)*(n^-1);
 sum_isotropic=0;
 for i=1:1:n
  sum_isotropic=sum_isotropic+sqrt(ones(1,m)*((data(:,i)-Center).^2));
 end;
 if sum_isotropic ~=0
  sf=(n*sqrt(2))/sum_isotropic;
 else
  sf=1;
 end;
 K=diag(sf*ones(1,m));
 K=[K -sf*Center; zeros(1,m) 1];
