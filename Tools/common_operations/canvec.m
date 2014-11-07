function e = canvec( i, N )
% e = canvec( i, N )
% Return i-th canonical vector of N dimensional space
if nargin==1
    N = 3; % By default 3 dimensional space
end
I = eye(N);
e = I(:,i);
end