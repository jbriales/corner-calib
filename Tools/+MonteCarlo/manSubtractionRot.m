function inc_eps = manSubtractionRot( X, X0 )

R  = reshape( X, 3,3, [] );
R0 = reshape( X0, 3,3,1 );
N = size(X,2);
inc_eps = zeros(3,N);
for i=1:N
    inc_eps(:,i) = logmap( R(:,:,i) * R0' );
end

end

function u = logmap(R)
%LOGMAP The Rodrigues' formula for rotation matrices.
[~,~,V] = svd(R - eye(3));
u = V(:,3);

c2 = trace(R) - 1; % cos(theta)*2
s2 = u' * [R(3,2)-R(2,3) , R(1,3)-R(3,1) , R(2,1)-R(1,2)]';
theta = atan2(s2,c2);

u = u * theta;
end