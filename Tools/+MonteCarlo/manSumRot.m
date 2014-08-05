function X = manSumRot( X0, inc_eps )

N = size( inc_eps, 2 );
R = zeros(3,3,N);
R0 = reshape( X0(:), 3, 3 ); % To ensure 3x3 rotation matrix
for i=1:N
    R(:,:,i) = expmap( inc_eps(:,i) ) * R0;
end
X = reshape( R, 9, [] );

end