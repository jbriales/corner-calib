function [c, l] = uncalibrateData( c,l, K )
% % [N, c, A_co, hL] = pts2calibratedData( P, K )
% % Input:
% %   x - 1x5 array with: c point (in pixels), 3 angles of direction
% %   K - intrinsic calibration matrix 
% % Output:
% %   N  - 2x3 array with direction of lines
% %   c  - 2D corner point
% %   A_co - covariance matrix with uncertainty of output data
% %   L_P2 - 3x3 array whose columns are (hom.) calibrated lines

for k=1:3
    l{k} = K' \ l{k};
    l{k} = l{k} / norm(l{k}(1:2));
end
c = makeinhomogeneous( K * makehomogeneous(c) );

end