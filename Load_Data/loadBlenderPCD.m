function xy = loadBlenderPCD( file )
% xy = loadBlenderPCD( file )
% Returns a 2xN array with xy coordinates of scanned points
pts = double(loadpcd( file ));
R_s_b = [ 0 0 -1
         -1 0  0
          0 1  0 ];
pts = R_s_b * pts(1:3,:);
xy  = pts(1:2,:);
end