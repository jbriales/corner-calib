% path = fullfile( pwd, 
path = '/home/jesus/Corner_Datasets/GT2';
scans = loadLidarData('Rawlog', path, '_LASER_LMS');

XY = [scans.xy];
XY = reshape( XY, 2, 361, [] );
xy_gt = mean( XY, 3 );

